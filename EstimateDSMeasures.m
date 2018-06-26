%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates SNR,CR,HET,DEN,NUM,BLE from centroids                         %
% A voxel is considered FG if distance from the closest voxel  < TH_FG    %
% A voxel is considered BG if distance from the closest voxel  > TH_BG    %
% For each time point:                                                    %
% - SNR = |avg(FG) - avg(BG)| / |std(BG)|                                 %
% - CR = avg(FG) / avg(BG)                                                %
% - HET = std(FG) / |avg(FG) - avg(BG)|                                   %
% - DEN = minimum distance between two spots                              %
% - NUM = number of spots                                                 %
% Total:                                                                  %
% - mean and std (over time) of SNR, CR, HET, DEN, NUM                    %
% Notes:                                                                  %
% - Requires imaging data as Imaris files and tracks in the LTDB format   %
% - SNR, CR, HET are estimated according with [1]                         %
% - Distance transform is computed using [2]                              %
% - Imaris files are interpreted using [3]                                %
%                                                                         %
% 1. Ulman et al, Nature Methods (2017)                                   %
% 2. https://ch.mathworks.com/matlabcentral/fileexchange/15455-3d-        %
%    euclidean-distance-transform-for-variable-data-aspect-ratio          %
% 3. https://github.com/PeterBeemiller/ImarisReader                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

%% Settings - Edit here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LTDB_PATH_IMS = 'G:\LTDB_r2\GT_IMS\';           %Imaris files
LTDB_PATH_TRACKS = 'G:\LTDB_r2\GT_TRACKS\';     %LTDB tracks
LTDB_PATH_MEASURES = 'G:\LTDB_r2\DS_MEASURES\'; %Output path

SAMPLING_DISTANCE = 1;  %[um] Regular grid distance for downsampling
TH_DISTANCE_BG = 20;    %TH_BG [um] Multiple of SAMPLING_DISTANCE
TH_DISTANCE_FG = 4;     %TH_FG [um] Multiple of SAMPLING_DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Included libraries
addpath '.\libs\ImarisReader-master';
addpath '.\libs\bwdistsc';

%Processing all the files in the folder
files = dir([LTDB_PATH_TRACKS, '*.xls']);

NUM_ENTRIES = size(files, 1);
overall = zeros(NUM_ENTRIES,21); % Contains the final average values
overall_names = cell(NUM_ENTRIES,1);
overall_tranames = cell(NUM_ENTRIES,1);

ifiles = 0;

% For each file
for file = files'
    ifiles = ifiles+1;
    [~,fname,fext] = fileparts(file.name);
    disp(['Reading ', fname]);
    fns = split(fname, '_');
    fn_ims = [fns{1}, '_GT.ims'];
    % Read the corresponding Imaris file
    try
    curr_ims = ImarisReader([LTDB_PATH_IMS, fn_ims]);
    catch
        disp(['Error initializing in ', fn_ims]);
        continue;
    end
    
    % Get dataset properties
    xmin = curr_ims.DataSet.ExtendMinX;
    ymin = curr_ims.DataSet.ExtendMinY;
    zmin = curr_ims.DataSet.ExtendMinZ;
    
    xmax = curr_ims.DataSet.ExtendMaxX;
    ymax = curr_ims.DataSet.ExtendMaxY;
    zmax = curr_ims.DataSet.ExtendMaxZ;
    
    W = curr_ims.DataSet.SizeX;
    H = curr_ims.DataSet.SizeY;
    D = curr_ims.DataSet.SizeZ;    
    C = curr_ims.DataSet.SizeC;    
    T = curr_ims.DataSet.SizeT;    
    
    [X,Y,Z] = meshgrid(1:W,1:H,1:D);
    
    vx = (xmax-xmin)/W;
    vy = (ymax-ymin)/H;
    vz = (zmax-zmin)/D;
    dt = seconds(datetime(curr_ims.DataSet.Timestamps{2}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS') - datetime(curr_ims.DataSet.Timestamps{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS'));
    
    if(vx ~= vy)
        disp(['--- ', file.name, ' ERROR: vx=', num2str(vx), ' vy=', num2str(vy), ' vz=', num2str(vz), ' USING vy = vx = mean(vx,vy)']);
        vy = mean([vx, vy]);
        vx = vy;
    end
    
    if((xmin > 0) || (ymin > 0) || (zmin > 0))
        disp(['--- ', file.name, ' WARNING: xmin=', num2str(xmin), ' ymin=', num2str(ymin), ' zmin=', num2str(zmin)]);
    end
    
    % Loac the entire dataset into RAM
    tot_td = 0;
    num_tracks = 0;
    try
        zstack = curr_ims.DataSet.GetData(); %WHDCTs
    catch
        disp(['ERROR reading ', fn_ims]);
        continue
    end
    
    % Read the tracks
    spots_IXYZT = xlsread([LTDB_PATH_TRACKS,file.name]);
    spots_IXYZT_um = spots_IXYZT;
    
    spots_IXYZT(3:end, 2) = (spots_IXYZT(3:end, 2))./spots_IXYZT(1, 2);
    spots_IXYZT(3:end, 3) = (spots_IXYZT(3:end, 3))./spots_IXYZT(1, 3);
    spots_IXYZT(3:end, 4) = spots_IXYZT(3:end, 4)./spots_IXYZT(1, 4);
    
    spots_IXYZT = spots_IXYZT(3:end,:);
    spots_IXYZT_um = spots_IXYZT_um(3:end,:);
    
    traxmin = min(spots_IXYZT(:,2));
    traymin = min(spots_IXYZT(:,3));
    trazmin = min(spots_IXYZT(:,4));
    disp(['--- track minimum: ', num2str([traxmin, traymin, trazmin])]);
    
    % Compute the measures
    SNR_T = zeros(T, C)+NaN;
    CR_T = zeros(T, C)+NaN;
    HET_T = zeros(T, C)+NaN;
    DEN_T = zeros(T,1)+NaN;
    NUM_T = zeros(T,1)+NaN;
    FG_T = zeros(T, C)+NaN;
    
    dsx = floor(SAMPLING_DISTANCE/vx)+1; %downsampling factor along x
    dsy = floor(SAMPLING_DISTANCE/vy)+1; %downsampling factor along x
    dsz = floor(SAMPLING_DISTANCE/vz)+1; %downsampling factor along x

    [x_sampling, y_sampling, z_sampling] = meshgrid(1:dsx:W, 1:dsy:H, 1:dsz:D);
    x_sampling = x_sampling(:);
    y_sampling = y_sampling(:);
    z_sampling = z_sampling(:);

    h = waitbar(0, ['Computing DS measures of ', file.name]);
    for tf = 1:T
        waitbar(tf/T, h);
        spots_idx = find(spots_IXYZT(:,5) == tf);
        if(isempty(spots_idx))
            disp(['No spots in tf ', num2str(tf)]);
            continue;
        end
        
        spots_spots_dist = pdist2(spots_IXYZT_um(spots_idx,2:4), spots_IXYZT_um(spots_idx,2:4));
        d = logical(eye(size(spots_spots_dist, 1)));
        spots_spots_dist(d) = 9999;
        DEN = min(spots_spots_dist(:));
        if((DEN == 9999) || (DEN == 0))
            DEN = NaN;
        end
        NUM = size(spots_spots_dist,1);   
        DEN_T(tf) = DEN;
        NUM_T(tf) = NUM;
        
        %% Obtain spots values
        x_sp = floor(spots_IXYZT(spots_idx,2))+1;
        y_sp = floor(spots_IXYZT(spots_idx,3))+1;
        z_sp = floor(spots_IXYZT(spots_idx,4))+1;   
        l_sp = floor(spots_IXYZT(spots_idx,1));      
        
        %% Compute the 3d euclidean distance transformation from spots
        BW = false(W,H,D);
        for ii = 1:numel(spots_idx)
            if ((x_sp(ii) < 1) || (x_sp(ii) > W) || (y_sp(ii) < 1) || (y_sp(ii) > H) || (z_sp(ii) < 1) || (z_sp(ii) > D))
                disp(['ERROR: position: ', num2str([x_sp(ii), y_sp(ii), z_sp(ii)]), ' T: ', num2str(tf)]);
            end
            if(x_sp(ii) < 1)
                x_sp(ii) = 1;
            end
            
            if(y_sp(ii) < 1)
                x_sp(ii) = 1;
            end
            
            if(z_sp(ii) < 1)
                x_sp(ii) = 1;
            end
            
            if(x_sp(ii) > W)
                x_sp(ii) = W;
            end
            
            if(y_sp(ii) > H)
                y_sp(ii) = H;
            end
            
            if(z_sp(ii) > D)
                z_sp(ii) = D;
            end
            
            BW(x_sp(ii), y_sp(ii), z_sp(ii)) = 1;
        end
        E = bwdistsc(BW, [vx, vy, vz]);
        
        %% Obtain background and freground coordinates with downsampling        
        is_bg = false(numel(x_sampling),1);
        is_fg = false(numel(x_sampling),1);
        
        for ii = 1:numel(x_sampling)
            is_bg(ii) =  E(x_sampling(ii), y_sampling(ii), z_sampling(ii)) > TH_DISTANCE_BG;
            is_fg(ii) =  E(x_sampling(ii), y_sampling(ii), z_sampling(ii)) < TH_DISTANCE_FG;
        end

        x_bg = x_sampling(is_bg);
        y_bg = y_sampling(is_bg);
        z_bg = z_sampling(is_bg);
        num_samples_bg = numel(x_bg);
        
        x_fg = x_sampling(is_fg);
        y_fg = y_sampling(is_fg);
        z_fg = z_sampling(is_fg);
        num_samples_fg = numel(x_fg);
        
        if((num_samples_fg < 3) || (num_samples_bg < 3))
            disp([' --- not enough points in tf ', num2str(tf)]);
            continue
        end
        
        BG_values = zeros(num_samples_bg, C);
        for ii = 1:num_samples_bg
            BG_values(ii,:) = zstack(x_bg(ii), y_bg(ii), z_bg(ii), :, tf);
        end
        
        FG_values = zeros(num_samples_fg, C);
        for ii = 1:num_samples_fg
            FG_values(ii,:) = zstack(x_fg(ii), y_fg(ii), z_fg(ii), :, tf);
        end
            
        FG_avg = mean(FG_values); %FG mean of columns (for each channel)
        FG_std = std(FG_values);  %FG std of columns (for each channel)
        BG_avg = mean(BG_values); %BG mean of columns (for each channel)
        BG_std = std(BG_values);  %BG std of columns (for each channel)

        SNR = abs(FG_avg - BG_avg) ./ abs(BG_std);
        CR = FG_avg ./ BG_avg;
        HET = FG_std ./ abs(FG_avg - BG_avg);
        
        SNR_T(tf,:) = SNR;
        CR_T(tf,:) = CR;
        HET_T(tf,:) = HET;
        FG_T(tf,:) = FG_avg;
        
    end
    close(h)

    DEN_T(DEN_T == 9999) = NaN;
    mean_SNR = nanmean(SNR_T);
    mean_CR = nanmean(CR_T);
    mean_HET = nanmean(HET_T);
    mean_DEN = nanmean(DEN_T);
    std_DEN = nanstd(DEN_T);
    mean_NUM = nanmean(NUM_T);
    std_NUM = nanstd(NUM_T);
    
    rate_BLEACH = zeros(1,C);
    
    Headings = {};
    for cc = 1:C
        Headings{end+1} = ['SNR C', num2str(cc)];
    end
    for cc = 1:C
        Headings{end+1} = ['CR C', num2str(cc)];
    end
    for cc = 1:C
        Headings{end+1} = ['HET C', num2str(cc)];
    end
    
    Headings{end+1} = 'DEN';
    Headings{end+1} = 'NUM';
    
    
    Values = [SNR_T, CR_T, HET_T, DEN_T, NUM_T];
    
    Totals = [mean_SNR, mean_CR, mean_HET, mean_DEN, mean_NUM];
    
    xlsData = [Headings; num2cell(Totals); num2cell(Values)];

    if(exist([LTDB_PATH_MEASURES, file.name, '_DS_MEASURES_', '.xls'], 'file'))
        delete([LTDB_PATH_MEASURES, file.name, '_DS_MEASURES_', '.xls']);
    end
    xlswrite([LTDB_PATH_MEASURES, file.name, '_DS_MEASURES_', '.xls'], xlsData);
    disp(' ---- DONE ---- ');
    
    overall(ifiles,1:numel(mean_SNR)) = mean_SNR;
    overall(ifiles,5:5+numel(mean_CR)-1) = mean_CR;
    overall(ifiles,9) = mean_DEN;
    overall(ifiles,10) = std_DEN;
    overall(ifiles,11) = mean_NUM;
    overall(ifiles,12) = std_NUM;
    overall(ifiles,13) = W;
    overall(ifiles,14) = H;
    overall(ifiles,15) = D;
    overall(ifiles,16) = C;
    overall(ifiles,17) = T;
    overall(ifiles,18) = vx;
    overall(ifiles,19) = vz;
    overall(ifiles,20) = dt;
    overall(ifiles,21) = numel(unique(spots_IXYZT(:,1)));
    overall_names{ifiles} = file.name;
end

%% Writing overall measures for the entire dataset
overall = overall(1:ifiles, :);
overall_names = overall_names(1:ifiles, :);

Headings2 = {'NAME'};
Headings2{end+1} = 'SNR CH0';
Headings2{end+1} = 'SNR CH1';
Headings2{end+1} = 'SNR CH2';
Headings2{end+1} = 'SNR CH3';
Headings2{end+1} = 'CR CH0';
Headings2{end+1} = 'CR CH1';
Headings2{end+1} = 'CR CH2';
Headings2{end+1} = 'CR CH3';
Headings2{end+1} = 'DEN avg';
Headings2{end+1} = 'DEN std';
Headings2{end+1} = 'NUM avg';
Headings2{end+1} = 'NUM std';
Headings2{end+1} = 'W';
Headings2{end+1} = 'H';
Headings2{end+1} = 'D';
Headings2{end+1} = 'C';
Headings2{end+1} = 'T';
Headings2{end+1} = 'dxy';
Headings2{end+1} = 'dz';
Headings2{end+1} = 'dt';
Headings2{end+1} = 'N.TRACKS';

xlsData2 = [Headings2; overall_names, num2cell(overall)];

if(exist([LTDB_PATH_MEASURES, 'overall.xls'], 'file'))
    delete([LTDB_PATH_MEASURES, 'overall.xls']);
end
xlswrite([LTDB_PATH_MEASURES, 'overall.xls'], xlsData2);