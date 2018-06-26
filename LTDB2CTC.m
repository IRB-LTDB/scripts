%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exports LTDB tracking ground truth in the Cell Tracking Fromat[1,2]     %
% Approximates segmentation ground truth by thresholding the distance     %
% transform. i.e. it creates 3d spheres centered in centroid i            %
% whose radius is the minimum  between TH_DIST_TRANSFORM and the distance %
% to the closer centroid j / 2                                            %
% Notes:                                                                  %
% - Requires imaging data as Imaris files and tracks in the LTDB format   %
% - Distance transform is computed using [3]                              %
% - Imaris files are interpreted using [4]                                %
% 1. Ulman et. al. Nature Methods, 2017                                   %
% 2. Maska et. al. Bioinformatics, 2014                                   %
% 3. https://ch.mathworks.com/matlabcentral/fileexchange/15455-3d-        %
%    euclidean-distance-transform-for-variable-data-aspect-ratio          %
% 4. https://github.com/PeterBeemiller/ImarisReader                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Settings - Edit here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LTDB_PATH_IMS = 'G:\LTDB_r2\GT_IMS\';                                     %
LTDB_PATH_TRACKS = 'G:\LTDB_r2\GT_TRACKS\';                               %
LTDB_PATH_CTC_GT = 'G:\LTDB_r2\format_ctc\GT\';                           %
CURR_OP = 'GT';                                                           %
EN_DEBUG = false;                                                         %
TH_DIST_TRANSFORM = 10;                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Included libraries
addpath '.\libs\ImarisReader-master';
addpath '.\libs\bwdistsc';

%Processing all the files in the folder
files = dir([LTDB_PATH_TRACKS, '*',CURR_OP,'.xls']);

%% Export as CTC format
for file = files'
    fns = split(file.name,['_',CURR_OP]);
        
    curr_ctc_dir = [LTDB_PATH_CTC_GT,num2str(fns{1},'%02d'),'_GT\TRA'];
    if(exist(curr_ctc_dir, 'dir'))
        rmdir(curr_ctc_dir, 's')
    end
    mkdir(curr_ctc_dir);
    dlmwrite([curr_ctc_dir,'\',file.name,'.txt'], ' ');
    
    [~,fname,fext] = fileparts(file.name);
    disp(['Reading ', fname]);
    fns = split(fname, '_');
    fn_ims = [fns{1}, '_GT.ims'];
    try
    curr_ims = ImarisReader([LTDB_PATH_IMS, fn_ims]);
    catch
        disp(['Error initializing in ', fn_ims]);
        continue;
    end
    
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
    
    if(vx ~= vy)
        disp([file.name, ' ERROR: vx=', num2str(vx), ' vy=', num2str(vy), ' vz=', num2str(vz), ' USING vy = vx = mean(vx,vy)']);
        vy = mean([vx, vy]);
        vx = vy;
    end
    
    tot_td = 0;
    num_tracks = 0;    
    
    spots_IXYZT = xlsread([LTDB_PATH_TRACKS,file.name]);
    spots_IXYZT_um = spots_IXYZT;
    
    spots_IXYZT(3:end, 2) = spots_IXYZT(3:end, 2)./spots_IXYZT(1, 2);
    spots_IXYZT(3:end, 3) = spots_IXYZT(3:end, 3)./spots_IXYZT(1, 3);
    spots_IXYZT(3:end, 4) = spots_IXYZT(3:end, 4)./spots_IXYZT(1, 4);
    
    spots_IXYZT = spots_IXYZT(3:end,:);
    spots_IXYZT_um = spots_IXYZT_um(3:end,:);
    tic;
    h = waitbar(0, ['Processing ', file.name]);
    for tf = 1:T
        waitbar(tf/T, h);
        zstack_label = zeros(W,H,D, 'uint16');
        
        spots_idx = find(spots_IXYZT(:,5) == tf);
        if(isempty(spots_idx))
            for cz=1:D
                zstack_label(1,1,1) = 1;
                imwrite(zstack_label(:,:,cz), [curr_ctc_dir,'\man_track',num2str(tf-1,'%03d'),'.tif'],'WriteMode','append');
            end
            continue
        end
        
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
        
        for ii = 1:numel(x_sp)
            cx_sp = x_sp(ii);
            cy_sp = y_sp(ii);
            cz_sp = z_sp(ii);
            cl_sp = l_sp(ii);
            BW = false(W,H,D);
            BW(cx_sp, cy_sp, cz_sp) = 1;
            
            CE = bwdistsc(BW, [vx, vy, vz]);
            
            lbidx = ((abs(CE - E) < 0.1) & (CE < TH_DIST_TRANSFORM));
            zstack_label(lbidx) = cl_sp+1;
        end
        
        zstack_label(1,1,1) = 1;
        for cz=1:D
            imwrite(zstack_label(:,:,cz), [curr_ctc_dir,'\man_track',num2str(tf-1,'%03d'),'.tif'],'WriteMode','append');
        end
        %disp(['WRITTEN ',curr_ctc_dir,'\man_track',num2str(tf-1,'%03d'),'.tif']);
    end
    toc;
    close(h);
    track_ids = unique(spots_IXYZT(:,1));
    
    aog = zeros(numel(track_ids), 4, 'uint16');
    
    aog(1,1) = 1;
    aog(1,2) = 0;
    aog(1,3) = T-1;
    aog(1,4) = 0;
    for ii = 1:numel(track_ids)
        curr_idx = find(spots_IXYZT(:,1) == track_ids(ii));
        curr_track_start = min(spots_IXYZT(curr_idx,5))-1;
        curr_track_end = max(spots_IXYZT(curr_idx,5))-1;
        %curr_track_id = max(spots_IXYZT(curr_idx,1));
        aog(ii+1,1) = track_ids(ii)+1;        
        aog(ii+1,2) = curr_track_start;
        aog(ii+1,3) = curr_track_end;
    end
    
    dlmwrite([curr_ctc_dir,'\man_track.txt'], aog, 'delimiter', ' ', 'newline', 'pc');    
end