ttc = 0;
min_track_id = min(curr_TrackID);


% Find the TrackIDs
curr_ids = zeros(size(curr_coordinatesXYZ, 1), 1)-1;
curr_edges = curr_edges + 1;
for kk = 1:size(curr_edges, 1)
    edge_start = curr_edges(kk,1);
    edge_end = curr_edges(kk,2);
   
    %if two tracklets are connected
    if(curr_ids(edge_start)>0 && curr_ids(edge_end)>0)
        conn_tracklets_idx = (curr_ids == curr_ids(edge_start));
        curr_ids(conn_tracklets_idx) = curr_ids(edge_end);
    end
    
    %if a point is connected
    if(curr_ids(edge_start) < 0)
        if(curr_ids(edge_end) < 0)
            ttc = ttc + 1;
            curr_ids(edge_start) = ttc;
        else
            curr_ids(edge_start) = curr_ids(edge_end);
        end
    end
    curr_ids(edge_end) = curr_ids(edge_start);
end

% Move the origin in (0,0,0) (if videos have different origins)
curr_coordinatesXYZ(:,1) = curr_coordinatesXYZ(:,1) - xmin;
curr_coordinatesXYZ(:,2) = curr_coordinatesXYZ(:,2) - ymin;
curr_coordinatesXYZ(:,3) = curr_coordinatesXYZ(:,3) - zmin;

% A contains track ID, X, Y, Z, time
A = [double(curr_ids), double(curr_coordinatesXYZ), double(curr_indicesT)+1];

%% Glitch correction and continuity check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_tracks_init = numel(unique(A(:,1)));
T = sizet;

M = zeros(num_tracks_init, T, 3)+NaN; % each row is a track. each column a time point. pages: (x,y,z) coordinates
M_original_id = zeros(num_tracks, 1);
num_tracks = 0;
track_ids = unique(A(:,1));
for curr_track_id = track_ids'
   num_tracks = num_tracks + 1;
   curr_track_idx = find(A(:,1) == curr_track_id);
   curr_track_pos = A(curr_track_idx,:);

   for curr_pos = curr_track_pos'
       M(num_tracks, curr_pos(5), :) = curr_pos(2:4);
   end
   M_original_id(num_tracks) = curr_track_id;
end   


for ii = 1:num_tracks
	tbegin = find(~isnan(M(ii,:,1)), 1);
    tend = find(~isnan(M(ii,:,1)), 1, 'last');
    if (sum(~isnan(M(ii,:,1))) ~= ((tend - tbegin)+1))
        disp([filename,' ', curr_name, ' Track ', num2str(ii), ' DISCONTINUOUS']);
    end
    
    for tt = 2:T-1
        p = squeeze(M(ii,tt-1,:))';
        c = squeeze(M(ii,tt,:))';
        q = squeeze(M(ii,tt+1,:))';
        
        pq = pdist2(p,q);
        pc = pdist2(p,c);
        cq = pdist2(c,q);
        
        if ((min(pc,cq) > 10) & (pq < 10))
            M(ii,tt,:) = mean([p;q]);
            disp([filename, ' - ', curr_name, ' - CORRECTED GLITCH: Track ', num2str(ii), ' - time ', num2str(tt), ' ', num2str(min(pc,cq)), 'um']);
        end
    end
end

B = [];
for ii = 1:num_tracks
    for tt = 1:T
        if(~isnan(M(ii,tt,1)))
            B = [B; ii, M(ii,tt,1), M(ii,tt,2), M(ii,tt,3), tt];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Non-merged tracks check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:num_tracks
    curr_track = M(ii,:,:);
    curr_track_duration = sum(~isnan(M(ii,:,1)));
    curr_track_start = find(~isnan(M(ii,:,1)), 1);
    curr_track_end = find(~isnan(M(ii,:,1)), 1, 'last');
    
    D = sqrt(sum((M - curr_track).^2, 3));
    O = D < TH_DIST_OVERLAP; %boolean matrix 1 is an overlapping
    O_count = sum(O, 2);
    IS_OV = O_count > min(curr_track_duration, TH_TIME_OVERLAP);
    if(sum(IS_OV) > 1)
        for jj = 1:numel(IS_OV)
            if (IS_OV(jj)) && (ii ~= jj)
                disp([filename, ' - ', curr_name, ' - WARNING OVERLAP: Track xyz:', num2str(ii), '(', num2str(M(ii,curr_track_start,:)), ' Tstart:', num2str(curr_track_start) ,' Tdur:', num2str(curr_track_duration) ') overlaps with track ', num2str(M_original_id(jj)), ' for ', num2str(O_count(jj)), ' instants ']),
            end
        end
    end
    
    % Short tracks
    if(curr_track_duration <= 2)
        disp([filename, ' - ', curr_name, ' - WARNING SHORT: Track xyz:', num2str(ii), '(', num2str(M(ii,curr_track_start,:)), ' Tstart:', num2str(curr_track_start) ,' Tdur:', num2str(curr_track_duration) ')']),
    end
    
    % Broken tracks
    end_pos = squeeze(M(ii,curr_track_end,:));
    if(curr_track_end < T && num_tracks > 1)
        D_af_end = pdist2(end_pos', squeeze(M(:,curr_track_end+1,:)));
        idx_af_end = find(~isnan(D_af_end) & (D_af_end < TH_DIST_OVERLAP) & isnan(squeeze(M(:,curr_track_end,1))'));
        if(sum(idx_af_end > 1))
            disp([filename, ' - ', curr_name, ' - WARNING BROKEN: Track xyz:', num2str(M_original_id(ii)), '(', num2str(M(ii,curr_track_start,:)), ' Tstart:', num2str(curr_track_start) ,' Tdur:', num2str(curr_track_duration) ')']),
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Wrongly recognized coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:num_tracks
    curr_track_duration = sum(~isnan(M(ii,:,1)));
    curr_track_start = find(~isnan(M(ii,:,1)), 1);
    curr_track_end = find(~isnan(M(ii,:,1)), 1, 'last');
    curr_track = M(ii,:,:);
    for tt = curr_track_start:curr_track_end
        if((M(ii,tt,1) == M(ii,tt,2)) && (M(ii,tt,2) == M(ii,tt,3)))
            disp([filename, ' - ', curr_name, ' ERROR: Tracks exported with wrong Imaris version ']);
            return;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Identical points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:num_tracks
    curr_track_duration = sum(~isnan(M(ii,:,1)));
    curr_track_start = find(~isnan(M(ii,:,1)), 1);
    curr_track_end = find(~isnan(M(ii,:,1)), 1, 'last');
    curr_track = M(ii,:,:);
    D = sqrt(sum((M - curr_track).^2, 3));
    O = (D <= TH_EXACT);
    O(ii,:) = 0;
    
    if sum(O(:)) > 0
        disp([filename, ' - ', curr_name, ' - WARNING EXACT: Track xyz:', num2str(ii), '(', num2str(M(ii,curr_track_start,:)), ' Tstart:', num2str(curr_track_start) ,' Tdur:', num2str(curr_track_duration) ')']),
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [0, xVoxelSize, yVoxelSize, zVoxelSize, dt; ...
     0  0           0           0           0;  ...
     B];

lcn = numel(curr_name);
curr_name = curr_name(1:min([lcn,10]));

if(exist([filename, '_', curr_name, '.xls'], 'file'))
    delete([filename, filename, '_', curr_name, '.xls']);
end
xlswrite([filename, '_', curr_name, '.xls'], A);