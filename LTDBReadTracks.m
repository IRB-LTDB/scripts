function [spots_IXYZT, channel_visibility, voxel_size] = LTDBReadTracks(filename)
%LTDBREADTRACKS Reads the tracks in the LTDB CSV Format
%   Returns:
%   - spots_IXYZT: Nx5 matrix (Track_ID[int], x [um], y[um], z[um], t[frame n.])
%   - channel_visibility: 1x5 matrix (bool) 1 if visible in the
%   corresponding channel
%   - voxel size (x,y,z [um], t[time instant])
    fileContent = dlmread(filename,';', 0, 1);
    voxel_size = fileContent(1, 2:end);
    
    fileContent = dlmread(filename,';', 1, 0);
    channel_visibility = logical(fileContent(1,:));
    
    fileContent = dlmread(filename,';', 2, 0);
    spots_IXYZT = fileContent(3:end, 1:5);
end
