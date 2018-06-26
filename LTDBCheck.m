%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks for common tracking errors and Exports spots/surfaces coordinates%
% from Imaris IMS files                                                   %
% Checks for                                                              %
% - DISCONTINUOUS tracks                                                  %
% - GLITCHES                                                              %
% - EXACT two spots (useful when a cell has been annotated two times      %
% in the same frame)                                                      %
% OVERLAP (same tracks by 2 operators or close cells)                     %
% These checks are included in doExportTracks.m                           %
% Requires ImarisReader library                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all
clc

%% Settings - Edit here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LTDB_PATH_IMS = 'G:\LTDB_r2\GT_IMS\';
LTDB_PATH_TRACKS = 'G:\LTDB_r2\GT_TRACKS\';
TH_DIST_OVERLAP = 10;  %[um]
TH_TIME_OVERLAP = 10;  % time instants
TH_EXACT = 1; %[um]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Libraries included
addpath '.\libs\ImarisReader-master';

%files = dir([LTDB_PATH_IMS, '*.ims']);
[FileName,PathName,FilterIndex] = uigetfile('*.ims');
LTDB_PATH_IMS = [PathName,'\'];
LTDB_PATH_TRACKS = [PathName, '\exported_tracks\'];
filename = [PathName,FileName];

tot_td = 0;
num_tracks = 0;

%Uncomment for batch processing
% for file = files'
%    filename = file.name;
    curr_ims = ImarisReader(filename);
    
    xmin = curr_ims.DataSet.ExtendMinX;
    ymin = curr_ims.DataSet.ExtendMinY;
    zmin = curr_ims.DataSet.ExtendMinZ;
    
    xmax = curr_ims.DataSet.ExtendMaxX;
    ymax = curr_ims.DataSet.ExtendMaxY;
    zmax = curr_ims.DataSet.ExtendMaxZ;
    
    sizex = curr_ims.DataSet.SizeX;
    sizey = curr_ims.DataSet.SizeY;
    sizez = curr_ims.DataSet.SizeZ;    
    sizec = curr_ims.DataSet.SizeC;    
    sizet = curr_ims.DataSet.SizeT;
    
    xVoxelSize = (xmax-xmin)/sizex;
    yVoxelSize = (ymax-ymin)/sizey;
    zVoxelSize = (zmax-zmin)/sizez;
    
    dt = seconds(datetime(curr_ims.DataSet.Timestamps{2}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS') - datetime(curr_ims.DataSet.Timestamps{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS'));
    
    num_spots = curr_ims.Scene.NumberOfSpots;
    
    num_surfaces = curr_ims.Scene.NumberOfSurfaces;
    
    if numel(curr_ims.Surfaces) ~= num_surfaces
        num_surfaces = numel(curr_ims.Surfaces);
        disp(['ERROR: Misisng surface in ', filename, '\n']);
    end
    
    if numel(curr_ims.Spots) ~= num_spots
        num_spots = numel(curr_ims.Spots);
        disp(['ERROR: Missing spots in ', filename, '\n']);
    end
    
    for ii = 1:num_spots
       curr_spots = curr_ims.Spots(ii);
       curr_name = curr_spots.Name; 
       curr_coordinatesXYZ = curr_spots.GetPositions;
       curr_indicesT = curr_spots.GetIndicesT;
       curr_edges = curr_spots.GetTrackEdges;
       curr_TrackID = curr_spots.GetTrackIDs;
       
       doExportTracks
    end
    
    for ii = 1:num_surfaces
       curr_surfaces = curr_ims.Surfaces(ii);
       curr_name = curr_surfaces.Name; 
       curr_coordinatesXYZ = curr_surfaces.GetPositions;
       curr_indicesT = curr_surfaces.GetIndicesT;
       curr_edges = curr_surfaces.GetTrackEdges;
       curr_TrackID = curr_surfaces.GetTrackIDs;
       
       doExportTracks
    end
%end
%Uncomment for batch processing

