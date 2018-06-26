%
%  Imports spots from a CSV file in the LTDB format.
%  
%  Revision: RC1 20180626
%
%  Author: Pizzagalli Diego Ulisse (1,2)
%          1. Institute for Research in Biomedicine - Bellinzona (CH)
%          2. Institute of Computational Science,
%             Universita della Svizzera Italiana - Lugano (CH)
%
%  INSTALLATION ON COMPUTERS WITH MATLAB ALREADY INSTALLED:
%                1. Copy this file into an XTensions folder
%                   (e.g C:\Program Files\Bitplane\Imaris [ver]\XT\matlab).
%                2. Restart Imaris and you can find this function in the 
%                   Image Processing menu with the "SVMColoc" name.
%                NOTE: Tested with MATLAB r2012b, r2013a
%
%  INSTALLATION ON COMPUTERS WITHOUT MATLAB:
%                1. Copy the compiled files SVMColor.exe and SVMColoc.xml
%                   into the compiled XTensions folder
%                   (e.g C:\Program Files\Bitplane\Imaris [ver]\XT\rtmatlab).
%                2. Restart Imaris and you can find this function in the 
%                   Image Processing menu with the "SVMColoc" name.
%                NOTE: You need the Matlab Compiler Runtime to be installed
%                      Please contact Bitplane for help.
%
%    <CustomTools>
%      <Menu>
%        <Item name="IRB --> Imports spots from CSV" icon="I"
%        tooltip="Imports spots (and relative tracks) from a CSV file">
%          <Command>MatlabXT::XTLTDBImportTracks(%i)</Command>
%        </Item>
%      </Menu>
%    </CustomTools>


function XTLTDBImportTracks(aImarisApplicationID)
    if isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
        vImarisApplication = aImarisApplicationID;
    else
        javaaddpath ImarisLib.jar;
        vImarisLib = ImarisLib;
        if ischar(aImarisApplicationID)
            aImarisApplicationID = round(str2double(aImarisApplicationID));
        end
        vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
    end
    
    COL_TID = 1;
    COL_X   = 2;
    COL_Y   = 3;
    COL_Z   = 4;
    COL_T   = 5;
    
    
    
    [FileName,PathName] = uigetfile({'*.csv;'}, ...
    'Pick CSV file with tracks');
    
    fn = [PathName,FileName];
    
    fileContent = dlmread(fn,';', 0, 1);
    voxel_size = fileContent(1, 2:end);
    
    fileContent = dlmread(fn,';', 1, 0);
    channel_visibility = logical(fileContent(1,:));
    
    fileContent = dlmread(fn,';', 2, 0);
    spots_IXYZT = fileContent(3:end, 1:5);
    
    vSpots = vImarisApplication.GetFactory.CreateSpots;
    aPositionsXYZ=spots_IXYZT(:,[COL_X,COL_Y,COL_Z]);
    aIndicesT=spots_IXYZT(:, COL_T)-1';
    aRadii=aIndicesT*0+1;
    vSpots.Set(aPositionsXYZ,aIndicesT,aRadii);
    spots_ids = unique(spots_IXYZT(:,COL_TID));
    aEdges = [];
    for ii = 1:numel(spots_ids)
        curr_id = spots_ids(ii);
        idx = find(spots_IXYZT(:,COL_TID) == curr_id);
        t = spots_IXYZT(idx, COL_T);
        [~, I] = sort(t);
        for jj = 1:numel(I)-1
            aEdges = [aEdges; idx(I(jj)), idx(I(jj+1))];
        end
    end
    vSpots.SetTrackEdges(aEdges-1);
    vImarisApplication.GetSurpassScene.AddChild(vSpots, -1);
    return;