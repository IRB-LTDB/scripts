function [zstack] = LTDBReadImages(folderpath, channel_to_read, showbar)
%Reads an LTDB video from tif files (single channel).
%Requires the path to the folder containing all the files, the number of
%channel to read (starting from 1), and a boolean parameter to display the
%progress bar.
%returns a 4D zstack (XYZT) as uint16.
%Example: zstack = LTDBReadImages([LTDB_TIFF_PATH , 'LTDB001'], 1, 1);
    
    if ((folderpath(end) == '/') || (folderpath(end) == '\'))
        folderpath = folderpath(1:end-1);
    end
    
    files = dir([folderpath,'/*.tif']);
    
    Z = 0;
    C = 0;
    T = 0;
    
    for file = files'
        currfn = file.name;
        parts = strsplit(currfn,'_');
        curr_z = parts{2};
        curr_c = parts{3};
        curr_t = parts{4};
        
        if(~isequal(curr_t(1), 'T') || ~isequal(curr_z(1), 'Z') || ~isequal(curr_c(1), 'C'))
            disp(['Warning: ignored ',file.name]);
            continue
        end
        
        curr_z = str2num(curr_z(2:end));
        curr_c = str2num(curr_c(2:end));
        curr_t = str2num(curr_t(2:end-4));
        
        Z = max(Z, curr_z);
        C = max(C, curr_c);
        T = max(T, curr_t);
    end
    
    Z = Z+1;
    C = C+1;
    T = T+1;
    Itemp = imread([folderpath,'/',currfn]);
    H = size(Itemp,1);
    W = size(Itemp,2);
    
    zstack = zeros(H,W,Z,T,'uint16');
    
    if(showbar)
        h = waitbar(0, 'Reading dataset');
        totfiles = numel(files');
    end
    countfiles = 0;
    for file = files'
        currfn = file.name;
        parts = strsplit(currfn,'_');
        currZ = parts{2}; currZ = str2num(currZ(2:end))+1;
        currC = parts{3}; currC = str2num(currC(2:end))+1;        
        currT = parts{4}; currT = str2num(currT(2:end-4))+1;
        if (currC == channel_to_read) 
            if (currT <= T) && (currZ <= Z) && (currC <= C)
                Itemp = imread([folderpath,'/',currfn]);
                if(size(Itemp) == size(zstack(:,:,currZ,currT)))
                    zstack(:,:,currZ,currT) = Itemp;
                else
                    errordlg(['Error in ', file.name]);
                    return;
                end
            else
                disp 'Warning: there are files not included in the size of the datset';
            end
        end
        countfiles = countfiles + 1;
        if(showbar && (totfiles > 0))
            waitbar(countfiles / totfiles, h);
        end
    end
    if(showbar)
        close(h);
    end
    
end