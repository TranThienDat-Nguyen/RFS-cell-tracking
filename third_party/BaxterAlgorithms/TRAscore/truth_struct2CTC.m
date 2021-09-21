function truth_struct2CTC(truth, output_path, width, height)
    folname = fullfile(output_path, 'truth', 'TRA') ; 
    if ~isfolder(folname)
        mkdir(folname) ; 
    end

    man_track = zeros(0,4) ; 
    counted_L = [] ; 
    for k = 1 : truth.K
        img_path = fullfile(folname, ['man_track', sprintf('%03d', k-1), '.tif']) ; 
        curr_img = zeros(height, width) ; 
        for n = 1 : truth.N(k)
            curr_lab = truth.track_list{k}(n) ; 
            location = round(truth.X{k}([1 3], n)) ; 
            curr_img(location(2),location(1)) = curr_lab ; 
            if ~ismember(curr_lab, counted_L)
                counted_L = [counted_L, curr_lab] ; 
                % check if this is a spawned track
                [isspawn, spawn_loc] = ismember(curr_lab, truth.newSpawn{k}) ; 
                man_track = [man_track ; [curr_lab, k-1, 0, 0] ] ;
                if isspawn
                    % try to find the parent 
                    if mod(spawn_loc,2) == 0 
                        idx = spawn_loc/2 ; 
                    else
                        idx = ceil(spawn_loc/2) ; 
                    end
                    parent = truth.spawnNextTimeTrack{k-1}(idx) ; 
                    man_track(end,4) = parent ; 
                end
            else
                % update the dead time
                temp_idx = find(man_track(:,1) == curr_lab) ; 
                man_track(temp_idx, 3) = k-1 ; 
            end
        end
        imwrite(uint16(curr_img), img_path) ;
    end
    temp_indc = find(man_track(:,3) == 0) ; 
    man_track(temp_indc,3) = man_track(temp_indc,2) ; 
    man_file_path = fullfile(folname, 'man_track.txt') ; 
    dlmwrite(man_file_path, man_track, 'delimiter', ' ') ; 
end
        
