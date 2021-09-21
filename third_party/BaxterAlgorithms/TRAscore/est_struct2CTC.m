function est_struct2CTC(est, output_path, width, height)
    K = length(est.N) ; 
    x_lim = width ; 
    y_lim = height ; 
    res_track = zeros(0,4) ; 
    counted_L = [] ; 
    % count labels
    [~, strlabel]= countestlabels(est.L) ; 
    fol_name = fullfile(output_path, 'res') ; 
    if ~isfolder(fol_name)
        mkdir(fol_name) ; 
    end
    for k = 1 : K
        img_path = fullfile(fol_name, ['mask', sprintf('%03d', k-1), '.tif']) ; 
        curr_img = zeros(y_lim, x_lim) ; 
        for n = 1 : est.N(k)
            curr_lab = sprintf('%d-' , est.L{k}{n});
            curr_idx = find(strcmp(strlabel, curr_lab)) ; 
            location = round(est.X{k}([1 3],n)) ; 
            if location(1)<1 
                location(1) = 1 ; 
            elseif location(1)>x_lim
                location(1) = x_lim ; 
            end
            if location(2)<1
                location(2) = 1 ; 
            elseif location(2)>y_lim
                location(2) = y_lim ; 
            end
            % write the current track index to the pixel location
            curr_img(location(2), location(1)) = curr_idx ; 
            if isempty(find(strcmp(counted_L, curr_lab),1)) 
                res_track = [res_track ; [curr_idx , k-1, 0, 0] ] ; 
                counted_L = [counted_L, {curr_lab}] ; 
                % then find the parent
                if length(est.L{k}{n})>2
                    parent_lab = sprintf('%d-', est.L{k}{n}(1:end-2)) ; 
                    parent_idx = find(strcmp(strlabel, parent_lab)) ; 
                    res_track(end, 4) = parent_idx ; 
                end
            else
                temp_idx = find(res_track(:,1) == curr_idx) ; 
                res_track(temp_idx, 3) = k-1 ; % update the dead time
            end
        end
        imwrite(uint16(curr_img), img_path) ; 
    end
    temp_indc = find(res_track(:,3) == 0) ; 
    res_track(temp_indc,3) = res_track(temp_indc,2) ; 
    filename = fullfile(fol_name, 'res_track.txt') ; 
    dlmwrite(filename, res_track, 'delimiter', ' ') ; 
end

function [count, c]= countestlabels(est_L)
    strlabel = [] ; 
    for k=1:length(est_L)
        tmp_str_label = cell(1 , length(est_L{k}) ) ;
        for lidx = 1 : length(est_L{k})
            tmp_str_label{lidx}= sprintf('%d-' , est_L{k}{lidx}');
        end
        strlabel = [strlabel, tmp_str_label] ;
    end
    [c,~,~]= unique(strlabel ,'rows');
    count=size(c,2);
end
