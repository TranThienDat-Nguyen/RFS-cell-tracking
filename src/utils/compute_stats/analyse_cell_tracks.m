function [tracks_img , trans_tracks_img , transdis] = analyse_cell_tracks (xy_pos , labels ,K, backgroundImg)
    labelcount = countlabels (labels) ;
    colorarray = makecolorarray(labelcount);
    % Convert labels vectors to strings
    labelspool = [] ; 
    for k = 1 : K
        for t = 1 : length(labels{k})
            temp_label = sprintf('%.0f-' ,cell2mat( labels{k}(t) ) );
            labelspool = [labelspool {temp_label}];
        end
    end
    uniquepool = unique(labelspool) ;
    % Extract the tracks
    labels_out.labels = cell(length(uniquepool) , 1) ;
    labels_out.time = zeros(length(uniquepool) , 1);
    legends = [] ;
%     h = figure('Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    h = figure ; 
    imshow(backgroundImg{end});
    hold all ;
    LW = 2 ; 
    for t = 1 : length(uniquepool)
        temp_time = [] ;
        labels_out(t).labels = uniquepool{t}(1:end-1) ;
        legend_str = uniquepool{t}(1:end-1) ; 
        legends = [legends {legend_str}] ;
        track = [] ;
        IgotColor = false ;
        for k = 2 : K
            curr_labels = [] ;
            % Convert current time labels to strings
            for tk = 1 : length(labels{k})
                temp_label = sprintf('%.0f-' ,cell2mat( labels{k}(tk) ) );
                curr_labels = [curr_labels {temp_label}];
            end
            % Compare the checking track with the current time labels
            tidx = find( ismember(curr_labels,cell2mat(uniquepool(t))) );
            if ~isempty(tidx)
                temp_time = [temp_time k];
                track = [track xy_pos{k}(:,tidx)] ;
                if ~IgotColor
                    color_opt = colorarray.rgb(assigncolor(labels{k}(:,tidx)),:) ;
                    IgotColor = true ;
                end
            end      
        end
        labels_out(t).time = temp_time ;
        if ~isempty(track)
            plot(track(1,:) , track(2,:) , 'Color' , color_opt, 'LineWidth', LW) ;
        end
%         legend(legends) ;
%         set(legend, 'NumColumns' ,4);
        set(gca,'FontSize',18);
    end
    hold off ;
    tracks_img = getframe(gca) ;
    close(h) ;
    % Get translated tracks
    legends = [] ;
    h = figure('Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
    hold all ;
    transdis = zeros(length(uniquepool) , 2);
    for t = 1 : length(uniquepool)
        IgotFirstpos = false ;
        Fisrtpos = [0 ; 0] ;
        legend_str = uniquepool{t}(1:end-1) ; 
        legends = [legends {legend_str}] ;
        track = [0 ; 0] ;
        IgotColor = false ;
        for k = 2 : K
            curr_labels = [] ;
            % Convert current time labels to strings
            for tk = 1 : length(labels{k})
                temp_label = sprintf('%.0f-' ,cell2mat( labels{k}(tk) ) );
                curr_labels = [curr_labels {temp_label}];
            end
            % Compare the checking track with the current time labels
            tidx = find( ismember(curr_labels,cell2mat(uniquepool(t))) );
            if ~isempty(tidx)
                if ~IgotFirstpos
                    Firstpos = xy_pos{k}(:,tidx);
                    IgotFirstpos = true ;
                else
                    track = [track xy_pos{k}(:,tidx)-Firstpos];
                end
 
                if ~IgotColor
                    color_opt = colorarray.rgb(assigncolor(labels{k}(:,tidx)),:) ;
                    IgotColor = true ;
                end
            end      
        end
        if ~isempty(track)
            plot(track(1,:) , track(2,:) , 'Color' , color_opt, 'LineWidth', LW) ;
        end
%         legend(legends) ;
%         set(legend, 'NumColumns' ,3)
        set(gca,'FontSize',18);
        transdis(t,1) = track(1,end) ;
        transdis(t,2) = track(2,end) ;
    end
    hold off ;
    trans_tracks_img = getframe(gca) ;
    close(h) ;
    
function ca= makecolorarray(nlabels)
lower= 0.1;
upper= 0.9;
rrr= rand(1,nlabels)*(upper-lower)+lower;
ggg= rand(1,nlabels)*(upper-lower)+lower;
bbb= rand(1,nlabels)*(upper-lower)+lower;
ca.rgb= [rrr; ggg; bbb]';
ca.lab= cell(nlabels,1);
ca.cnt= 0;   
end

function idx= assigncolor(label)
    str= sprintf('%.0f-' ,label{1});
    str= str(1:end-1);
    tmp= strcmp(str,colorarray.lab);
    if any(tmp)
        idx= find(tmp);
    else
        colorarray.cnt= colorarray.cnt + 1;
        colorarray.lab{colorarray.cnt}= str;
        idx= colorarray.cnt;
    end
end

function count= countlabels(labels)
    strlabel= [];
    for k = 1 : K
        tmp_str_label = cell(1 , length(labels{k}) ) ;
        for lidx = 1 : length(labels{k})
            tmp_str_label{lidx}= sprintf('%.0f,' ,labels{k}{lidx}');
        end
        strlabel = [strlabel tmp_str_label] ;
    end
    [c,~,~]= unique(strlabel ,'rows');
    count=size(c,2);
end
end