function plot_truth(model, truth)
    %% Plots properties
    LW = 2 ; % linewidth
    FS = 14 ; % font size
    marker_size = 80 ;
    %% extract the tracks (grab 2-D matrices for OSPA computation)
    [X_track,k_birth,k_death,truth_labels]= extract_truth_tracks(truth);
    truth_colorarray= distinguishable_colors(length(truth_labels));
    %% plot ground truth
    spawn_id = [] ; 
    for k = 1 : truth.K
        spawn_id = [spawn_id, truth.newSpawn{k}] ; 
    end
       % plot the figure
    figure('Name', 'Ground truth')
    for t = 1 : length(truth_labels)
        true_lab= truth_labels(t) ; 
        Pt = X_track(model.xy_pos, k_birth(t):k_death(t), t) ; 
        plot(Pt(1,:), Pt(2,:), 'LineWidth', LW, 'Color',  truth_colorarray(t, :), 'LineWidth', LW) ; 
        hold on ; 
        if ismember(true_lab,spawn_id)
            marker_style = '^' ; 
        else
            marker_style = 's' ; 
        end
        s = scatter(Pt(1,1),  Pt(2,1), marker_size, marker_style, 'MarkerFaceColor', truth_colorarray(t, :), ...
                'MarkerEdgeColor', [0 0 0]) ;
        if strcmp(marker_style,'^')
            p1 = s ; % spawn
        else
            p2 = s ; % birth
        end
        d = scatter(Pt(1,end),  Pt(2,end), 'd', marker_style, 'MarkerFaceColor', truth_colorarray(t, :), ...
                'MarkerEdgeColor', [0 0 0]) ;
        p3 = d ; % dead case
    end
    legend([p2 p1 p3], {'Birth', 'Mitosis', 'Dead'}, 'FontSize', FS, 'NumColumns',3) ;  
    xlim(model.range_c(1,:)) ; ylim(model.range_c(2,:)) ; 
    xlabel('x coordinate (pixel)','FontSize', FS) ; ylabel('y coordinate (pixel)','FontSize', FS) ; 
    grid on
    box on
end


function [X_track, k_birth, k_death, labels]= extract_truth_tracks(truth)
    K= length(truth.X); 
    x_dim= size(truth.X{K},1); 
    k=K-1; while x_dim==0, x_dim= size(truth.X{k},1); k= k-1; end

    labels = [] ; 
    for k = 1 : K
        labels = [labels, truth.track_list{k}] ; 
        labels = unique(labels) ; 
    end
    num_tracks = length(labels) ; 
    X_track= NaN(x_dim,K, num_tracks);
    k_birth = nan(1, num_tracks) ; 
    k_death = zeros(1, num_tracks) ; 

    for k=1:K
        if ~isempty(truth.X{k})
            for lidx = 1 : length(truth.track_list{k})
                [~, trkidx] = ismember(truth.track_list{k}(lidx), labels) ; 
                X_track(:,k,trkidx)= truth.X{k}(:,lidx);
                if isnan(k_birth(trkidx))
                    k_birth(trkidx) = k ; 
                end
                k_death(trkidx) = k ; 
            end
        end
    end
end
