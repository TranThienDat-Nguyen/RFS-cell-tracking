function [ospa_vals , ospa2_cell] = eval_results(model, truth, meas, est, cutoff, order, win_len, identifier, output_path)
    %% Extract the tracks (grab 2-D matrices for OSPA computation)
    [X_track,~,~,~]= extract_truth_tracks(truth);
    [Y_track,~,~,~]= extract_est_tracks(est);
    %% compute OSPA and OSPA(2)
        % compute OSPA
    disp('Computing OSPA distance ...')
    ospa_vals= zeros(truth.K,3);
    for k=1:truth.K
        [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},model.xy_pos),get_comps(est.X{k},model.xy_pos),cutoff,order);
    end
        % compute OSPA^(2)
    win_len_comp= [win_len,truth.K];
    ospa2_cell = cell(1,length(win_len_comp));
    for i = 1:length(win_len_comp)
        ospa2_cell{i} = compute_ospa2(X_track(model.xy_pos,:,:),Y_track(model.xy_pos,:,:),cutoff,order,win_len_comp(i));
    end
        % compute TRA
    disp('Computing TRA score ...')
    % convert from MATLAB structures to CTC format
    truth_struct2CTC(truth, output_path, model.range_c(1,2), model.range_c(2,2)) ; 
    est_struct2CTC(est, output_path, model.range_c(1,2), model.range_c(2,2)) ;     
    gtPath = fullfile(output_path, 'truth') ; 
    resPath = fullfile(output_path, 'res') ; 
    [TRA_score, ~, ~, ~] = PerformanceTRA(gtPath, resPath, truth.K, cutoff) ; % note: -1 since CTC convention time start at 0

    %% Plot ospa results
    figure('Name', ['OSPA_resutls_', identifier]); hold on;
    subplot(3,1,1); plot(1:meas.K,ospa_vals(:,1),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA Dist');
    subplot(3,1,2); plot(1:meas.K,ospa_vals(:,2),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA Loc');
    subplot(3,1,3); plot(1:meas.K,ospa_vals(:,3),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA Card');
    xlabel('Time');

    figure('Name', ['OSPA2_resutls_', identifier]); hold on;
    subplot(3,1,1);
    plot(1:truth.K,ospa2_cell{1}(1,:),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA(2) Dist');
    subplot(3,1,2);
    plot(1:truth.K,ospa2_cell{1}(2,:),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA(2) Loc');
    subplot(3,1,3);
    plot(1:truth.K,ospa2_cell{1}(3,:),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 cutoff]); ylabel('OSPA(2) Card');
    xlabel('Time');

    %% Tabulate the results
    tab_var_names = {'OSPA(2) Overall', 'OSPA(2) Localization', 'OSPA(2) Cardinality', 'TRA Score'} ; 
    tab_ospa2 = table(ospa2_cell{2}(1,end), ...
        ospa2_cell{2}(2,end), ospa2_cell{2}(3,end), TRA_score, 'VariableNames', tab_var_names) ; 
    tab_ospa2 = table(tab_ospa2, 'VariableNames', {['Overall OSPA(2) error and TRA score: ', identifier]} ) ; 
    disp(tab_ospa2) ; 

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

function [Y_track, k_birth, k_death, labels]= extract_est_tracks(est)
    K= length(est.X); 
    x_dim= size(est.X{K},1); 
    k=K-1; while x_dim==0, x_dim= size(est.X{k},1); k= k-1; end
    labels = [] ; 
    for k = 1 : K
        temp_labs = cell(1, length(est.L{k})) ; 
        for ll = 1 : length(est.L{k})
            temp_labs{ll} = sprintf('%.0f,' , est.L{k}{ll});
        end
        labels = [labels, temp_labs] ; 
        labels = unique(labels) ; 
    end
    num_tracks = length(labels) ; 
    Y_track= NaN(x_dim,K, num_tracks);
    k_birth = nan(1, num_tracks) ; 
    k_death = zeros(1, num_tracks) ; 

    for k=1:K
        if ~isempty(est.X{k})
            for lidx = 1 : length(est.L{k})
                [~, trkidx] = ismember(sprintf('%.0f,' , est.L{k}{lidx}), labels) ; 
                Y_track(:,k,trkidx)= est.X{k}(:,lidx);
                if isnan(k_birth(trkidx))
                    k_birth(trkidx) = k ; 
                end
                k_death(trkidx) = k ; 
            end
        end
    end
end


function Xc= get_comps(X,c)
    if isempty(X)
        Xc= [];
    else
        Xc= X(c,:);
    end
end