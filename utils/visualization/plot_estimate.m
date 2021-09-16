function plot_estimate(model, truth, meas, est, identifier)
    %% extract the tracks (grab 2-D matrices for OSPA computation)
    [X_track,k_birth,k_death,truth_labels]= extract_truth_tracks(truth);
    [Y_track,l_birth,l_death,est_labels]= extract_est_tracks(est);
    % generate unique colors
    est_colorarray= distinguishable_colors(length(est_labels));
    %% Plot measurements and tracking results
    figure('Name', ['Track_estimates_', identifier]); hold on;

    %plot x measurement
    subplot(211); box on; 

    for k=1:meas.K
        if ~isempty(meas.Z{k})
            hlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(1,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
        end   
    end

    %plot x track
    for i=1:length(truth_labels)
        Px= X_track(:,k_birth(i):1:k_death(i),i); Px=Px(model.xy_pos,:);
        hline1= line(k_birth(i):1:k_death(i),Px(1,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
    end

    %plot x estimate
    for t=1:size(Y_track,3)
        hline2= line(1:truth.K,Y_track(1,:,t),'LineStyle','none','Marker','.','Markersize',8,'Color',est_colorarray(t,:));
    end

    %plot y measurement
    subplot(212); box on;

    for k=1:meas.K
        if ~isempty(meas.Z{k})
            yhlined= line(k*ones(size(meas.Z{k},2),1),meas.Z{k}(2,:),'LineStyle','none','Marker','x','Markersize',5,'Color',0.7*ones(1,3));
        end
    end

    %plot y track
    for i=1:truth.total_tracks
            Py= X_track(:,k_birth(i):1:k_death(i),i); Py=Py(model.xy_pos,:);
            yhline1= line(k_birth(i):1:k_death(i),Py(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
    end

    %plot y estimate
    for t=1:size(Y_track,3)
        hline2= line(1:truth.K,Y_track(3,:,t),'LineStyle','none','Marker','.','Markersize',8,'Color',est_colorarray(t,:));
    end

    subplot(211); xlabel('Time'); ylabel('x-coordinate (m)');
    set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(1,:));
    legend([hline2 hline1 hlined],'Estimates          ','True tracks','Measurements');

    subplot(212); xlabel('Time'); ylabel('y-coordinate (m)');
    set(gca, 'XLim',[1 truth.K]); set(gca, 'YLim',model.range_c(2,:));


    %% Plot cardinality
    figure('Name', ['Cardinality_', identifier]); 
    box on; hold on;
    stairs(1:meas.K,truth.N,'k'); 
    plot(1:meas.K,est.N,'k.');

    grid on;
    legend(gca,'True','Estimated');
    set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 max([max(truth.N), max(est.N)])+1]);
    xlabel('Time'); ylabel('Cardinality');


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