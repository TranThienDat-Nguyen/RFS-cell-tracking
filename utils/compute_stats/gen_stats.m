function gen_stats(setting, output_path, id)
%% Load glmb density and compute corresponding statistics
    disp('Generating statistics of tracking results ...')
    file_info = dir(fullfile(output_path, 'glmb', '*.mat')) ; 
    K = length(file_info) ; 
    % Initialize heatmaps
    pos_heatmap = zeros(setting.dataset_img_size(1,2),setting.dataset_img_size(2,2)) ;
    vel_shift = setting.intensity_max_speed ;
    vel_heatmap = zeros(2*vel_shift + 1 ,2*vel_shift + 1)  ;
    % Initialize cardinality
    mean_cell_num = zeros(K ,1) ;
    std_cell_num = zeros(K ,1) ;
    mean_survivals_num = zeros(K ,1) ;
    std_survivals_num = zeros(K ,1) ;
    mean_births_num = zeros(K ,1) ;
    std_births_num = zeros(K ,1) ;
    mean_spawns_num = zeros(K ,1) ;
    std_spawns_num = zeros(K ,1) ;
    % Compute cardinality stats
    for k = 2 : K  % only start from 2 since we are using measurement driven birth
        temp_load = load(fullfile(output_path, 'glmb', [sprintf('%03d',k), '.mat']));
        glmb = temp_load.glmb ; 
        % Check the labels
        survival_flag = zeros(length(glmb.tt),1) ;
        birth_flag = zeros(length(glmb.tt),1) ;
        spawn_flag = zeros(length(glmb.tt),1) ;
        for tidx = 1 : length(glmb.tt)
            if (glmb.tt{tidx}.l(end-1) == (k-1)) && (length(glmb.tt{tidx}.l) == 2)
                birth_flag(tidx) = 1 ;
            elseif glmb.tt{tidx}.l(end-1) == k
                spawn_flag(tidx) = 1 ;
            else
                survival_flag(tidx) = 1 ;
            end
         end
        flag_mat = [survival_flag birth_flag spawn_flag] ;

        % Loop through hypotheses to get info
        tt_weights = zeros(length(glmb.tt),1) ; % weights of tracks
        cell_num = zeros (length(glmb.w) ,1 );
        survivals_num = zeros (length(glmb.w) ,1 );
        births_num = zeros (length(glmb.w) ,1 );
        spawns_num = zeros (length(glmb.w) ,1 );
        for hidx = 1 : length(glmb.w)
             tt_weights(glmb.I{hidx}) = tt_weights(glmb.I{hidx}) + glmb.w(hidx) ;
             cell_num(hidx) = length(glmb.I{hidx}) ;
             for Iidx = 1 : length(glmb.I{hidx})
                 flagidx = find(flag_mat(glmb.I{hidx}(Iidx),:) == 1) ;
                 if flagidx == 1
                     survivals_num(hidx) = survivals_num(hidx) + 1 ;
                 elseif flagidx == 2
                     births_num(hidx) = births_num(hidx) + 1;
                 elseif flagidx == 3
                     spawns_num(hidx) = spawns_num(hidx) + 1;
                 end
             end

        end
        % Extract card statistic (note: cdn start from 0 -> max)
            % cell_num
        cell_num_cdn = zeros(1 , max(cell_num) + 1) ;
        for cdn = 1 : max(cell_num) + 1
            for i = 1 : length(cell_num)
                if cell_num(i) == cdn-1
                    cell_num_cdn(cdn) = cell_num_cdn(cdn) + glmb.w(i) ;
                end
            end
        end
        mean_cell_num(k) = find(cell_num_cdn==max(cell_num_cdn)) - 1 ;
        std_cell_num(k) = std(0:max(cell_num) , cell_num_cdn) ;
            % survivals_num
        survivals_num_cdn = zeros(1 , max(survivals_num) + 1) ;
        for cdn = 1 : max(survivals_num) + 1
            for i = 1 : length(survivals_num)
                if survivals_num(i) == cdn-1
                    survivals_num_cdn(cdn) = survivals_num_cdn(cdn) + glmb.w(i) ;
                end
            end
        end
        mean_survivals_num(k) = find(survivals_num_cdn==max(survivals_num_cdn)) - 1 ;
        std_survivals_num(k) = std(0:max(survivals_num) , survivals_num_cdn) ;
            % births_num
        births_num_cdn = zeros(1 , max(births_num) + 1) ;
        for cdn = 1 : max(births_num) + 1
            for i = 1 : length(births_num)
                if births_num(i) == cdn-1
                    births_num_cdn(cdn) = births_num_cdn(cdn) + glmb.w(i) ;
                end
            end
        end
        mean_births_num(k) = find(births_num_cdn==max(births_num_cdn)) - 1 ;
        std_births_num(k) = std(0:max(births_num) , births_num_cdn) ;
            % spawns_num
        spawns_num_cdn = zeros(1 , max(spawns_num) + 1) ;
        for cdn = 1 : max(spawns_num) + 1
            for i = 1 : length(spawns_num)
                if spawns_num(i) == cdn-1
                    spawns_num_cdn(cdn) = spawns_num_cdn(cdn) + glmb.w(i) ;
                end
            end
        end
        mean_spawns_num(k) = find(spawns_num_cdn==max(spawns_num_cdn)) - 1 ;
        std_spawns_num(k) = std(0:max(spawns_num) , spawns_num_cdn) ;

        % compute values of heatmaps (use the Gaussians)
        for tidx = 1 : length(glmb.tt)
            x = round(glmb.tt{tidx}.m) ;
            all_x_pos = x([1 3] , :) ;
            all_P_pos = glmb.tt{tidx}.P([1 3],[1 3],:) ;
            all_x_vel = x([2 4] , :) ;
            all_P_vel = glmb.tt{tidx}.P([2 4],[2 4],:) ;
            for cmpidx = 1 : size(x,2) 
                % position intensity
                x_pos = all_x_pos(:,cmpidx) ; 
                P_pos = all_P_pos(:,:,cmpidx) ;
                all_w = glmb.tt{tidx}.w ; 
                P_pos = (P_pos+P_pos')/2 ; 
                pos_rad = round(crr(P_pos)) ;
                % velocity intensity
                x_vel = all_x_vel(:,cmpidx) ; 
                P_vel = all_P_vel(:,:,cmpidx) ; 
                P_vel = (P_vel+P_vel')/2 ; 
                vel_rad = round(crr(P_vel)) ; 
                % sum values to the heatmaps
                    % position
                i_min = max(1,x_pos(1)-pos_rad) ; i_max = min(size(pos_heatmap,1),x_pos(1)+pos_rad) ; 
                j_min = max(1,x_pos(2)-pos_rad) ; j_max = min(size(pos_heatmap,2),x_pos(2)+pos_rad) ; 
                for i = i_min : i_max
                    for j = j_min : j_max
                        pos_heatmap(i,j) = pos_heatmap(i,j) + mvnpdf([i;j] , x_pos , P_pos)*tt_weights(tidx)*all_w(cmpidx); 
                    end
                end
                    % velocity
                i_min = max(1,x_vel(1)-vel_rad+vel_shift+1) ; i_max = min(size(vel_heatmap,1),x_vel(1)+vel_rad+vel_shift+1) ; 
                j_min = max(1,x_vel(2)-vel_rad+vel_shift+1) ; j_max = min(size(vel_heatmap,2),x_vel(2)+vel_rad+vel_shift+1) ; 
                for i = i_min : i_max
                    for j = j_min : j_max
                        vel_heatmap(i,j) = vel_heatmap(i,j) + mvnpdf([i-(vel_shift+1);j-(vel_shift+1)] , x_vel , P_vel)*tt_weights(tidx)*all_w(cmpidx); 
                    end
                end
            end
        end 
    end
    % decrease the resolution of position heatmap
    step = setting.intensity_location_step ; 
    new_pos_heatmap = zeros(round(size(pos_heatmap,1)/step),round(size(pos_heatmap,2)/step)) ; 
    count_row = 1 ;
    count_col = 1 ;
    for i = 1 : step : size(pos_heatmap,1)
        for j = 1 : step : size(pos_heatmap,2)
            if i+step-1 > size(pos_heatmap,1)
                end_row = size(pos_heatmap,1) ; 
            else
                end_row = i+step-1 ; 
            end

            if j+step-1 > size(pos_heatmap,2)
                end_col = size(pos_heatmap,2) ; 
            else
                end_col = j+step-1 ; 
            end

            new_pos_heatmap(count_row , count_col) = sum(sum(pos_heatmap(i:end_row,j:end_col))) ; 
            count_col = count_col + 1 ;
        end
        count_row = count_row + 1 ;
        count_col = 1 ;
    end
    
    %% Plotting the results
        % Plots properties
    FS = 12 ;
    fs = 12 ; 
        % Plot the cardinality
    figure('Name', ['Cardinality_statistics_', id]); 
    subplot(2,2,1); box on; grid on ; hold on;
    plot(1:K,mean_cell_num,'LineWidth',2);
    plot(1:K,mean_cell_num+std_cell_num,'k--');
    plot(1:K,mean_cell_num-std_cell_num,'k--');
    set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 round(1.1*max(mean_cell_num+std_cell_num))],'FontSize', fs); ylabel('Total','FontSize', FS); 
    title('Cells cardinality statistics','FontSize', 16);
    legend({'Mean values' , '3 sigma bound'},'FontSize', FS) ;
    set(gca, 'FontSize', FS);

    subplot(2,2,2); box on; grid on ; hold on;
    plot(1:K,mean_survivals_num,'LineWidth',2);
    plot(1:K,mean_survivals_num+std_survivals_num,'k--');
    plot(1:K,mean_survivals_num-std_survivals_num,'k--');
    set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 round(1.1*max(mean_survivals_num+std_survivals_num))],'FontSize', fs); ylabel('Survivals','FontSize', FS);
    set(gca, 'FontSize', FS);

    subplot(2,2,3);box on; grid on ; hold on;
    plot(1:K,mean_births_num,'LineWidth',2);
    plot(1:K,mean_births_num+std_births_num,'k--');
    plot(1:K,mean_births_num-std_births_num,'k--');
    set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 round(1.1*max(mean_births_num+std_births_num))],'FontSize', fs); ylabel('Instantaneous births','FontSize', FS);
    set(gca, 'FontSize', FS);

    subplot(2,2,4);box on; grid on ; hold on;
    plot(1:K,mean_spawns_num/2,'LineWidth',2);
    plot(1:K,mean_spawns_num/2+std_spawns_num,'k--');
    plot(1:K,mean_spawns_num/2-std_spawns_num,'k--');
    set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 round(1.1*max(mean_spawns_num/2+std_spawns_num))],'FontSize', fs); ylabel('Mitotic cells','FontSize', FS);
    set(gca, 'FontSize', FS); 

    xlabel('Time step','FontSize', FS);
    title('Cardinality statistics')
    
        % Plot the position heatmaps
    figure('Name', ['Position_intensity_', id])
    heatmap(new_pos_heatmap,'Colormap',jet) ; 
    set(gca,'FontSize', fs) ; 
    temp_ax = gca ; 
    temp_ax.XDisplayLabels = nan(size(temp_ax.XDisplayData));
    temp_ax.YDisplayLabels = nan(size(temp_ax.YDisplayData));
    title('Cell intensity in position space')
        % Plot the velocity heatmap
    figure('Name', ['Velocity_intensity_', id])
    vel_heat_axis = linspace(-vel_shift, vel_shift, vel_shift*2+1) ; 
    heatmap(vel_heat_axis,vel_heat_axis, flipud(vel_heatmap), 'Colormap', jet) ; % flip the velocity heatmap
    for i = 1 : length(vel_heat_axis)
        if i == 1 || mod(i-1,10) == 0 
            vel_tick_label{i} = num2str(vel_heat_axis(i)) ;
        else
            vel_tick_label{i} = '' ; 
        end
    end
    set(gca,'XDisplayLabels',vel_tick_label) ; 
    set(gca,'YDisplayLabels',vel_tick_label(end:-1:1)) ;  % flip the y-axis labels to match the flipped velocity heatmap
    set(gca,'FontSize', fs) ; 
    title('Cell intensity in velocity space')
end   