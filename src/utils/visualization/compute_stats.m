function [location_intensity, velocity_intensity, card] = compute_stats(setting,  output_path, K)
%% Initialize arrays to store data
pos_heatmap = zeros(setting.dataset_img_size(1,2), setting.dataset_img_size(2,2)) ;
vel_shift = setting.intensity_max_speed ;
vel_heatmap = zeros(2*vel_shift + 1 ,2*vel_shift + 1)  ;

mean_cell_num = zeros(K ,1) ;
std_cell_num = zeros(K ,1) ;

mean_survivals_num = zeros(K ,1) ;
std_survivals_num = zeros(K ,1) ;

mean_births_num = zeros(K ,1) ;
std_births_num = zeros(K ,1) ;

mean_spawns_num = zeros(K ,1) ;
std_spawns_num = zeros(K ,1) ;
%% Compute statistics from posterior GLMB
for k = 2 : K
    % load the glmb posterior
    load( fullfile(output_path, ['glmb_', sprintf('%03d', k), '.mat'])  ) ; 
    % Check the labels
     survival_flag = zeros(length(glmb_update.tt),1) ;
     birth_flag = zeros(length(glmb_update.tt),1) ;
     spawn_flag = zeros(length(glmb_update.tt),1) ;
     for tidx = 1 : length(glmb_update.tt)
         if glmb_update.tt{tidx}.l(end-1) < k
             survival_flag(tidx) = 1 ;
         else
             if length(glmb_update.tt{tidx}.l) > 2
                 spawn_flag(tidx) = 1 ;
             else
                 birth_flag(tidx) = 1 ;
             end
         end
     end
     flag_mat = [survival_flag, birth_flag, spawn_flag] ;
        % count each type of cells
     tt_weights = zeros(length(glmb_update.tt),1) ; % weights of tracks
     cell_num = zeros (length(glmb_update.w) ,1 );
     survivals_num = zeros (length(glmb_update.w) ,1 );
     births_num = zeros (length(glmb_update.w) ,1 );
     spawns_num = zeros (length(glmb_update.w) ,1 );
    for hidx = 1 : length(glmb_update.w)
         tt_weights(glmb_update.I{hidx}) = tt_weights(glmb_update.I{hidx}) + glmb_update.w(hidx) ;
         cell_num(hidx) = length(glmb_update.I{hidx}) ;
         for Iidx = 1 : length(glmb_update.I{hidx})
             flagidx = find(flag_mat(glmb_update.I{hidx}(Iidx),:) == 1) ;
             if flagidx == 1
                 survivals_num(hidx) = survivals_num(hidx) + 1 ;
             elseif flagidx == 2
                 births_num(hidx) = births_num(hidx) + 1;
             elseif flagidx == 3
                 spawns_num(hidx) = spawns_num(hidx) + 1;
             else
                 disp('invalid') ;
                 keyboard
             end
         end
    end
    % Extract cardinality info (note: cdn start from 0 -> max)
    % total cell_num
    cell_num_cdn = zeros(1 , max(cell_num) + 1) ;
    for cdn = 1 : max(cell_num) + 1
        for i = 1 : length(cell_num)
            if cell_num(i) == cdn-1
                cell_num_cdn(cdn) = cell_num_cdn(cdn) + glmb_update.w(i) ;
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
                survivals_num_cdn(cdn) = survivals_num_cdn(cdn) + glmb_update.w(i) ;
            end
        end
    end
    mean_survivals_num(k) = find(survivals_num_cdn==max(survivals_num_cdn)) - 1 ;
    std_survivals_num(k) = std(0:max(survivals_num) , survivals_num_cdn) ;
    %births_num
    births_num_cdn = zeros(1 , max(births_num) + 1) ;
    for cdn = 1 : max(births_num) + 1
        for i = 1 : length(births_num)
            if births_num(i) == cdn-1
                births_num_cdn(cdn) = births_num_cdn(cdn) + glmb_update.w(i) ;
            end
        end
    end
    mean_births_num(k) = find(births_num_cdn==max(births_num_cdn)) - 1 ;
    std_births_num(k) = std(0:max(births_num) , births_num_cdn) ;
    %spawns_num
    spawns_num_cdn = zeros(1 , max(spawns_num) + 1) ;
    for cdn = 1 : max(spawns_num) + 1
        for i = 1 : length(spawns_num)
            if spawns_num(i) == cdn-1
                spawns_num_cdn(cdn) = spawns_num_cdn(cdn) + glmb_update.w(i) ;
            end
        end
    end
    mean_spawns_num(k) = find(spawns_num_cdn==max(spawns_num_cdn)) - 1 ;
    std_spawns_num(k) = std(0:max(spawns_num) , spawns_num_cdn) ;
        
    % Extract the PHD for heatmaps
%     for tidx = 1 : length(glmb_update.tt)
%         x = round(glmb_update.tt{tidx}.m) ;
%         x_pos = x([1 3] , :) ;
%         x_vel = x([2 4] , :) ;
%         if x_pos(1)<=size(pos_heatmap,1) && x_pos(2)<=size(pos_heatmap,2) && x_pos(1)>0 && x_pos(2)>0
%             pos_heatmap(x_pos(1) , x_pos(2)) = pos_heatmap(x_pos(1) , x_pos(2)) + tt_weights(tidx) ;  
%         end
% 
%         if x_vel(1)+vel_shift+1<=size(vel_heatmap,1) && x_vel(2)+vel_shift+1<=size(vel_heatmap,2) && x_vel(1)+vel_shift+1>0 && x_vel(2)+vel_shift+1>0
%             vel_heatmap(x_vel(1)+vel_shift+1 , x_vel(2)+vel_shift+1) = vel_heatmap(x_vel(1)+vel_shift+1 , x_vel(2)+vel_shift+1) + tt_weights(tidx) ;
%         end
%     end  

    for tidx = 1 : length(glmb_update.tt)
        x = round(glmb_update.tt{tidx}.m) ;
        all_x_pos = x([1 3] , :) ;
        all_P_pos = glmb_update.tt{tidx}.P([1 3],[1 3],:) ;
        all_x_vel = x([2 4] , :) ;
        all_P_vel = glmb_update.tt{tidx}.P([2 4],[2 4],:) ;
        for cmpidx = 1 : size(x,2) 
            x_pos = all_x_pos(:,cmpidx) ; 
            P_pos = all_P_pos(:,:,cmpidx) ;
            all_w = glmb_update.tt{tidx}.w ; 
            if ~issymmetric(P_pos)
                disp(P_pos);
                P_pos(2,2) = P_pos(1,1) ; 
            end
            P_pos = (P_pos+P_pos')/2 ; 
            pos_rad = round(crr(P_pos)) ; 
            x_vel = all_x_vel(:,cmpidx) ; 
            P_vel = all_P_vel(:,:,cmpidx) ; 
            P_vel = (P_vel+P_vel')/2 ; 
            vel_rad = round(crr(P_vel)) ; 
            % update the intensity
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
% sum heatmap to bin
step = setting.intensity_location_step ; 
location_intensity = zeros(round(size(pos_heatmap,1)/step),round(size(pos_heatmap,2)/step)) ; 
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
        
        location_intensity(count_row , count_col) = sum(sum(pos_heatmap(i:end_row,j:end_col))) ; 
        count_col = count_col + 1 ;
    end
    count_row = count_row + 1 ;
    count_col = 1 ;
end
location_intensity = flipud(location_intensity)/K ; 
velocity_intensity = vel_heatmap ; 
card = cat(2, mean_cell_num, std_cell_num, mean_survivals_num, std_survivals_num, ...
    mean_births_num, std_births_num, mean_spawns_num,  std_spawns_num);

     
% % Make plots
% % Plot the heatmap
% figure
% heatmap(location_intensity,'Colormap',jet) ; 
% set(gca,'FontSize', 22) ; 
% 
% figure
% 
% vel_heat_axis = linspace(-50, 50, 101) ; 
% h = heatmap(vel_heat_axis,vel_heat_axis,vel_heatmap,'Colormap',jet) ;
% for i = 1 : length(vel_heat_axis)
%     if i == 1 || mod(i-1,10) == 0 
%         vel_tick_label{i} = num2str(vel_heat_axis(i)) ;
%     else
%         vel_tick_label{i} = '' ; 
%     end
% end
% set(gca,'XDisplayLabels',vel_tick_label) ; 
% set(gca,'YDisplayLabels',vel_tick_label) ; 
% set(gca,'FontSize', 18) ; 
% % plot velocity
% x_vel = sum(vel_heatmap,1)/88 ;
% y_vel = sum(vel_heatmap,2)/88 ;
% LW = 2 ;
% FS = 18 ; 
% figure
% vx = plot(-vel_shift:vel_shift, x_vel, 'Color', 'b', 'LineWidth',LW) ; hold on;
% vy = plot(-vel_shift:vel_shift, y_vel, 'Color', 'r', 'LineWidth',LW) ;
% hold off
% grid on
% legend([vx,vy], {'x-velocity', 'y-velocity'}, 'FontSize', FS)
% set(gca,'FontSize',18) ; 
% 
% 
% 
% FS = 20 ; 
% figure; cardinality= gcf; 
% subplot(2,2,1); box on; grid on ; hold on;
% plot(1:K,mean_cell_num,'LineWidth',2);
% plot(1:K,mean_cell_num+std_cell_num,'k--');
% plot(1:K,mean_cell_num-std_cell_num,'k--');
% set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 80],'FontSize', 14); ylabel('Total','FontSize', FS); 
% title('Cells cardinality statistics','FontSize', 16);
% legend({'Mean values' , '3 sigma bound'},'FontSize', FS) ;
% set(gca, 'FontSize', FS);
% 
% subplot(2,2,2); box on; grid on ; hold on;
% plot(1:K,mean_survivals_num,'LineWidth',2);
% plot(1:K,mean_survivals_num+std_survivals_num,'k--');
% plot(1:K,mean_survivals_num-std_survivals_num,'k--');
% set(gca, 'XLim',[1 K]); set(gca, 'YLim',[0 80],'FontSize', 14); ylabel('Survivals','FontSize', FS);
% set(gca, 'FontSize', FS);
% 
% subplot(2,2,3);box on; grid on ; hold on;
% plot(1:K,mean_births_num,'LineWidth',2);
% plot(1:K,mean_births_num+std_births_num,'k--');
% plot(1:K,mean_births_num-std_births_num,'k--');
% set(gca, 'XLim',[1 K]); set(gca, 'YLim',[-5 10],'FontSize', 14); ylabel('Instantaneous births','FontSize', FS);
% set(gca, 'FontSize', FS);
% 
% subplot(2,2,4);box on; grid on ; hold on;
% plot(1:K,mean_spawns_num/2,'LineWidth',2);
% plot(1:K,mean_spawns_num/2+std_spawns_num,'k--');
% plot(1:K,mean_spawns_num/2-std_spawns_num,'k--');
% set(gca, 'XLim',[1 K]); set(gca, 'YLim',[-3 10],'FontSize', 14); ylabel('Mitotic cells','FontSize', FS);
% 
% set(gca, 'FontSize', FS); 
% xlabel('Time step','FontSize', FS);
% % Analysing est
% % for k = 2 : length(glmb_out)
% %     [est.X{k},est.N(k),est.L{k}]= extract_estimates(glmb_out(k),model);
% % end
% est = est_ah ; 
% x_vel = [] ;
% for k= 2: K
%     for t = 1 : size(est.X{k},2)
%         x = est.X{k}(:,t) ;
%         x_vel = [x_vel x([2 4],1)] ;
%     end
% end
% x_vel = x_vel' ;
% % figure
% % boxplot(x_vel , 'Labels',{'x-axis','y-axis'},'Whisker',1,'FontSize', 14) ;
% % title('Velocity') ;
% % Plot cell tracks
% % Grab estimated x,y pos
% % pos_xy = [1 3];
% % xy_pos = cell( length(est.X) , 1 );
% % for k = 1 : length(est.X)
% %     if ~isempty( est.X{k} )
% %         xy_pos{k} = est.X{k}(pos_xy , :) ;
% %     end
% % end 
% % backgroundImg = cell(1,88) ; 
% % for k = 1 : K
% %     path = ['data/Yu_Suk_Experiments/2017-07-06_M33/', sprintf('%03d', k-1), '.jpg'] ; 
% %     backgroundImg{k} = imread(path) ; 
% % end
%     
% % [tracks_img , trans_tracks_img , transdis] = analyse_cell_tracks(xy_pos , est.L(1:88) , 88, backgroundImg);
% % figure
% % imshow(tracks_img.cdata) ;
% % title('Trajectories');
% % 
% % figure
% % imshow(trans_tracks_img.cdata) ;
% % title('Translated trajectories');
% % 
% % figure
% % h = boxplot(transdis , 'Labels',{'x-axis','y-axis'},'Whisker',1) ;
% % set(h,{'linew'},{3});
% % set(gca, 'FontSize', 30);
% 
% 
% % title('Final displacement','FontSize', 18) ;
% % make the ancestry plots
% % read the ancestry info in
% unique_labs = [] ; 
% for k = 1 : length(est.L)
%     for lidx = 1 : length(est.L{k})
%         unique_labs = [unique_labs , {sprintf('%.0f-' , est.L{k}{lidx})}] ;
%         unique_labs = unique(unique_labs) ; 
%     end
% end
% nodes = form_nodes(unique_labs) ; 
% % clean nodes before tree plot .. 
% % nodes_plot = cell(1,0) ;
% % del_indc = [] ; 
% % nodes_temp = nodes ; 
% % for ii = 1 : length(nodes)-1
% %     if nodes(ii) == 0 && nodes(ii+1) == 0
% %         del_indc = [del_indc , ii+1] ; 
% %     end
% % end
% % nodes(del_indc) = [] ; 
% % if nodes(end) == 0 
% %     nodes(end) =[] ; 
% % end
% my_treeplot({nodes}, [0 0 0]) ;
% 
% % generate a tree for tree plot
% % find fami
% fam_ances = zeros(2,0) ; 
% fam = cell(1,0) ;
% % init_ances = cell(1,0) ; 
% for k = 1 : length(est.L)
%     for t = 1 : length(est.L{k})
%         lab = est.L{k}{t} ; 
%         ances_lab = lab([1:2],1) ; 
%         curr_fam_idx = find_vector_in_matrix(ances_lab, fam_ances) ; 
%         
%         if curr_fam_idx == 0
% %             init_ances(:,fam_idx) = ances_lab ; 
%             fam_ances = [fam_ances, ances_lab] ; 
%             fam{length(fam)+1} = cell(1,0) ; 
%         else
%             temp_idx = find_vector_in_cell(lab, fam{curr_fam_idx}) ; 
%             if temp_idx == 0
%                 fam{curr_fam_idx}{length(fam{curr_fam_idx})+1} = lab ; 
%             end
%         end
%     end
% end
% 
% nodes_fam = cell(1,length(fam)) ; 
% for fidx = 1 : length(fam)
%     fam_temp = fam{fidx} ; 
%     for n = 1 : length(fam{fidx})
%         fam_temp{n} = sprintf('%.0f-' , fam_temp{n});
%     end
%     nodes_fam{fidx} = form_nodes(fam_temp) ; 
% end
% colord = distinguishable_colors(length(nodes_fam)) ; 
% my_treeplot(nodes_fam, colord) ; 
% 
% 
% function nodes = form_nodes(unique_labs)
% nodes = zeros(1, length(unique_labs)) ;
% if length(unique_labs)>1
%     for lidx = 1 : length(unique_labs)
%         % split the strings
%         temp = unique_labs{lidx} ; 
%         temp(end) = [] ; 
%         temp = split(temp,'-') ; 
%         if length(temp) == 2
%             nodes(lidx) = 0 ;
%         else
%             L = length(temp) - 2 ;
%             target_str =[] ; 
%             for ii = 1 : L
%                 target_str = [target_str, temp{ii}, '-'] ; 
%             end
%             cnt = 1 ; 
%             not_found = 1 ;
%             while cnt<=length(unique_labs) && not_found
%                 if strcmp(target_str, unique_labs(cnt))
%                     not_found = 0 ; 
%                     nodes(lidx) = cnt ; 
%                 else
%                     cnt = cnt + 1 ;
%                 end
%             end
%         end
%     end
% end
% end
% 
% function idx = find_vector_in_matrix(v,m)
% idx = 0 ;
% L = size(m,2) ;
% ii = 1 ; 
% while ii<=L && idx==0
%     if all(v-m(:,ii) == 0)
%         idx = ii ; 
%     end
%     ii = ii + 1 ;
% end
% end
% 
% function idx = find_vector_in_cell(v,c)
% idx = 0 ; 
% L = length(c) ; 
% ii = 1 ; 
% while ii<=L && idx==0
%     if length(v) == length(c{ii})
%         if all(v-c{ii}==0)
%             idx = ii ; 
%         end
%     end
%     ii = ii + 1 ; 
% end
% end
%     
        
            


        
 

 

     
     