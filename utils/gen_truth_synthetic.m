% function truth = gen_truth(model)
% this function generate synthetic data for conveyor belt motion model
K = 100 ;
% range_c = [1024 , 768] ; 
range_c = [1 , 1000] ; 
birth_region = [100 900 ; 100 900] ; 
P_on = 0.3 ; 
pSp = 0.03 ; 
init_cells = 10 ;
new_birth_rate = 0.03 ; 
d_s = 10 ; % spawning distance
force_vel_range = [-5 5] ; 
location_std = 3 ; 
spawn_location_std = 0 ; 
vel_std = 5 ; 
pS = 0.99 ; 
overlap_lim = 0 ; 
birth_region_scale = [0.1 , 0.9 ; 0.2 , 0.8] ; 
A = [1 1 ; 0 1] ; 
F = [A , zeros(2) ; zeros(2) , A] ; 
% generate cells for the initial frame
truth.img_size = range_c ; 
truth.X = cell(K,1) ; 
truth.track_list = cell(K,1) ; 
truth.N = zeros(K,1) ; 
truth.spawnNextTimeTrack = cell(K,1) ; 
truth.newSpawn = cell(K,1) ; 
truth.newBirths= cell(K,1) ; 
truth.K = K ; 
for cidx =1 : init_cells
    truth.N(1) = init_cells ;
    truth.track_list{1} = 1:init_cells ; 
    
    truth.X{1} = zeros(4,init_cells) ; 
%     truth.X{1}(1,:) = round(birth_region_scale(1,1)*range_c(1) + ((birth_region_scale(1,2)-birth_region_scale(1,1))*range_c(1)*rand(1,init_cells))) ; 
    truth.X{1}(1,:) = birth_region(1,1) + (birth_region(1,2)-birth_region(1,1))*rand(1,init_cells) ; 
    
    truth.X{1}(2,:) = force_vel_range(1) + (force_vel_range(2)-force_vel_range(1))*rand(1,init_cells) ; 
%     truth.X{1}(3,:) = round(birth_region_scale(2,1)*range_c(2) + ((birth_region_scale(2,2)-birth_region_scale(2,1))*range_c(2)*rand(1,init_cells))) ;
    truth.X{1}(3,:) = birth_region(2,1) + (birth_region(2,2)-birth_region(2,1))*rand(1,init_cells) ; 
    truth.X{1}(4,:) = force_vel_range(1) + (force_vel_range(2)-force_vel_range(1))*rand(1,init_cells) ; 
    truth.newBirths{1} = 1 : init_cells ;
end
curr_max_lab = init_cells ;

for k =2 : K
    % generate new tracks by looping through previous tracks
    truth.N(k) = 0 ; 
    for ot = 1 : length(truth.track_list{k-1})
        deadFlag = true ;
        spawnFlag = false ; 
        if rand<pS % if this track survive
            deadFlag = false ; 
            if rand<pSp % if spawn
                accept = false ; 
                while ~accept
                    accept_1 = false ;
                    accept_2 = false ;
                    spawnFlag = true ; 
                    delta_angle = 180*rand ; 
                    x_shift = d_s * cosd(delta_angle) ; 
                    y_shift = d_s * sind(delta_angle) ;
                    sp1 = zeros(4,1) ; sp2 = zeros(4,1) ; 
                    sp1([1,3]) = round(truth.X{k-1}([1,3],ot) + [x_shift ; y_shift] + spawn_location_std * randn(2,1)) ; 
                    sp1([2,4]) = force_vel_range(1)*ones(2,1) + (force_vel_range(2)-force_vel_range(1))*rand(2,1) ; 
                    sp2([1,3]) = round(truth.X{k-1}([1,3],ot) - [x_shift ; y_shift] + spawn_location_std * randn(2,1)) ;
                    sp2([2,4]) = force_vel_range(1)*ones(2,1) + (force_vel_range(2)-force_vel_range(1))*rand(2,1) ; 
                    if truth.N(k)==0
                        accept_1 = true ; 
                    else
                        all_ok = false(1,truth.N(k)) ; 
                        for nt = 1 : truth.N(k)
                            if norm(sp1([1,3])-truth.X{k}([1,3],nt))<overlap_lim || norm(sp2([1,3])-truth.X{k}([1,3],nt))<overlap_lim 
                                break ; 
                            else
                                all_ok(nt) = true ; 
                            end
                        end
                        accept_1 = all(all_ok) ; 
                    end

                    if sp1(1)>range_c(1) && sp1(1)<range_c(2) && sp1(3)>range_c(1) && sp1(3)<range_c(2)
                        if sp2(1)>range_c(1) && sp2(1)<range_c(2) && sp2(3)>range_c(1) && sp2(3)<range_c(2)
                            accept_2 = true ; 
                        end
                    end
%                     if sp1(1)<1 || sp1(1)>512 || sp1(3)<1 || sp1(3)>512
%                         keyboard
%                     end
%                     if sp2(1)<1 || sp2(1)>512 || sp2(3)<1 || sp2(3)>512
%                         keyboard
%                     end
                    accept = accept_1 & accept_2 ; 
                end

            else % not spawn
                if rand>P_on % take the free diffusion
                    s = truth.X{k-1}(:,ot) ; 
                    s([1,3]) = round(truth.X{k-1}([1,3],ot) + location_std * randn(2,1)) ; 
                else % take the force field
                    s = F*truth.X{k-1}(:,ot); 
%                     s([1,3]) = round(s([1,3]) + location_std * randn(2,1)) ; 
                end
                if s(1)<range_c(1) || s(1)>range_c(2) || s(3)<range_c(1) || s(3)>range_c(2)
                    deadFlag = true ;
                end
            end
        end
        if ~deadFlag
            if spawnFlag
                truth.spawnNextTimeTrack{k-1}=[truth.spawnNextTimeTrack{k-1} , truth.track_list{k-1}(ot)] ; 
                % if parent is outside tracking region I consider these
                % guys as new births not spawns
                if (truth.X{k-1}(1,ot)<range_c(1) || truth.X{k-1}(1,ot)>range_c(2) || ...
                    truth.X{k-1}(3,ot)<range_c(1) || truth.X{k-1}(3,ot)>range_c(2))
                    truth.newBirths{k} = [truth.newBirths{k} , curr_max_lab + 1 , curr_max_lab + 2] ;
                    % remove parent from spawnNextTimeTrack list
                    truth.spawnNextTimeTrack{k-1}(end) = [] ;
                else
                    truth.newSpawn{k} = [truth.newSpawn{k} , curr_max_lab + 1 , curr_max_lab + 2] ;
                end
                truth.X{k} = [truth.X{k} , sp1 , sp2] ; 
                truth.track_list{k} = [truth.track_list{k} , curr_max_lab + 1 , curr_max_lab + 2] ; 
                truth.N(k) = truth.N(k) + 2 ;
                if sp1(1)<1 || sp1(1)>range_c(2) || sp1(3)<1 || sp1(3)>range_c(2)
                    keyboard
                end
                if sp2(1)<1 || sp2(1)>range_c(2) || sp2(3)<1 || sp2(3)>range_c(2)
                    keyboard
                end
                curr_max_lab = curr_max_lab + 2 ; 
            else
%                 truth.X{k} = [truth.X{k} , round(s)] ; 
                truth.X{k} = [truth.X{k} , s] ; 
                if s(1)<1 || s(1)>range_c(2) || s(3)<1 || s(3)>range_c(2)
                    keyboard
                end
                truth.track_list{k} = [truth.track_list{k} , truth.track_list{k-1}(ot)] ; 
                truth.N(k) = truth.N(k) + 1 ; 
            end
        end
    end
    
    % generate new births
    nbirths = poissrnd(new_birth_rate) ; 
    for bt = 1 : nbirths
        accept = false ;
        while ~accept
            X_births = zeros(4,1) ; 
%             X_births(1,:) = round(birth_region_scale(1,1)*range_c(1) + ((birth_region_scale(1,2)-birth_region_scale(1,1))*range_c(1)*rand)) ;
            X_births(1,1) = birth_region(1,1) + (birth_region(1,2)-birth_region(1,1))*rand ; 
%             X_births(3,:) = round(birth_region_scale(2,1)*range_c(2) + ((birth_region_scale(2,2)-birth_region_scale(2,1))*range_c(2)*rand)) ;
            X_births(3,1) = birth_region(2,1) + (birth_region(2,2)-birth_region(2,1))*rand ; 
            X_births([2,4]) = force_vel_range(1)*ones(2,1) + (force_vel_range(2)-force_vel_range(1))*rand(2,1) ;
            all_ok = false(1, truth.N(k)) ; 
            for nt = 1 : truth.N(k) 
                if norm(X_births-truth.X{k}(:,nt)) < overlap_lim
                    break ;
                else
                    all_ok(nt) = true ; 
                end
            end
            accept = all(all_ok) ; 
        end
        truth.X{k} = [truth.X{k} , round(X_births)] ; 
    end
                    
    truth.track_list{k} = [truth.track_list{k} , curr_max_lab+1 : curr_max_lab+nbirths] ; 
    truth.newBirths{k} = [truth.newBirths{k} , curr_max_lab+1 : curr_max_lab+nbirths] ; 
    truth.N(k) = truth.N(k) + nbirths ; 
    curr_max_lab = curr_max_lab+nbirths ;  
%     disp(num2str(k)) ; 
end
truth.total_tracks = curr_max_lab ; 
disp(num2str(truth.N(K))) ; 
% correct_spawn = [] ; 
% for k = 1 : K
%     del_indc = [] ; 
%     for t =1 : truth.N(k)
%         if truth.X{k}(1,t)>range_c(1) || truth.X{k}(2,t)>range_c(2) || any(truth.X{k}(:,t)<0)
%             del_indc = [del_indc,t] ; 
%         end
%     end
%     disp(['k = ',num2str(k),' , delete: ' , num2str(length(del_indc)) , ' tracks']) ; 
%     truth.N(k) = truth.N(k) - length(del_indc) ; 
%     truth.X{k}(:,del_indc) = [] ; 
%     del_track_labs =  truth.track_list{k}(del_indc) ; 
%     [spawns,~,ib] = intersect(del_track_labs , truth.newSpawn{k}) ; 
%     if ~isempty(ib)
%         truth.newSpawn{k}(ib) = [] ; 
%     end
%     [births,~,ib] = intersect(del_track_labs , truth.newBirths{k}) ; 
%     if ~isempty(ib)
%         truth.newBirths{k}(ib) =[] ;
%     end
%     if ~isempty(del_indc)
%         truth.track_list{k}(del_indc) = [] ;
%     end
% end
            
            
            
