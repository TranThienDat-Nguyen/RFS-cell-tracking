% function truth = gen_truth(model)
% this function generate synthetic data for conveyor belt motion model
K = 100 ;
% range_c = [1024 , 768] ; 
range_c = [1 , 1000] ; 
birth_region = [100 900 ; 100 900] ; 
P_on = 0.9 ; 
pSp = 0.03 ; 
init_cells = 5 ;
new_birth_rate = 0.03 ; 
d_s = 10 ; % spawning distance
force_vel_range = [-5 5] ; 
location_std = 3 ; 
spawn_location_std = 0 ; 
vel_std = 5 ; 
pS = 0.95 ; 
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

% Initial birth
nbirths = 21;
% we assume that daughter cells have same velocity

parent(:,1) = [0;0] ; 
xstart(:,1)  = [ 110; 8; 100; 2.5 ];  
tbirth(1)  = 1; 
tdead(1) = 9 ; 

parent(:,2) = [0;0] ; 
xstart(:,2)  = [ 100; 10; 505; -6 ];  
tbirth(2)  = 1; 
tdead(2)  = 14;

parent(:,3) = [0;0] ; 
xstart(:,3)   = [ 890; -7.5; 495; -4.5 ]; 
tbirth(3)   = 1; 
tdead(3)  = 39;

parent(:,4) = [1;1] ; 
xstart(:,4) = [0 ; 0 ; 0 ; 0] ; % just padding 
tbirth(4) = 10 ; 
tdead(4) =  19   ;

parent(:,5) = [1;2] ; 
xstart(:,5) = [0 ; 0 ; 0 ; 0] ; % just padding 
tbirth(5) = 10 ; 
tdead(5) =  70   ;

parent(:,6) = [2;1] ; 
xstart(:,6) = [0 ; 0 ; 0 ; 0] ; 
tbirth(6) = 15 ; 
tdead(6) =  39 ; 

parent(:,7) = [2;2] ; 
xstart(:,7) = [0 ; 0 ; 0 ; 0] ; 
tbirth(7) = 15 ; 
tdead(7) =  80 ; 

parent(:,8) = [0;0] ; 
xstart(:,8) = [ 893; -8.5; 907; -1 ] ; 
tbirth(8) = 20 ; 
tdead(8) =  74 ; 

parent(:,9) = [4;1] ; 
xstart(:,9) = [0 ; 0 ; 0 ; 0] ; 
tbirth(9) = 20 ; 
tdead(9) =  60 ; 

parent(:,10) = [4;2] ; 
xstart(:,10) = [0 ; 0 ; 0 ; 0] ; 
tbirth(10) = 20 ; 
tdead(10) =  100 ; 

parent(:,11) = [0;0] ; 
xstart(:,11) = [ 908; -8.5; 106; 0 ]; 
tbirth(11) = 30 ; 
tdead(11) =  100 ; 

parent(:,12) = [0;0] ; 
xstart(:,12) = [ 97; 6; 909; -7 ]; 
tbirth(12) = 30 ; 
tdead(12) =  100 ; 

parent(:,13) = [0;0] ; 
xstart(:,13) = [ 100; 8; 100; 8.5]; 
tbirth(13) = 30 ; 
tdead(13) =  84 ; 

parent(:,14) = [3;1] ; 
xstart(:,14) = [0 ; 0 ; 0 ; 0] ; 
tbirth(14) = 40 ; 
tdead(14) =  70 ; 

parent(:,15) = [3;2] ; 
xstart(:,15) = [0 ; 0 ; 0 ; 0] ; 
tbirth(15) = 40 ; 
tdead(15) =  60 ; 

parent(:,16) = [6;1] ; 
xstart(:,16) = [0 ; 0 ; 0 ; 0] ; 
tbirth(16) = 40 ; 
tdead(16) =  80 ; 

parent(:,17) = [6;2] ; 
xstart(:,17) = [0 ; 0 ; 0 ; 0] ; 
tbirth(17) = 40 ; 
tdead(17) =  80 ; 


parent(:,18) = [8;1] ; 
xstart(:,18) = [0 ; 0 ; 0 ; 0] ; 
tbirth(18) = 75 ; 
tdead(18) =  100 ; 

parent(:,19) = [8;2] ; 
xstart(:,19) = [0 ; 0 ; 0 ; 0] ; 
tbirth(19) = 75 ; 
tdead(19) =  100 ; 

parent(:,20) = [13;1] ; 
xstart(:,20) = [0 ; 0 ; 0 ; 0] ; 
tbirth(20) = 85 ; 
tdead(20) =  100 ; 

parent(:,21) = [13;2] ; 
xstart(:,21) = [0 ; 0 ; 0 ; 0] ; 
tbirth(21) = 85 ; 
tdead(21) =  100 ; 




truth.total_tracks = length(tdead); 
for ii = 1 : truth.total_tracks
    targetstate = xstart(:,ii);
    if parent(1,ii) == 0 
        truth.newBirths{tbirth(ii)} = [truth.newBirths{tbirth(ii)} , ii] ; 
        % propagate normally
        for k = tbirth(ii) : min(tdead(ii),truth.K)
            if rand>P_on
                targetstate([1,3]) = targetstate([1,3]) + location_std * randn(2,1) ;
            else      
                targetstate = F*targetstate ; 
                targetstate([1,3]) = targetstate([1,3]) + location_std * randn(2,1) ; 
            end
            truth.X{k} = [truth.X{k} targetstate];
            truth.track_list{k} = [truth.track_list{k} ii];
            truth.N(k) = truth.N(k) + 1 ;
        end
        last_state(:,ii) = targetstate ; 
    else
        truth.spawnNextTimeTrack{tbirth(ii)-1} = [truth.spawnNextTimeTrack{tbirth(ii)-1} parent(1,ii)] ; 
        truth.newSpawn{tbirth(ii)} = [truth.newSpawn{tbirth(ii)} ii] ; 
        
        parent_state = last_state(:,parent(1,ii)) ; 
        delta = rad2deg(atan2(parent_state(4),parent_state(3)))+90 ; 
        dx = d_s * cosd(delta) ; 
        dy = d_s * sind(delta) ; 
        targetstate = parent_state ; 
        if parent(2,ii) == 1
            targetstate = targetstate + [dx; 0 ; dy; 0] ; 
        elseif parent(2,ii) == 2
            targetstate = targetstate - [dx; 0 ; dy; 0] ; 
        end
        k = tbirth(ii) ; 
        truth.X{k} = [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} ii];
        truth.N(k) = truth.N(k) + 1 ;
        for k = tbirth(ii)+1 : min(tdead(ii),truth.K)
            
            if rand>P_on
                targetstate([1,3]) = targetstate([1,3]) + location_std * randn(2,1) ; 
            else
                targetstate = F*targetstate ; 
                targetstate([1,3]) = targetstate([1,3]) + location_std * randn(2,1) ; 
            end
            truth.X{k} = [truth.X{k} targetstate];
            truth.track_list{k} = [truth.track_list{k} ii];
            truth.N(k) = truth.N(k) + 1 ;
        end
        last_state(:,ii) = targetstate ; 
    end
end

% plot the true plot
figure ; hold on ; 
colord = rand(truth.total_tracks,3) ; 
for k = 1 : truth.K
    for t = 1 : length(truth.track_list{k})
        id = truth.track_list{k}(t) ; 
        scatter(truth.X{k}(1,t),  truth.X{k}(3,t), 50, 'MarkerFaceColor', colord(id, :), ...
            'MarkerEdgeColor', colord(id, :)) ; 
    end
end
figure
stairs(1:truth.K,truth.N)
                           
