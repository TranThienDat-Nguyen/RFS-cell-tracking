%% Info
% Run the breast cancer cell tracking experiment in the paper
%% Prepare workspace
clear all
clc
addpath(genpath(pwd))
disp('Preparing workspace ...')
clear all ;
dbstop if error ;
warning('off','all') ;

%% Dataset setting
    % dataset settings
setting.dataset_name = 'MDA-MB-231' ; 
setting.dataset_img_size = [1 2560; 1 1920];
setting.kstart = 1 ; 
setting.kend = 88 ; 
%% (Optional) Only for synthetic cell tracking experiment

%% Dynamic setting
    % robust filtering options
setting.est_clutter = true; 
setting.est_Pd = true ; 
setting.detection_rate = 0.9 ; % use this value if setting.est_clutter = false
setting.lambda_c = 30 ;        % use this value if setting.est_Pd = false

    % dnamic settings: cardinality distribution given for each mode
% ----------------------------------------
%  \   Cardinality |  0  |  1  |  2  |
% Mode|  Survival  |  -     -     -  |   
%     |  Mitosis   |  -     -     -  |
% ----------------------------------------
setting.card_dist_mode = [0.03, 0.94, 0.03  ; 
                          0.05, 0.05, 0.90] ;
                      
    % dynamic settings: birth
setting.init_birth_rate = 0.03 ; 
setting.birth_rate = 0.01 ; 
setting.birth_loc_std = 10 ; 
setting.birth_vel_std = 10 ; 
setting.birth_init_st = [80 10] ; 
    % dynamic setting: survival transition
setting.survival_rate = 0.99;
setting.survival_w_CV_RDW = [0.3 0.7] ; % weights for CV vs RDW models
setting.survival_CV_vel_std = 5 ;       % noise standard deviation of the CV model
setting.survival_RDW_loc_std = 20 ;      % noise standard deviation for RDW model
    % dynamc setting: division
setting.spawn_init_st = [80, 10] ; 
setting.spawn_rate = 0.03 ;
setting.spawn_dist = 40 ;        % distance of spawns from parent
setting.spawn_loc_std = 10 ;      % noise on the location of the spawn
setting.spawn_vel_std = 10 ;      % noise on the velocity of the spawn
    % dynamic: clutter
setting.clutter_method = 'data_association' ; % 'data_association' or 'bootstrap' (bootstrap not implemented yet)
setting.clutter_PB = 0.5 ;                    % clutter detection rate
setting.clutter_PS = 0.5 ;                    % clutter survival rate
setting.clutter_PD = 0.9 ;                    % clutter detection rate
setting.clutter_gen = 30 ;                    % number of cluter generator
setting.clutter_init_st = [8 1] ;             % this setting only for 'bootstrap'
setting.clutter_Lc_birth = 1 ;                % this setting only for 'bootstrap'
setting.load_clutter = false ;                 % this setting only for 'data_association'
    % measurement model
setting.meas_std = 10 ;                        % noise standard deviation
%% Analyzing setting
setting.store_glmb = true ;                   % store the glmb density at each time step
setting.store_glmb_cap = 100 ;                % maximum number of components to be stored
setting.store_glmb_prune_thres = 0.01 ;       % only store glmb components greater than this threshold
setting.visualize = true ;                    % generate tracking movie
setting.compute_stats = true ;                % generate statistics of tracking results
setting.intensity_location_step = 20 ;        % sum the intensity over a set of pixels
setting.intensity_max_speed = 50 ;            % maximum speed to plot velocity heatmap (speed step is set to 1)
%% Evaluation setting (only when ground truth is available)
setting.noisy_eval = true;       % generate plots while evaluating
setting.ospa_c = 25 ;            % ospa cut-off
setting.ospa_order = 1 ;         % ospa norm-order
setting.ospa_win_length = 20 ;   % ospa(2) window length 
%% Generate model
disp('Generating model ...')
model = gen_model(setting) ; 
model.kalman_predict_multi_object = @kalman_predict_multi_object ; 
model.kalman_predict_mixture = @kalman_predict_mixture ; 
%% Generating clutter table
if strcmp(setting.clutter_method, 'data_association')
    if setting.load_clutter % provide your location of clutter tables
        disp('Loading clutter lookup table ...')
        init_clutter_table_location = fullfile('clut_tab_1000_1000', ...
            'init_clutter_table_z_1_100_cl_0_100.mat') ; 
        clutter_table_location = fullfile('clut_tab_1000_1000', ...
            'clutter_table_z_1_100_cl_0_100.mat') ; 
        load(init_clutter_table_location) ; 
        load(clutter_table_location) ; 
    else
        disp('Generating clutter lookup table ... (may take long!!!)')
        N_max_clt_births_first = 100 ; % initial birth clutter
        N_max_cl = 100 ; % maximum number of clutter
        N_top_cl = 100 ; % maximum number of clutter components (hypotheses)
        N_max_z = 100 ;  % maximum number of measurements
        N_min_z = 0 ;    % minimum number of measurements
        [init_clutter_table , clutter_table] = gen_clutter_lookup(model,  ...
                    N_max_clt_births_first, N_max_cl, N_top_cl, N_max_z, N_min_z) ; 
    end
end
%% Filter parameters
filter.H_upd = 50000 ;                 % requested number of updated components
filter.H_max = 50000 ;                 % cap on number of posterior components
filter.H_min = 0 ;                     % minimum number of components
filter.H_clutter = 100 ;                 % number of clutter components (for clutter data association method)
filter.hyp_threshold = 1e-5 ;          % pruning threshold for components
filter.P_G = 0.999999 ;                           %gate size in percentage
filter.gamma = chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag = 1;
filter.L_max= 50;                   % limit on number of Gaussians in each track
filter.elim_threshold= 1e-3;        % pruning threshold for Gaussians in each track
filter.merge_threshold= 4;          % merging threshold for Gaussians in each track


%% Main experiment
%% Initialize output folder
output_path = fullfile('output', setting.dataset_name) ; 
if ~exist(output_path,'dir')
    mkdir(output_path) ; 
end 
img_path = fullfile('datasets', 'MDA_MB_231') ;  % path to the images
%% Loading measurements
    load(fullfile('data', 'meas_exp3.mat')) ;  % load the meas data
%% Prediction Approximation
disp('Running prediction approximation filter ...')
kstart = 2 ; 
meas_ru = eps * ones(size(meas.Z{1},2),1) ;
tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
est.filter= filter;
est.X= cell(length(meas.Z),1);
est.N= zeros(length(meas.Z),1);
est.L= cell(length(meas.Z),1);
est.pD= cell(length(meas.Z),1);
est.nspawns= zeros(length(meas.Z),1);
est.nbirths= zeros(length(meas.Z),1);
est.nclutter= zeros(length(meas.Z),1);
est.nclutter_count= zeros(length(meas.Z),1);
est.cpu_time = zeros(length(meas.Z),1);
    % Initialise prior
glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
glmb_update.w= 1;               %vector of GLMB component/hypothesis weights
glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
glmb_update.I_old= {[]}; 
glmb_update.n= 0;               %vector of GLMB component/hypothesis cardinalities
glmb_update.cdn= 1;
glmb_update.total_births= 0;
glmb_update.total_spawns= 0;
glmb_update.N_clutter= 0;
glmb_update.N_clutter_count= 0;
% stuff for the estimate
glmb_update.est_X_pD = {zeros(model.x_dim+1,0)} ;
glmb_update.est_L = {[]} ;
glmb_update.est_others = {zeros(5,0)} ;
if setting.store_glmb
    store_glmb(setting, glmb_update, output_path, 1) ; 
end
for k = kstart:meas.K
    tic; 
    if k == kstart
        [glmb_update,meas_ru] = prediction_approximation(model, filter, glmb_update, tt_birth, 1 , meas , k,init_clutter_table) ; 
    else
        [glmb_update,meas_ru] = prediction_approximation(model, filter, glmb_update, tt_birth, 1 , meas , k,clutter_table) ;
    end
    tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru), 200) ;
    est_temp = extract_tracks(glmb_update, model) ; 
    est.X = est_temp.X ; 
    est.N = est_temp.N ; 
    est.L = est_temp.L ; 
    est.pD = est_temp.pD ; 
    est.nspawns = est_temp.nspawns ; 
    est.nbirths = est_temp.nbirths ; 
    est.nclutter = est_temp.nclutter ; 
    est.nclutter_count = est_temp.nclutter_count ; 
    H_posterior = length(glmb_update.w) ; 
    disp([ ' time= ' , num2str(k), ...
           ' #est_card = ' , num2str(est_temp.N(k)), ...
           ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
           ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
           ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
           ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
           ' #Hypo = ' , num2str(length(glmb_update.w))]) ;                                     
    est.cpu_time(k) = toc ; 
    est.H_pos(k) = H_posterior ; 
    % store the glmb 
    if setting.store_glmb
        store_glmb(setting, glmb_update, output_path, k) ; 
    end
end
% Analysing results
identifier = 'MDA-MB-231' ; 
plot_lineage(est, identifier) ; 
movie_out = gen_cell_movie(img_path, 'jpg', model, est) ; implay(movie_out) ;
gen_stats(setting, output_path, identifier) ; 
