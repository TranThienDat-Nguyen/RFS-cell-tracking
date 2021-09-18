%% Info
% Run the simulation experiment in the paper
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
setting.dataset_name = 'simulation' ; 
setting.dataset_img_size = [1 1000; 1 1000];
setting.kstart = 1 ; 
setting.kend = 100 ; 
%% (Optional) Only for simulation setting
setting.sim_pD = 0.9 ; 
setting.sim_lambda_c = 30 ;
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
setting.card_dist_mode = [0.01, 0.98, 0.01  ; 
                          0.01, 0.10, 0.89] ;
                      
    % dynamic settings: birth
setting.init_birth_rate = 0.03 ; 
setting.birth_rate = 0.01 ; 
setting.birth_loc_std = 10 ; 
setting.birth_vel_std = 10 ; 
setting.birth_init_st = [80 10] ; 
    % dynamic setting: survival transition
setting.survival_rate = 0.99;
setting.survival_w_CV_RDW = [0.7 0.3] ; % weights for CV vs RDW models
setting.survival_CV_vel_std = 1 ;       % noise standard deviation of the CV model
setting.survival_RDW_loc_std = 3 ;      % noise standard deviation for RDW model
setting.survival_init_st = [80 10] ;    % initial beta parameters for objects
    % dynamc setting: division
setting.spawn_init_st = [80 , 10] ; 
setting.spawn_rate = 0.03 ;
setting.spawn_dist = 10 ;        % distance of spawns from parent
setting.spawn_loc_std = 3 ;      % noise on the location of the spawn
setting.spawn_vel_std = 3 ;      % noise on the velocity of the spawn
    % dynamic: clutter
setting.clutter_method = 'data_association' ; % 'data_association' or 'bootstrap' (bootstrap not implemented yet)
setting.clutter_PB = 0.5 ;                    % clutter detection rate
setting.clutter_PS = 0.9 ;                    % clutter survival
setting.clutter_PD = 0.9 ;                    % clutter detection rate
setting.clutter_gen = 30 ;                    % number of cluter generator
setting.clutter_init_st = [8 1] ;           % initial beta parameters for unknown pD clutter
setting.clutter_Lc_birth = 1 ;                % this setting only for 'bootstrap'
setting.load_clutter = true ;                 % this setting only for 'data_association'
    % measurement model
setting.meas_std = 5 ;                        % noise standard deviation
%% Analyzing setting
setting.store_glmb = false ; 
setting.store_glmb_cap = 100 ; 
setting.store_glmb_prune_thres = 0.01 ; 
setting.visualize = false ; % generate tracking movie
setting.compute_stats = false ;
setting.intensity_location_step = 10 ;
setting.intensity_max_speed = 20 ; 
%% Evaluation setting (only when ground truth is available)
setting.ospa_c = 100 ;           % ospa cut-off
setting.ospa_order = 1 ;         % ospa norm-order
setting.ospa_win_length = 20 ;   % ospa(2) window length
setting.noisy_eval = true ;      % generate plots when evaluating
%% Initialize output folder
output_path = fullfile('output', setting.dataset_name) ; 
if ~exist(output_path,'dir')
    mkdir(output_path) ; 
end  
%% Generate model
disp('Generating model ...')
model = gen_model(setting) ; 
% choose spawning with parent bearing model for this experiment
model.kalman_predict_multi_object = @kalman_predict_multi_object_parent_bearing ; 
model.kalman_predict_mixture = @kalman_predict_mixture_parent_bearing ; 
model.nws = 1 ; model.wsp = 1 ; 
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
        [init_clutter_table , clutter_table] = gen_clutter_lookup_clean(model,  ...
                    N_max_clt_births_first, N_max_cl, N_top_cl, N_max_z, N_min_z) ; 
    end
end
%% Filter parameters
filter.H_upd = 30000 ;                 % requested number of updated components
filter.H_max = 30000 ;                 % cap on number of posterior components
filter.H_min = 0 ;                     % minimum number of components
filter.H_clutter = 2 ;                 % number of clutter components (for clutter data association method)
filter.hyp_threshold = 1e-5 ;          % pruning threshold for components
filter.P_G = 0.999999 ;                           %gate size in percentage
filter.gamma = chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag = 1;
filter.run_flag = 'noisy' ;         % display tracking info
filter.L_max= 50;                   % limit on number of Gaussians in each track
filter.elim_threshold= 1e-3;        % pruning threshold for Gaussians in each track
filter.merge_threshold= 4;          % merging threshold for Gaussians in each track

%% Load truth
disp('Load truth ...')
load('data/truth_exp1.mat') % load the truth data

%% Main filtering
disp('Running filter ...')
kstart = 2 ;
sim_pD = setting.sim_pD ; 
sim_lambda_c = setting.sim_lambda_c ; 
model.pD = sim_pD ; 
model.lambda_c = sim_lambda_c ; 
for r = 1:1
    %% Generating measurements
    meas = gen_meas_sim(truth, sim_pD, sim_lambda_c, model.sigma_d, model.range_c) ;  
	%% Initialisation
    %output variables
    est_pa(r).X= cell(length(meas.Z),1);
    est_pa(r).N= zeros(length(meas.Z),1);
    est_pa(r).L= cell(length(meas.Z),1);
    est_pa(r).pD= cell(length(meas.Z),1);
    est_pa(r).nspawns= zeros(length(meas.Z),1);
    est_pa(r).nbirths= zeros(length(meas.Z),1);
    est_pa(r).nclutter= zeros(length(meas.Z),1);
    est_pa(r).nclutter_count= zeros(length(meas.Z),1);
    est_pa(r).cpu_time = zeros(length(meas.Z),1);
    est_pa(r).H_pos = zeros(length(meas.Z),1);
    est_pa(r).filter= filter;
   
    %output variables
    est_ua(r).X= cell(length(meas.Z),1);
    est_ua(r).N= zeros(length(meas.Z),1);
    est_ua(r).L= cell(length(meas.Z),1);
    est_ua(r).pD= cell(length(meas.Z),1);
    est_ua(r).nspawns= zeros(length(meas.Z),1);
    est_ua(r).nbirths= zeros(length(meas.Z),1);
    est_ua(r).nclutter= zeros(length(meas.Z),1);
    est_ua(r).nclutter_count= zeros(length(meas.Z),1);
    est_ua(r).cpu_time = zeros(length(meas.Z),1);
    est_ua(r).H_pos = zeros(length(meas.Z),1);
    est_ua(r).filter= filter;
    
    %output variables
    est_ef(r).X= cell(length(meas.Z),1);
    est_ef(r).N= zeros(length(meas.Z),1);
    est_ef(r).L= cell(length(meas.Z),1);
    est_ef(r).pD= cell(length(meas.Z),1);
    est_ef(r).nspawns= zeros(length(meas.Z),1);
    est_ef(r).nbirths= zeros(length(meas.Z),1);
    est_ef(r).nclutter= zeros(length(meas.Z),1);
    est_ef(r).nclutter_count= zeros(length(meas.Z),1);
    est_ef(r).cpu_time = zeros(length(meas.Z),1);
    est_ef(r).H_pos = zeros(length(meas.Z),1);
    est_ef(r).filter= filter;
    
    %% Run main filtering
        %% Approximation prediction strategy
    disp('Running prediction approximation filter ...')
    meas_ru = eps * ones(size(meas.Z{1},2),1) ;
    tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
    glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
    glmb_update.w= 1;               %vector of GLMB component/hypothesis weights
    glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
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
    for k = kstart:meas.K
        tic; 
        if k == kstart
            [glmb_update,meas_ru] = prediction_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,init_clutter_table) ; 
        else
            [glmb_update,meas_ru] = prediction_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,clutter_table) ;
        end
        tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru)) ;
        est_temp = extract_tracks(glmb_update, model) ; 
        est_pa(r).X = est_temp.X ; 
        est_pa(r).N = est_temp.N ; 
        est_pa(r).L = est_temp.L ; 
        est_pa(r).pD = est_temp.pD ; 
        est_pa(r).nspawns = est_temp.nspawns ; 
        est_pa(r).nbirths = est_temp.nbirths ; 
        est_pa(r).nclutter = est_temp.nclutter ; 
        est_pa(r).nclutter_count = est_temp.nclutter_count ; 
        H_posterior = length(glmb_update.w) ; 
%         display_diaginfo(k,est(r),filter,H_posterior);
        disp([ ' time= ' , num2str(k), ...
               ' #est_card = ' , num2str(est_temp.N(k)), ...
               ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
               ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
               ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
               ' #truth_N = ' , num2str(truth.N(k)) , ...
               ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
               ' #Hypo = ' , num2str(length(glmb_update.w))]) ;                                     
        est_pa(r).cpu_time(k) = toc ; 
        est_pa(r).H_pos(k) = H_posterior ; 
        % store the glmb 
        if setting.store_glmb
            store_glmb(glmb_update, fullfile(output_path, 'glmb_exact'), k) ; 
        end
    end
    
        %% Approximation update strategy
    disp('Running update approximation filter ...')
    meas_ru = eps * ones(size(meas.Z{1},2),1) ;
    tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
    glmb_update.tt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
    glmb_update.w= 1;               %vector of GLMB component/hypothesis weights
    glmb_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
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
    for k = kstart:meas.K
        tic; 
        if k == kstart
            [glmb_update,meas_ru] = update_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,init_clutter_table) ; 
        else
            [glmb_update,meas_ru] = update_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,clutter_table) ;
        end
        tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru)) ;
        est_temp = extract_tracks(glmb_update, model) ; 
        est_ua(r).X = est_temp.X ; 
        est_ua(r).N = est_temp.N ; 
        est_ua(r).L = est_temp.L ; 
        est_ua(r).pD = est_temp.pD ; 
        est_ua(r).nspawns = est_temp.nspawns ; 
        est_ua(r).nbirths = est_temp.nbirths ; 
        est_ua(r).nclutter = est_temp.nclutter ; 
        est_ua(r).nclutter_count = est_temp.nclutter_count ; 
        H_posterior = length(glmb_update.w) ; 
%         display_diaginfo(k,est_ju(r),filter,H_posterior);
        disp([ ' time= ' , num2str(k), ...
               ' #est_card = ' , num2str(est_temp.N(k)), ...
               ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
               ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
               ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
               ' #truth_N = ' , num2str(truth.N(k)) , ...
               ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
               ' #Hypo = ' , num2str(length(glmb_update.w))]) ;                                     
        est_ua(r).cpu_time(k) = toc ; 
        est_ua(r).H_pos(k) = H_posterior ; 
    end
    
        %% Exact filtering strategy
    disp('Running exact filter ...')
    meas_ru = eps * ones(size(meas.Z{1},2),1) ;
    tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
    mixture_update.jtt= cell(0,1);      %track table for GLMB (cell array of structs for individual tracks)
    mixture_update.w= 1;               %vector of GLMB component/hypothesis weights
    mixture_update.I= {[]};            %cell of GLMB component/hypothesis labels (labels are indices/entries in track table)
    mixture_update.n= 0;               %vector of GLMB component/hypothesis cardinalities
    mixture_update.cdn= 1;
    mixture_update.total_births= 0;
    mixture_update.total_spawns= 0;
    mixture_update.N_clutter= 0;
    mixture_update.N_clutter_count= 0;
        % stuff for the estimate
    mixture_update.est_X_pD = {zeros(model.x_dim+1,0)} ;
    mixture_update.est_L = {[]} ;
    mixture_update.est_others = {zeros(5,0)} ;
    for k = kstart:meas.K
        tic; 
        if k == kstart
            [mixture_update,meas_ru] = exact_filter(model , filter , mixture_update , tt_birth , 1 , meas , k,init_clutter_table) ; 
        else
            [mixture_update,meas_ru] = exact_filter(model , filter , mixture_update , tt_birth , 1 , meas , k,clutter_table) ;
        end
        tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru)) ;
        est_temp = extract_tracks(mixture_update, model) ; 
        est_ef(r).X = est_temp.X ; 
        est_ef(r).N = est_temp.N ; 
        est_ef(r).L = est_temp.L ; 
        est_ef(r).pD = est_temp.pD ; 
        est_ef(r).nspawns = est_temp.nspawns ; 
        est_ef(r).nbirths = est_temp.nbirths ; 
        est_ef(r).nclutter = est_temp.nclutter ; 
        est_ef(r).nclutter_count = est_temp.nclutter_count ; 
        H_posterior = length(mixture_update.w) ; 
%         display_diaginfo(k,est_ju(r),filter,H_posterior);
        disp([ ' time= ' , num2str(k), ...
               ' #est_card = ' , num2str(est_temp.N(k)), ...
               ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
               ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
               ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
               ' #truth_N = ' , num2str(truth.N(k)) , ...
               ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
               ' #Hypo = ' , num2str(length(mixture_update.w))]) ;                                     
        est_ef(r).cpu_time(k) = toc ; 
        est_ef(r).H_pos(k) = H_posterior ; 
    end
    
    
end
%% Plotting results
disp('Plotting results ...')
identifier = 'PA' ; 
plot_estimate(model, truth, meas, est_pa, identifier) ; 
eval_results(model, truth, meas, est_pa, setting.ospa_c, setting.ospa_order, setting.ospa_win_length, identifier, output_path) ; 

identifier = 'UA' ; 
plot_estimate(model, truth, meas, est_ua, identifier) ; 
eval_results(model, truth, meas, est_ua, setting.ospa_c, setting.ospa_order, setting.ospa_win_length, identifier, output_path) ; 

identifier = 'EF' ; 
plot_estimate(model, truth, meas, est_ef, identifier) ; 
eval_results(model, truth, meas, est_ef, setting.ospa_c, setting.ospa_order, setting.ospa_win_length, identifier, output_path) ; 
%% END


