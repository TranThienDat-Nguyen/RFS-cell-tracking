%% Info
% Run the synthetic cell tracking experiment in the paper
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
setting.dataset_name = 'synthetic' ; 
setting.dataset_img_size = [1 1000; 1 1000];
setting.kstart = 1 ; 
setting.kend = 100 ; 
%% (Optional) Only for synthetic cell tracking experiment
setting.ccd_noise = [0.05, 0.1, 0.2, 0.3, 0.4] ; 
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
                          0.01, 0.09, 0.90] ;
                      
    % dynamic settings: birth
setting.init_birth_rate = 0.03 ; 
setting.birth_rate = 0.01 ; 
setting.birth_loc_std = 10 ; 
setting.birth_vel_std = 10 ; 
setting.birth_init_st = [80 10] ; 
    % dynamic setting: survival transition
setting.survival_rate = 0.99;
setting.survival_w_CV_RDW = [0.3 0.7] ; % weights for CV vs RDW models
setting.survival_CV_vel_std = 1 ;       % noise standard deviation of the CV model
setting.survival_RDW_loc_std = 3 ;      % noise standard deviation for RDW model
    % dynamc setting: division
setting.spawn_init_st = [80, 10] ; 
setting.spawn_rate = 0.03 ;
setting.spawn_dist = 10 ;        % distance of spawns from parent
setting.spawn_loc_std = 3 ;      % noise on the location of the spawn
setting.spawn_vel_std = 3 ;      % noise on the velocity of the spawn
    % dynamic: clutter
setting.clutter_method = 'data_association' ; % 'data_association' or 'bootstrap' (bootstrap not implemented yet)
setting.clutter_PB = 0.5 ;                    % clutter detection rate
setting.clutter_PS = 0.9 ;                    % clutter survival rate
setting.clutter_PD = 0.9 ;                    % clutter detection rate
setting.clutter_gen = 30 ;                    % number of cluter generator
setting.clutter_init_st = [8 1] ;             % this setting only for 'bootstrap'
setting.clutter_Lc_birth = 1 ;                % this setting only for 'bootstrap'
setting.load_clutter = true ;                 % this setting only for 'data_association'
    % measurement model
setting.meas_std = 2 ;                        % noise standard deviation
%% Analyzing setting
setting.store_glmb = true ;                   % store the glmb density at each time step
setting.store_glmb_cap = 100 ;                % maximum number of components to be stored
setting.store_glmb_prune_thres = 0.01 ;       % only store glmb components greater than this threshold
setting.visualize = true ;                    % generate tracking movie
setting.compute_stats = true ;                % generate statistics of tracking results
setting.intensity_location_step = 10 ;        % sum the intensity over a set of pixels
setting.intensity_max_speed = 20 ;            % maximum speed to plot velocity heatmap (speed step is set to 1)
%% Evaluation setting (only when ground truth is available)
setting.noisy_eval = true;       % generate plots while evaluating
setting.ospa_c = 25 ;           % ospa cut-off
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
            'init_clutter_table_z_1_500_cl_0_500.mat') ; 
        clutter_table_location = fullfile('clut_tab_1000_1000', ...
            'clutter_table_z_1_500_cl_0_500.mat') ; 
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
filter.L_max= 50;                   % limit on number of Gaussians in each track
filter.elim_threshold= 1e-3;        % pruning threshold for Gaussians in each track
filter.merge_threshold= 4;          % merging threshold for Gaussians in each track

%% Load the truth and generate synthetic measurements
disp('Loading truth ...')
load('data/truth_exp2.mat') % load the truth data

%% Main experiment

for s = 1:5
    %% Initialize output folder
    output_path = fullfile('output', [setting.dataset_name, '_', num2str(s)]) ; 
    if ~isfolder(output_path)
        mkdir(output_path) ; 
    end 
    img_path = fullfile(output_path, 'img') ;  % path of the synthesized images
    disp('Generating synthesized measurements ...')
    if ~isfolder(img_path)
        mkdir(img_path) ; 
    end
    gen_synthesized_image(model, setting.ccd_noise(s), img_path, truth) ; 
    meas = gen_meas_synthetic(model, img_path) ;
%% Loading measurements
%     load(fullfile('data', 'meas_exp2', ['meas_', sprintf('%d', s), '.mat'])) ;  % load the meas data
    %% Prediction Approximation
    disp('Running prediction approximation filter ...')
    kstart = 2 ; 
    meas_ru = eps * ones(size(meas.Z{1},2),1) ;
    tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
    est_pa(s).filter= filter;
    est_pa(s).X= cell(length(meas.Z),1);
    est_pa(s).N= zeros(length(meas.Z),1);
    est_pa(s).L= cell(length(meas.Z),1);
    est_pa(s).pD= cell(length(meas.Z),1);
    est_pa(s).nspawns= zeros(length(meas.Z),1);
    est_pa(s).nbirths= zeros(length(meas.Z),1);
    est_pa(s).nclutter= zeros(length(meas.Z),1);
    est_pa(s).nclutter_count= zeros(length(meas.Z),1);
    est_pa(s).cpu_time = zeros(length(meas.Z),1);
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
        tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru)) ;
        est_temp = extract_tracks(glmb_update, model) ; 
        est_pa(s).X = est_temp.X ; 
        est_pa(s).N = est_temp.N ; 
        est_pa(s).L = est_temp.L ; 
        est_pa(s).pD = est_temp.pD ; 
        est_pa(s).nspawns = est_temp.nspawns ; 
        est_pa(s).nbirths = est_temp.nbirths ; 
        est_pa(s).nclutter = est_temp.nclutter ; 
        est_pa(s).nclutter_count = est_temp.nclutter_count ; 
        H_posterior = length(glmb_update.w) ; 
        disp([ ' time= ' , num2str(k), ...
               ' #est_card = ' , num2str(est_temp.N(k)), ...
               ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
               ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
               ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
               ' #truth_N = ' , num2str(truth.N(k)) , ...
               ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
               ' #Hypo = ' , num2str(length(glmb_update.w))]) ;                                     
        est_pa(s).cpu_time(k) = toc ; 
        est_pa(s).H_pos(k) = H_posterior ; 
        % store the glmb 
        if setting.store_glmb
            store_glmb(setting, glmb_update, output_path, k) ; 
        end
    end
    % Analysing results
    identifier = ['PA_', num2str(s)] ; 
    plot_cardinality(truth, est_pa(s), identifier) ;
    plot_lineage(est_pa(s), identifier) ; 
    eval_results(model, truth, meas, est_pa(s), setting.ospa_c, setting.ospa_order, setting.ospa_win_length, identifier, output_path) ; 
    movie_out = gen_cell_movie(img_path, 'tif', model, est_pa(s)) ; implay(movie_out) ;
    gen_stats(setting, output_path, identifier) ; 
    
    %% Approximation update strategy
    disp('Running update approximation filter ...')
    kstart = 2 ;
    meas_ru = eps * ones(size(meas.Z{1},2),1) ;
    tt_birth = gen_meas_driven_birth(model, meas, meas_ru, 1, 1:length(meas_ru));
    est_ua(s).filter= filter;
    est_ua(s).X= cell(length(meas.Z),1);
    est_ua(s).N= zeros(length(meas.Z),1);
    est_ua(s).L= cell(length(meas.Z),1);
    est_ua(s).pD= cell(length(meas.Z),1);
    est_ua(s).nspawns= zeros(length(meas.Z),1);
    est_ua(s).nbirths= zeros(length(meas.Z),1);
    est_ua(s).nclutter= zeros(length(meas.Z),1);
    est_ua(s).nclutter_count= zeros(length(meas.Z),1);
    est_ua(s).cpu_time = zeros(length(meas.Z),1);
    % Initialize prior
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
    for k = kstart:meas.K
        tic; 
        if k == kstart
            [glmb_update,meas_ru] = update_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,init_clutter_table) ; 
        else
            [glmb_update,meas_ru] = update_approximation(model , filter , glmb_update , tt_birth , 1 , meas , k,clutter_table) ;
        end
        tt_birth = gen_meas_driven_birth(model , meas, meas_ru, k , 1:length(meas_ru)) ;
        est_temp = extract_tracks(glmb_update, model) ; 
        est_ua(s).X = est_temp.X ; 
        est_ua(s).N = est_temp.N ; 
        est_ua(s).L = est_temp.L ; 
        est_ua(s).pD = est_temp.pD ; 
        est_ua(s).nspawns = est_temp.nspawns ; 
        est_ua(s).nbirths = est_temp.nbirths ; 
        est_ua(s).nclutter = est_temp.nclutter ; 
        est_ua(s).nclutter_count = est_temp.nclutter_count ; 
        H_posterior = length(glmb_update.w) ; 
%         display_diaginfo(k,est_ju(s),filter,H_posterior);
        disp([ ' time= ' , num2str(k), ...
               ' #est_card = ' , num2str(est_temp.N(k)), ...
               ' #est_spawns = ' , num2str(est_temp.nspawns(k)), ...
               ' #est_births = ' , num2str(est_temp.nbirths(k)), ...
               ' #est_clutter = ' , num2str(est_temp.nclutter_count(k)), ...
               ' #truth_N = ' , num2str(truth.N(k)) , ...
               ' #meas = ' , num2str(size(meas.Z{k},2)) , ...
               ' #Hypo = ' , num2str(length(glmb_update.w))]) ;                                     
        est_ua(s).cpu_time(k) = toc ; 
        est_ua(s).H_pos(k) = H_posterior ; 
        % store the glmb 
        if setting.store_glmb
            store_glmb(setting, glmb_update, output_path, k) ; 
        end
    end
    % Analysing results
    identifier = ['UA_', num2str(s)] ; 
    plot_cardinality(truth, est_ua(s), identifier) ;
    plot_lineage(est_ua(s), identifier) ; 
    eval_results(model, truth, meas, est_ua(s), setting.ospa_c, setting.ospa_order, setting.ospa_win_length, identifier, output_path) ; 
    movie_out = gen_cell_movie(img_path, 'tif', model, est_ua(s)) ; implay(movie_out) ;
    gen_stats(setting, output_path, identifier) ; 
end