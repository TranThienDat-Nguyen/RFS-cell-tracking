function model = gen_model(setting)
%% General setting
    model.unknown_clutter = setting.est_clutter; 
    model.unknown_pD = setting.est_Pd ;
    model.pD = setting.detection_rate ; 
    model.lambda_c = setting.lambda_c ; 
    model.card_dist_mode = setting.card_dist_mode ; 
%% Set up surveillance region
    model.range_c = setting.dataset_img_size; % region
%% General dynamic model
    model.x_dim = 4; % state vector dimension
    model.z_dim = 2; % observation vector dimension
    model.T = 1 ; % transition time
    model.xy_pos = [1,3] ; % position of x,y in the state vector
    model.x_loc = 1 ; 
    model.y_loc = 3 ; 
%% Set up birth dynamic
    model.pB_max = setting.birth_rate ;
    model.pB_init = setting.init_birth_rate ;
    model.pB = setting.birth_rate ;
    model.qB = 1 - model.pB ;
    model.lambda_b = 1 ; % expected number of birth per scan
    model.nwb = 1;
    model.wb = 1;
    temp = diag([setting.birth_loc_std , setting.birth_vel_std , setting.birth_loc_std , setting.birth_vel_std]) ; 
    model.Q_b = temp * temp';
%% Setup survival dynamic
    model.ws = setting.survival_w_CV_RDW ; 
    model.pS = setting.survival_rate;
    model.qS = 1 - model.pS ;
    % CV model
    model.A0= [ 1 model.T; 0 1 ];                         %transition matrix                     
    model.F(:,:,1)=  [ model.A0 zeros(2,2); zeros(2,2) model.A0 ];
    model.B0= [ (model.T^2)/2; model.T ];
    model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
    model.Q(:,:,1)= (setting.survival_CV_vel_std)^2 * model.B * model.B';   %process noise covariance
    % RDW model
    model.F(:,:,2) = diag([1 1 1 1]) ; 
    model.Q(:,:,2) = diag([setting.survival_RDW_loc_std^2, 0, setting.survival_RDW_loc_std^2, 0]) ; 
%% Setup spawning dynamic
    model.pSp = setting.spawn_rate ; % damping the spawn rate down with this factor 
    model.nws = 9; % number of Gaussian components in the spawn density
    temp = [1   0   1   0] ; 
    model.S = diag(temp) ;
    model.x_shift = zeros(1,model.nws) ; 
    model.y_shift = zeros(1,model.nws) ; 
    delta_angle = [20 40 60 80 100 120 140 160 180] ; % angles to generate spawns GM in degree (upper), for lower + 180
    model.wtsp = 1 ; % 1 transition matrix
    model.wsp = (1/model.nws)*ones(1,model.nws) ; 
    model.spawn_dist = setting.spawn_dist; 
    for cidx = 1 : model.nws
        model.x_shift(cidx) = setting.spawn_dist * cosd(delta_angle(cidx)) ; 
        model.y_shift(cidx) = setting.spawn_dist * sind(delta_angle(cidx)) ; 
    end
    temp = [setting.spawn_loc_std^2 , setting.spawn_vel_std^2 , setting.spawn_loc_std^2 , setting.spawn_vel_std^2 ] ;
    model.Q_s = diag(temp) ; 
    model.spawn_init_st = setting.spawn_init_st ; 
%% Setup measurement model
    model.init_st = setting.birth_init_st ; 
    if ~setting.est_Pd 
        model.pD = setting.detection_rate; % detection rate
        model.qD = 1 - model.pD ; % Misdetection rate
        model.unknown_pD = false ; 
    end
    % Observation model
    model.H = [ 1 0 0 0   ; 0 0 1 0];    %observation matrix
    model.sigma_d = setting.meas_std;
    model.D= diag([ model.sigma_d; model.sigma_d ]); 
    model.R= model.D*model.D';              %observation noise covariance
%% Setup clutter model
    model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1));
    if ~setting.est_clutter
        model.lambda_c = setting.lambda_c  ;
        model.unknown_clutter = false ;
    end
    % model for CPHD - robust bootstrap filtering
    if strcmp(setting.clutter_method, 'bootstrap')
        model.cphd.P_S = setting.survival_rate ; 
        model.cphd.Q_S= 1-setting.survival_rate;
        model.cphd.lambda_cb = 1;       % birth rate for clutter targets
        model.cphd.w_cb = 1;
        model.cphd.clutter_P_S = setting.clutter_PS;    % survival probability for clutter targets
        model.cphd.clutter_P_D = setting.clutter_PD;     % detection probability for clutter targets
        model.cphd.Lc_birth = setting.clutter_Lc_birth; 
        model.cphd.u_cb = setting.clutter_init_st(1); 
        model.cphd.v_cb = setting.clutter_init_st(2);
        model.cphd.clutter_Nt = setting.clutter_gen;     % number of clutter generators
    elseif strcmp(setting.clutter_method, 'data_association')
        model.clutter_P_S = setting.clutter_PS ;
        model.clutter_P_D = setting.clutter_PD ;
        model.clutter_P_B = setting.clutter_PB ;
        model.clutter_gen = setting.clutter_gen ;
    end
end




            
     
     