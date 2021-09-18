function meas = gen_meas_synthetic(model, img_path)
    dir_info = dir(fullfile(img_path, '*.tif')) ;
    K = length(dir_info) ; 
    img_sequence = zeros(model.range_c(1,2), model.range_c(2,2), K) ; 
    for ii = 1 : K
        img_sequence(:,:,ii) = imread(fullfile(img_path,dir_info(ii).name)) ; 
    end
    MPHD_P = [3 ,3] ; 
    AV_rate = 1; 
    [Temp_Spots,Spots,Temp_Background,~]=MPMD_Enhancement_Sequences(img_sequence,MPHD_P,AV_rate);
    centroids = cell(1, K) ; 
    Iz = cell(1, K) ; 
    for f=1:K
        [centroids{1,f},Iz{1,f}]=Spot_Detector_MPHD(Spots(:,:,f),Temp_Spots(:,:,f)+Temp_Background(:,:,f),Temp_Background(:,:,f),10,'Max',10,0);
    end
    % then form the meas structure
    meas.K = K ;
    meas.Z = cell(K, 1) ; 
    meas.Z_normal = cell(K,1) ; 
    meas.Z_spawn = cell(K, 1) ; 
    meas.Z_clutter = cell(K,1) ; 
    meas.Z_I = cell(K,1) ; 
    for k = 1 : K
        meas.Z{k} = centroids{k}' ;
        N = length(Iz{k}) ; 
        meas.Z_normal{k} = zeros(1, N) ; 
        meas.Z_spawn{k} = zeros(1, N) ;
        meas.Z_clutter{k} = zeros(1, N) ; 
        meas.Z_I{k} = Iz{k} ; 
        for n = 1 : length(Iz{k})
            meas.Z_normal{k}(1,n) = max(min(gauss2mf(Iz{k}(n), [5, 10, 5, 25]), 0.9), 0.1) ;
            meas.Z_spawn{k}(1,n) = max(min(smf(Iz{k}(n), [10, 40]), 0.9), 0.1) ; 
            meas.Z_clutter{k}(1,n) = max(min(smf(-Iz{k}(n)+15, [5 10]), 0.9), 0.1) ;
        end
    end
end
