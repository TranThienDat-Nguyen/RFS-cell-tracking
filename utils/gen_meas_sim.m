% generate measurements
function meas = gen_meas_fixed(truth, P_D, lambda_c, meas_noise, range_c)
meas.K = truth.K ; 
meas.Z = cell(truth.K,1) ; 
meas.Z_normal = cell(truth.K,1) ; 
meas.Z_spawn = cell(truth.K,1) ; 
meas.Z_clutter = cell(truth.K,1) ; 
for k = 1 : truth.K
    detectedidx= find( rand(truth.N(k),1) <= P_D );
    numDetected = length(detectedidx) ; 
    measNoise = diag(meas_noise*ones(1,2))  * randn(2 ,  numDetected) ; 
    meas.Z{k} = truth.X{k}([1,3],detectedidx)+measNoise ; 
    temp_tracklist = truth.track_list{k}(detectedidx) ;
    % Generate the rate
        spawntracks = truth.spawnNextTimeTrack{k} ; 
        if ~isempty(spawntracks)
            [~,spawntracksIdx] = ismember(spawntracks , temp_tracklist) ;
            spawntracksIdx(spawntracksIdx==0) = [] ; 
        else
            spawntracksIdx = [] ; 
        end
        normaltracksIdx = setdiff(1:numDetected , spawntracksIdx) ; 
        N_c = poissrnd(lambda_c) ;
%         N_c = clutterRate(k) ; 
        C = repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(2,N_c);
        meas.Z{k} = [meas.Z{k} , C] ; 
        numMeas = size(meas.Z{k},2) ; 
        [a , b] = grab_beta_param(0.9 , 0.1) ; 
        score_high = betarnd(a , b , [numMeas ,1]) ; 
%         score_high = 0.8 * ones(numMeas,1) ; 
        [a , b] = grab_beta_param(0.2 , 0.1) ; 
        score_low_1 = betarnd(a , b , [numMeas ,1]) ; 
%         score_low_1 = 0.1 * ones(numMeas,1) ; 
        [a , b] = grab_beta_param(0.4 , 0.1) ; 
        score_low_2 = betarnd(a , b , [numMeas ,1]) ; 
        
        [a , b] = grab_beta_param(0.1 , 0.1) ; 
        score_low_3 = betarnd(a , b , [numMeas ,1]) ; 
%         score_low_2 = 0.2 * ones(numMeas,1) ; 
        temp_Z_normal = zeros(1 , numMeas) ; 
        temp_Z_spawn = zeros(1 , numMeas) ; 
        temp_Z_clutter = zeros(1 , numMeas) ; 
        % set up rates for normal cells
        temp_Z_normal(normaltracksIdx) = score_high(normaltracksIdx) ;
        temp_Z_spawn(normaltracksIdx) = score_low_1(normaltracksIdx) ; 
        temp_Z_clutter(normaltracksIdx) = score_low_2(normaltracksIdx) ;
        % set up rates for spawning cells
        temp_Z_normal(spawntracksIdx) = score_low_1(spawntracksIdx) ;
        temp_Z_spawn(spawntracksIdx) = score_high(spawntracksIdx) ; 
        temp_Z_clutter(spawntracksIdx) = score_low_2(spawntracksIdx) ;
        % set up rates for clutter
        temp_Z_normal(numDetected+1 : end) = score_low_2(numDetected+1 : end) ;
        temp_Z_spawn(numDetected+1 : end) = score_low_3(numDetected+1 : end) ; 
        temp_Z_clutter(numDetected+1 : end) = score_high(numDetected+1 : end) ;
        meas.Z_normal{k} = zeros(1 , numMeas) ; 
        meas.Z_spawn{k} = zeros(1 , numMeas) ; 
        meas.Z_clutter{k} = zeros(1 , numMeas) ; 
        for midx = 1 : numMeas
%             temp_sum = temp_Z_normal(midx) + temp_Z_spawn(midx) ; 
            meas.Z_normal{k}(midx) = temp_Z_normal(midx) ; 
            meas.Z_spawn{k}(midx) = temp_Z_spawn(midx) ; 
            meas.Z_clutter{k}(midx) = temp_Z_clutter(midx) ; 
        end
end
end
    
function [a , b] = grab_beta_param(mu , sigma)
    a = -mu * (sigma^2 + mu^2 - mu)/sigma^2 ; 
    b = (sigma^2 + mu^2 - mu) * (mu - 1) / sigma^2 ; 
end
    
