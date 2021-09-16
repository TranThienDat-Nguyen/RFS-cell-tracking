% run prediction apprx to generate hypotheses .. 
% same as before exept doesnt have "ah" field and dont apply Gaussian
% cap/prune/merge
function hypo = gen_hypo(model, tt_birth, r_b,  tt_in, meas, N_clutter_in, sampling_factor, filter, k)
    % this funtion performs the sampling step ..
    numBirths = length(tt_birth) ; 
    numprevTT = length(tt_in) ;
    tt_surv = cell(numprevTT , 1) ;
    tt_sp1 = cell(numprevTT , 1) ;
    tt_sp2 = cell(numprevTT , 1) ; 
    parfor tabidx = 1 : numprevTT
        % Survivals
            % predict the parameters for pD beta
        ss = tt_in{tabidx}.st(1) ; 
        tt = tt_in{tabidx}.st(2) ; 
        beta_mean = ss / (ss + tt);
        beta_var = (ss*tt) / ((ss+tt)^2*(ss+tt+1)); 
        tt_surv{tabidx}.st = zeros(1,2) ;  
        tt_surv{tabidx}.st(1) = ( ( beta_mean * (1-beta_mean)/ beta_var) - 1 ) * beta_mean;
        tt_surv{tabidx}.st(2) = ( ( beta_mean * (1-beta_mean)/ beta_var - 1) ) * (1 - beta_mean); 
            % predict the positions
        [wtemp_predict , mtemp_predict , Ptemp_predict] = model.kalman_predict_multi_object(model , tt_in{tabidx}.w , tt_in{tabidx}.m , tt_in{tabidx}.P , [] , true) ; 
        tt_surv{tabidx}.w= wtemp_predict/sum(wtemp_predict);
        tt_surv{tabidx}.m= mtemp_predict;                                                                                   
        tt_surv{tabidx}.P= Ptemp_predict;   
            % predict other states
        tt_surv{tabidx}.l= tt_in{tabidx}.l ; 
        tt_surv{tabidx}.mode = [(1-model.pSp) , model.pSp] ;
        
            % First spawns
        tt_sp1{tabidx}.st = model.init_st ;
        [wtemp_predict , mtemp_predict , Ptemp_predict] = model.kalman_predict_multi_object(model , tt_in{tabidx}.w , tt_in{tabidx}.m , tt_in{tabidx}.P , 1 , true) ; 
        tt_sp1{tabidx}.w = wtemp_predict/sum(wtemp_predict) ; 
        tt_sp1{tabidx}.m= mtemp_predict{1};
        tt_sp1{tabidx}.P= Ptemp_predict{1};        
        tt_sp1{tabidx}.l= [tt_in{tabidx}.l;k;1] ;  
        tt_sp1{tabidx}.mode = [(1-model.pSp) , model.pSp] ;
        
            % Second spawns
        tt_sp2{tabidx}.st = model.init_st;
        tt_sp2{tabidx}.w = wtemp_predict/sum(wtemp_predict) ;
        tt_sp2{tabidx}.m= mtemp_predict{2};
        tt_sp2{tabidx}.P= Ptemp_predict{2};
        tt_sp2{tabidx}.l= [tt_in{tabidx}.l;k;2] ; 
        tt_sp2{tabidx}.mode = [(1-model.pSp) , model.pSp] ;
    end
    
    % restructuring tt_predict
    tt_predict = cell(numprevTT*3 + numBirths , 1) ;
    tt_predict(1:numBirths) = tt_birth ; 
    for tabidx = 1 : numprevTT
        tabssidx = numBirths+ (tabidx-1)*3+1 ; 
        tt_predict(tabssidx) = tt_surv(tabidx) ; 
        tt_predict(tabssidx+1) = tt_sp1(tabidx) ; 
        tt_predict(tabssidx+2) = tt_sp2(tabidx) ;
    end
    
    %% Update track tables
    % gating the measurements by location
    for tabidx=1:length(tt_predict)
        tt_predict{tabidx}.gatemeas= gate_meas_gms_idx(meas.Z{k},filter.gamma,model,tt_predict{tabidx}.m,tt_predict{tabidx}.P);
    end
    
    % compute update track tables
    m= size(meas.Z{k},2); 
    L_TT_P = length(tt_predict) ; 
    tt_update= cell(L_TT_P*(m+1),1) ; 
        % missed detection tracks
    model_unknown_pD = model.unknown_pD ; 
    model_pD = model.pD ;
    model_ws = model.ws ;
    model_wsp = model.wsp ;
    parfor tabidx= 1:L_TT_P  % parallelizable
        tt_update{tabidx}= tt_predict{tabidx};
        if model_unknown_pD
            tt_update{tabidx}.pD = tt_predict{tabidx}.st(1)/(tt_predict{tabidx}.st(1) + tt_predict{tabidx}.st(2));
        else
            tt_update{tabidx}.pD = model_pD ; 
        end
        % update s,t
        tt_update{tabidx}.st(2) = tt_update{tabidx}.st(2)+1 ; 
        tt_update{tabidx}.qz_mode = 1 ; 
        % fabricate the weight for model prediction 
        if tabidx<= numBirths
            tt_update{tabidx}.qz = 1 ; 
        elseif mod(tabidx-numBirths,3) == 1
            tt_update{tabidx}.qz = model_ws ; 
        elseif mod(tabidx-numBirths,3) == 2 || mod(tabidx-numBirths,3) == 0 
            tt_update{tabidx}.qz = model_wsp ; 
        end
            
    end

        %measurement updated tracks (all pairs) 
    log_allcostm_cell = cell(1,L_TT_P) ; 
    tt_update_para = cell(1,L_TT_P) ; 
    parfor tabidx= 1:L_TT_P % parallelizable
        tt_update_para{tabidx} = cell(1,m) ; 
        if model_unknown_pD
            curr_pD = tt_predict{tabidx}.st(1)/(tt_predict{tabidx}.st(1) + tt_predict{tabidx}.st(2));
        else
            curr_pD = model_pD ; 
        end
        log_allcostm_cell{tabidx} = log(zeros(1,m)) ; 
        for mmidx = 1 : length(tt_predict{tabidx}.gatemeas)
            emm= tt_predict{tabidx}.gatemeas(mmidx) ; 
            % update the detection rate
            tt_update_para{tabidx}{emm}.pD = curr_pD ;  
            tt_update_para{tabidx}{emm}.st = tt_predict{tabidx}.st ; 
            tt_update_para{tabidx}{emm}.st(1) = tt_update_para{tabidx}{emm}.st(1) + 1 ; 
            % update kinematic state
            [w_temp,m_temp,P_temp] = kalman_update_single_object(meas.Z{k}(:,emm),model, ...
                tt_predict{tabidx}.w, tt_predict{tabidx}.m, tt_predict{tabidx}.P); 
            tt_update_para{tabidx}{emm}.qz = sum(w_temp) ; 
            log_allcostm_cell{tabidx}(emm) = log(sum(w_temp)) + log(curr_pD);
            w_temp = w_temp/sum(w_temp) ; 
            w_temp= w_temp+eps;  
            tt_update_para{tabidx}{emm}.w = w_temp ; 
            tt_update_para{tabidx}{emm}.m = m_temp ; 
            tt_update_para{tabidx}{emm}.P = P_temp ;
            % update other states
            tt_update_para{tabidx}{emm}.l = tt_predict{tabidx}.l ;
            g_normal = meas.Z_normal{k}(emm) ; 
            g_spawn = meas.Z_spawn{k}(emm) ;
            temp = tt_predict{tabidx}.mode .* [g_normal g_spawn]; 
            tt_update_para{tabidx}{emm}.qz_mode = sum(temp) ; 
            log_allcostm_cell{tabidx}(emm)  = log_allcostm_cell{tabidx}(emm)  + log(sum(temp)) ; 
            tt_update_para{tabidx}{emm}.mode = temp/sum(temp);
        end
    end
    % restructuring
    log_allcostm= log(zeros(L_TT_P,m)); 
    for tabidx =1  : L_TT_P
        for emm = tt_predict{tabidx}.gatemeas
            stoidx= L_TT_P*emm + tabidx;
            tt_update{stoidx} = tt_update_para{tabidx}{emm} ; 
            log_allcostm(tabidx,emm) = log_allcostm_cell{tabidx}(emm) ; 
        end
    end
    %% Compute the rates
        % detection rate
    avpd = zeros(L_TT_P , 1) ; 
    for tidx =1 : L_TT_P % paralellizable
        avpd(tidx) = tt_update{tidx}.pD ; % using the legacy part pD (still the same for other parts)
    end
    avqd = 1 - avpd ; 
    log_avpd = log(avpd) ; 
        % cardinality distribution
    avps = ones(3*numprevTT , 1) ;
    avqs = zeros(3*numprevTT , 1) ;
    for tabidx =1 : numprevTT
        tidx = (tabidx-1)*3+1 ;  
        pSp = tt_in{tabidx}.mode(2) ;
        avqs(tidx)   = (1-pSp)*model.card_dist_mode(1,1) +  pSp*model.card_dist_mode(2,1); 
        avps(tidx)   = (1-pSp)*model.card_dist_mode(1,2) +  pSp*model.card_dist_mode(2,2); 
        avps(tidx+1) = (1-pSp)*model.card_dist_mode(1,3) +  pSp*model.card_dist_mode(2,3); 
    end
    avqs = [1-r_b ; avqs] ; 
    avps = [r_b  ; avps] ; 
    log_avps = log(avps) ; 
    log_avqs = log(avqs) ;
    
    % gated measurement  index matrix
    gatemeasidxs= zeros(length(tt_predict),m);
    for tabidx= 1:length(tt_predict)
        gatemeasidxs(tabidx,1:length(tt_predict{tabidx}.gatemeas))= tt_predict{tabidx}.gatemeas;
    end
    gatemeasindc= gatemeasidxs>0;

    %% Sampling for components
    cpreds= L_TT_P;
    nbirths= numBirths;
    nexists= length(tt_in);
    ntracks= nbirths + nexists*3;
    temp_tindices = 1:nexists*3 ;
    if model.unknown_clutter
        temp_clutter_rate = (N_clutter_in*model.clutter_P_S + model.clutter_gen*model.clutter_P_B)*model.clutter_P_D*model.pdf_c ; 
    else
        temp_clutter_rate = model.lambda_c*model.pdf_c ; 
    end
    tindices= [1:nbirths nbirths+temp_tindices];
    lselmask= false(length(tt_predict),m); lselmask(tindices,:)= gatemeasindc(tindices,:);
    mindices= unique_faster(gatemeasidxs(lselmask)); 
    lindices= length(tindices) ; 
    diag_indc = 1:lindices+1:lindices*lindices;
    
    dead_mat = -Inf*ones(lindices) ; 
    dead_mat(diag_indc) = log_avqs(tindices) ; 
    miss_mat = -Inf*ones(lindices) ; 
    miss_mat(diag_indc) = log_avps(tindices) + log_avpd(tindices) ; 
    log_const_jointcostm = [dead_mat, miss_mat] ;
    logcostm = [log_const_jointcostm ...
          repmat(log_avps(tindices),[1 length(mindices)]) + log_allcostm(tindices,mindices)] ;
    logcostm(:,2*length(tindices)+1:end) = logcostm(:,2*length(tindices)+1:end) - log(temp_clutter_rate) ; 
    neglogcostm= -logcostm;                         
    [uasses]= block_gibbs(neglogcostm,ceil(filter.H_upd*sampling_factor),nbirths); 
    uasses(uasses<=ntracks)= -Inf;
    uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;
    uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;   
    uasses(uasses>0)= mindices(uasses(uasses>0));
    %% Generating children hypotheses  
    hypo = cell(size(uasses,1),1) ; 
    for hidx=1:size(uasses,1)
        update_hypcmp_tmp= uasses(hidx,:)'; 
        update_hypcmp_idx=update_hypcmp_tmp ;
        update_hypcmp_idx= cpreds.*update_hypcmp_idx+tindices';
        % Construct the hypo from here
        hypo{hidx}.meas_ass = update_hypcmp_tmp ; 
            % select prediction model
        pred_mode_cost = zeros(nexists, length(model.ws)+length(model.wsp)) ; 
        log_pE = 0 ; % store log of existence probability
        n = 0 ; 
        % for births (dont need to choose models for births)
        nbirths = 0 ; 
        for tidx = 1 : numBirths
            if update_hypcmp_idx(tidx)>=0
                log_pE = log_pE + log_avps(tidx) ; 
                nbirths = nbirths + 1 ;
                n = n + 1 ;
            elseif update_hypcmp_idx(tidx)<0
                log_pE = log_pE + log_avqs(tidx) ; 
            end
        end
        % for survive + spawn objects 
        nspawns = 0 ;
        otidx = 1 ; % old track idx
        dead_indc = false(1, (ntracks-numBirths)/3) ; 
        for tidx = numBirths+1 : 3 : ntracks
            if update_hypcmp_idx(tidx)>=0 && update_hypcmp_idx(tidx+1)< 0
                temp = tt_update{update_hypcmp_idx(tidx)}.w ; 
                temp = reshape(temp, length(model.ws), []) ; temp = temp' ; 
                pred_mode_cost(otidx, 1:length(model.ws)) = sum(temp,1) ; 
                log_pE = log_pE + log_avps(tidx) ; 
                n = n + 1 ; 
            elseif update_hypcmp_idx(tidx)<0  && update_hypcmp_idx(tidx+1)>= 0
                temp = tt_update{update_hypcmp_idx(tidx+1)}.w.*tt_update{update_hypcmp_idx(tidx+2)}.w ; 
                temp = reshape(temp, length(model.wsp), []) ; temp = temp' ; 
                pred_mode_cost(otidx, length(model.ws)+1 : end) = sum(temp,1) ; 
                log_pE = log_pE +log_avps(tidx+1) ; 
                nspawns = nspawns + 1 ; 
                n = n + 2 ; 
            elseif update_hypcmp_idx(tidx)<0 && update_hypcmp_idx(tidx+1)<0 ...
                    && update_hypcmp_idx(tidx+2)<0
                log_pE = log_pE + log_avqs(tidx) ; 
                dead_indc(otidx) = true ; 
            else
                keyboard ; 
            end    
            otidx = otidx + 1 ; 
        end
        % run k-shortest path to generate combination of models
        pred_mode_cost(dead_indc,:) = [] ; 
        nexists_now = nexists - sum(dead_indc) ; 
        num_models = 1 ; 
        neglog_pred_cost = -log(pred_mode_cost) + eps ; 
        if size(neglog_pred_cost,1)>0
            [spath, ~] = kshortestwrap_pred_model(neglog_pred_cost, num_models) ; 
        else
            spath =cell(0,1) ; 
        end
        % post-processing spath to form pred_model
        pred_model = cell(length(spath),1) ; 
        for pidx = 1 : length(spath)
            if ~isempty(spath{pidx})
                pred_model{pidx} = zeros(2, nexists_now) ; 
                for tidx = 1 : nexists_now
                    if spath{pidx}(tidx)<=length(model.ws)
                        pred_model{pidx}(1,tidx) = spath{pidx}(tidx) ;  
                        pred_model{pidx}(2,tidx) = 1 ; % only 1 shifting survival (that is no shift!)
                    else
                        pred_model{pidx}(1,tidx) = -1 ; % only 1 spawning model
                        pred_model{pidx}(2,tidx) = spath{pidx}(tidx)-length(model.ws) ; 
                    end
                end
            end
        end
        hypo{hidx}.pred_model = pred_model ; 
        % grab the pD and s,t and mode so that I dont need to do it again in the
        % exact propagation
        existing_tt = update_hypcmp_idx(update_hypcmp_idx>0) ; 
        LE = length(existing_tt) ;
        pD = ones(LE, 1) ; % I will return this pD later
        pD_qD = ones(LE, 1) ; % store the pD or qD for weight calculation
        qz_mode = ones(LE, 1) ; 
        st = zeros(LE, 2) ; 
        mode = zeros(LE, 2) ; 
        labels = cell(LE,1) ; 
        for eidx = 1 : LE 
            pD(eidx) = tt_update{existing_tt(eidx)}.pD ; 
            if existing_tt(eidx)<=cpreds
                pD_qD(eidx) = 1-pD(eidx) ; 
            else
                pD_qD(eidx) = pD(eidx) ; 
            end
            qz_mode(eidx) = tt_update{existing_tt(eidx)}.qz_mode ; 
            st(eidx,:) = tt_update{existing_tt(eidx)}.st ;
            mode(eidx,:) = tt_update{existing_tt(eidx)}.mode ;
            labels{eidx} = tt_update{existing_tt(eidx)}.l ;
        end
        hypo{hidx}.pD = pD ; 
        hypo{hidx}.st = st ; 
        hypo{hidx}.mode = mode ;
        hypo{hidx}.log_pE_weight = log_pE ; 
        hypo{hidx}.n = n ; 
        hypo{hidx}.nbirths = nbirths ; 
        hypo{hidx}.nspawns = nspawns ;  
        hypo{hidx}.log_qD_pD_qz_mode_weight = sum(log(pD_qD)) + sum(log(qz_mode)); % log of everything that is independent among objects
        hypo{hidx}.l = labels ; 
    end
end