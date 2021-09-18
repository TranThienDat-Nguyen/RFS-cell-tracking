function [glmb_nextupdate, meas_ru] = update_approximation(model , filter , glmb_in , tt_birth , ~ , meas , k,clutter_table)
%% Predict track tables
    % create birth tracks 
    nbirths = length(tt_birth) ;
    r_b_normal = zeros(nbirths , 1) ;
    for bidx = 1 : nbirths
        r_b_normal(bidx)= tt_birth{bidx}.r_b  ; 
    end

    % create survivals and spawned tracks
    numprevTT = length(glmb_in.tt) ;
    tt_surv = cell(numprevTT , 1) ;
    tt_sp1 = cell(numprevTT , 1) ;
    tt_sp2 = cell(numprevTT , 1) ;
    
    tt_in = glmb_in.tt ; 
    parfor tabidx = 1 : numprevTT
        % Survivals
        % first predict the parameters for pD beta
        ss = tt_in{tabidx}.st(1) ; 
        tt = tt_in{tabidx}.st(2) ; 
        beta_mean = ss / (ss + tt);
        beta_var = (ss*tt) / ((ss+tt)^2*(ss+tt+1)); 
        tt_surv{tabidx}.st = zeros(1,2) ;  
        tt_surv{tabidx}.st(1) = ( ( beta_mean * (1-beta_mean)/ beta_var) - 1 ) * beta_mean;
        tt_surv{tabidx}.st(2) =( ( beta_mean * (1-beta_mean)/ beta_var - 1) ) * (1 - beta_mean); 
        [wtemp_predict , mtemp_predict , Ptemp_predict] = model.kalman_predict_multi_object(model , tt_in{tabidx}.w , tt_in{tabidx}.m , tt_in{tabidx}.P , [] , true) ; 
        [wtemp_predict , mtemp_predict , Ptemp_predict] = gaus_cap(wtemp_predict ,...
                            mtemp_predict ,Ptemp_predict , filter.L_max) ; 
        [wtemp_predict , mtemp_predict , Ptemp_predict] = gaus_merge(wtemp_predict , mtemp_predict , Ptemp_predict, filter.merge_threshold) ; 
        [wtemp_predict , mtemp_predict , Ptemp_predict] = gaus_prune(wtemp_predict , mtemp_predict , Ptemp_predict, filter.elim_threshold) ;
        tt_surv{tabidx}.w= wtemp_predict/sum(wtemp_predict);
        tt_surv{tabidx}.m= mtemp_predict;                                                                                   
        tt_surv{tabidx}.P= Ptemp_predict;   
        tt_surv{tabidx}.l= tt_in{tabidx}.l ; 
        tt_surv{tabidx}.ah = tt_in{tabidx}.ah ;
        tt_surv{tabidx}.mode = [(1-model.pSp) , model.pSp] ;
        
        % First spawns
        tt_sp1{tabidx}.st = model.init_st ;
        [wtemp_predict , mtemp_predict , Ptemp_predict] = model.kalman_predict_multi_object(model , tt_in{tabidx}.w , tt_in{tabidx}.m , tt_in{tabidx}.P , 1 , true) ; 
        tt_sp1{tabidx}.w = wtemp_predict/sum(wtemp_predict) ; 
        tt_sp1{tabidx}.m= mtemp_predict{1};
        tt_sp1{tabidx}.P= Ptemp_predict{1};        
        tt_sp1{tabidx}.l= [tt_in{tabidx}.l;k;1] ;  
        tt_sp1{tabidx}.ah = [] ; 
        tt_sp1{tabidx}.mode = [(1-model.pSp) , model.pSp] ;
        
        % Second spawns
        tt_sp2{tabidx}.st = model.init_st;
        tt_sp2{tabidx}.w = wtemp_predict/sum(wtemp_predict) ;
        tt_sp2{tabidx}.m= mtemp_predict{2};
        tt_sp2{tabidx}.P= Ptemp_predict{2};
        tt_sp2{tabidx}.l= [tt_in{tabidx}.l;k;2] ; 
        tt_sp2{tabidx}.ah = [] ;    
        tt_sp2{tabidx}.mode = [(1-model.pSp) , model.pSp] ;

    end

    % restructuring tt_pred
    tt_predict = cell(numprevTT*3 + nbirths , 1) ;
    tt_predict(1:nbirths) = tt_birth ; 
    for tabidx = 1 : numprevTT
        tabssidx = nbirths+ (tabidx-1)*3+1 ; 
        tt_predict(tabssidx) = tt_surv(tabidx) ; 
        tt_predict(tabssidx+1) = tt_sp1(tabidx) ; 
        tt_predict(tabssidx+2) = tt_sp2(tabidx) ;
    end
    
%% Update track tables
    %gating by tracks
    if filter.gate_flag
        for tabidx=1:length(tt_predict)
            tt_predict{tabidx}.gatemeas= gate_meas_gms_idx(meas.Z{k},filter.gamma,model,tt_predict{tabidx}.m,tt_predict{tabidx}.P);
        end
    else
        for tabidx=1:length(tt_predict)
            tt_predict{tabidx}.gatemeas= 1:size(meas.Z{k},2);
        end
    end
    % pre-calcs for the mode stuff
    m= size(meas.Z{k},2); 
    L_TT_P = length(tt_predict) ; 
    
    %create updated tracks (single target Bayes update)
    tt_update= cell(length(tt_predict)*(m+1),1) ; 
    
        %missed detection tracks (legacy tracks)
    for tabidx= 1:L_TT_P
        tt_update{tabidx}= tt_predict{tabidx};
        % Cap the Gaussian of the miss detected tracks
        [w_temp , m_temp , P_temp] = gaus_cap(tt_update{tabidx}.w ,...
                            tt_update{tabidx}.m ,tt_update{tabidx}.P , 100) ; 
        [w_temp , m_temp , P_temp] = gaus_merge(w_temp , m_temp , P_temp, filter.merge_threshold) ; 
        [w_temp , m_temp , P_temp] = gaus_prune(w_temp , m_temp , P_temp, 1e-4) ;
        tt_update{tabidx}.w = w_temp ; 
        tt_update{tabidx}.m = m_temp ; 
        tt_update{tabidx}.P = P_temp ; 
        tt_update{tabidx}.ah= [tt_update{tabidx}.ah; 0];
        % grab pD
        if model.unknown_pD
            tt_update{tabidx}.pD = tt_predict{tabidx}.st(1)/(tt_predict{tabidx}.st(1) + tt_predict{tabidx}.st(2));
        else
            tt_update{tabidx}.pD = model.pD ;
        end
        % update s,t
        tt_update{tabidx}.st(2) = tt_update{tabidx}.st(2)+1 ; 
        tt_update{tabidx}.logqz = 0 ; 
    end

         % measurement updated tracks (all pairs)
    
    log_allcostm_cell = cell(1,L_TT_P) ; 
    tt_update_para = cell(1,L_TT_P) ; 
    parfor tabidx= 1:L_TT_P
        tt_update_para{tabidx} = cell(1,m) ; 
        if model.unknown_pD
            curr_pD = tt_predict{tabidx}.st(1)/(tt_predict{tabidx}.st(1) + tt_predict{tabidx}.st(2));
        else
            curr_pD = model.pD ; 
        end
        log_allcostm_cell{tabidx} = log(zeros(1,m)) ; 
        for mmidx = 1 : length(tt_predict{tabidx}.gatemeas)
            emm= tt_predict{tabidx}.gatemeas(mmidx) ; 
            % detection state
            tt_update_para{tabidx}{emm}.pD = curr_pD ;  
            tt_update_para{tabidx}{emm}.st = tt_predict{tabidx}.st ; 
            tt_update_para{tabidx}{emm}.st(1) = tt_update_para{tabidx}{emm}.st(1) + 1 ; 
            % kinematic state
            [w_temp,m_temp,P_temp] = kalman_update_single_object(meas.Z{k}(:,emm),model, ...
                tt_predict{tabidx}.w, tt_predict{tabidx}.m, tt_predict{tabidx}.P);   %kalman update for this track and this measurement
            log_allcostm_cell{tabidx}(emm) = log(sum(w_temp)) + log(curr_pD);
            w_temp = w_temp/sum(w_temp) ; 
            [w_temp , m_temp , P_temp] = gaus_cap(w_temp , m_temp ,P_temp , filter.L_max) ; 
            [w_temp , m_temp , P_temp] = gaus_merge(w_temp , m_temp , P_temp, filter.merge_threshold) ; 
            [w_temp , m_temp , P_temp] = gaus_prune(w_temp , m_temp , P_temp, filter.elim_threshold) ;
            w_temp= w_temp+eps;  
            tt_update_para{tabidx}{emm}.w = w_temp ; 
            tt_update_para{tabidx}{emm}.m = m_temp ; 
            tt_update_para{tabidx}{emm}.P = P_temp ;
            % label
            tt_update_para{tabidx}{emm}.l = tt_predict{tabidx}.l ;
            % update the mode
            g_normal = meas.Z_normal{k}(emm) ; g_spawn = meas.Z_spawn{k}(emm) ;
            temp = tt_predict{tabidx}.mode .* [g_normal g_spawn]; 
            tt_update_para{tabidx}{emm}.mode = temp/sum(temp);
            tt_update_para{tabidx}{emm}.ah= [tt_predict{tabidx}.ah; emm];
            % store the likelihood for later use
            tt_update_para{tabidx}{emm}.logqz = log_allcostm_cell{tabidx}(emm)-log(curr_pD) ; 
        end
    end
    % convert to standard tt format
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
    avpd = zeros(length(tt_predict) , 1) ; 
    parfor tidx =1 : length(tt_predict)
        avpd(tidx) = tt_update{tidx}.pD ;
    end
    avqd = 1 - avpd ; 
    log_avpd = log(avpd) ; 
    log_avqd = log(avqd) ;
        % cardinality distribution
    avps = ones(3*numprevTT , 1) ;
    avqs = zeros(3*numprevTT , 1) ;
    for tabidx =1 : numprevTT
        tidx = (tabidx-1)*3+1 ;  
        %   real cost
        pSp = tt_in{tabidx}.mode(2) ;
        avqs(tidx)   = (1-pSp)*model.card_dist_mode(1,1) +  pSp*model.card_dist_mode(2,1); 
        avps(tidx)   = (1-pSp)*model.card_dist_mode(1,2) +  pSp*model.card_dist_mode(2,2); 
        avps(tidx+1) = (1-pSp)*model.card_dist_mode(1,3) +  pSp*model.card_dist_mode(2,3); 
    end
    avqs = [1-r_b_normal ; avqs] ; 
    avps = [r_b_normal  ; avps] ; 
    log_avps = log(avps) ; 
    log_avqs = log(avqs) ;

    % gated measurement  index matrix
    gatemeasidxs= zeros(length(tt_predict),m);
    for tabidx= 1:length(tt_predict)
        gatemeasidxs(tabidx,1:length(tt_predict{tabidx}.gatemeas))= tt_predict{tabidx}.gatemeas;
    end
    gatemeasindc= gatemeasidxs>0;
    
    %% Update multi-object prior components
    num_comp = length(glmb_in.w) ; 
    % initilization for parallelization
    glmb_nextupdate_w = cell(1,num_comp); 
    glmb_nextupdate_N_clutter = cell(1,num_comp) ; 
    glmb_nextupdate_N_clutter_count = cell(1,num_comp); 
    glmb_nextupdate_I = cell(1,num_comp);                                                                                              %hypothesis/component tracks (via indices to track table)
    glmb_nextupdate_I_old = cell(1,num_comp) ;
    glmb_nextupdate_n = cell(1,num_comp);   
    glmb_nextupdate_total_births = cell(1,num_comp) ;
    glmb_nextupdate_total_spawns = cell(1,num_comp) ;
    meas_ru_ind_cell = cell(1,num_comp) ; 
        
    % copy to small arrays to reduce parallel overhead
    num_new_comp = zeros(1, length(glmb_in.w)) ; 
    I_in = glmb_in.I ; 
    w_in = glmb_in.w ; 
    N_clutter_in = glmb_in.N_clutter ; 
    tt_update_new = cell(2000,1) ; % store new tt_update (initialize 2000 cells every time, avoid unknown size cell array)
    curr_max_tt_new = 2000 ; 
    tt_update_new_info = zeros(0,7) ; % store assignement info of new tt_udpate 
    L_tt_update = length(tt_update) ; 
    tt_idx = 1 ;  % idx to the tt_update_new
    cpreds= L_TT_P;
    ch_hypo = cell(num_comp,1) ; % store children components
    tindices_cell = cell(num_comp,1) ; 
    parfor pidx=1:num_comp % paralell
        nexists= length(I_in{pidx});
        ntracks= nbirths + nexists*3;
        %calculate best updated components
        temp_tindices = zeros(1 , nexists*3) ; 
        eidx = 1 ; 
        if model.unknown_clutter
            temp_clutter_rate = (N_clutter_in(pidx)*model.clutter_P_S + ...
                model.clutter_gen*model.clutter_P_B)*model.clutter_P_D*model.pdf_c ; 
        else
            temp_clutter_rate = model.lambda_c*model.pdf_c ;
        end
            
        for ttidx =1 : 3 : length(temp_tindices)
            temp_tindices(ttidx) = (I_in{pidx}(eidx)-1)*3 + 1 ; 
            temp_tindices(ttidx+1) = temp_tindices(ttidx) +1 ;
            temp_tindices(ttidx+2) = temp_tindices(ttidx) +2 ;
            eidx = eidx + 1 ;
        end
        tindices= [1:nbirths nbirths+temp_tindices];
        tindices_cell{pidx} = tindices ; 
        lselmask= false(length(tt_predict),m); lselmask(tindices,:)= gatemeasindc(tindices,:);
        mindices= unique_faster(gatemeasidxs(lselmask)); 
        lindices = length(tindices) ; 
        diag_indc = 1:lindices+1:lindices*lindices;
        dead_mat = -Inf*ones(lindices) ; 
        dead_mat(diag_indc) = log_avqs(tindices) ; 
        miss_mat = -Inf*ones(lindices) ; 
        miss_mat(diag_indc) = log_avps(tindices) + log_avqd(tindices) ; 
        log_const_jointcostm = [dead_mat, miss_mat] ;
        logcostm = [log_const_jointcostm ...
          repmat(log_avps(tindices),[1 length(mindices)]) + log_allcostm(tindices,mindices)] ;
        logcostm(:,2*length(tindices)+1:end) = logcostm(:,2*length(tindices)+1:end) - log(temp_clutter_rate) ; 
        neglogcostm= -logcostm;
        [uasses]= block_gibbs(neglogcostm,round(filter.H_upd*sqrt(w_in(pidx))/sum(sqrt(w_in))),nbirths); 
        uasses(uasses<=ntracks)= -Inf;
        uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;
        uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;
        uasses(uasses>0)= mindices(uasses(uasses>0));
        ch_hypo{pidx} = uasses; 
    end

    
    for pidx=1:num_comp 
        runidx =1 ; 
        uasses = ch_hypo{pidx} ; 
        tindices = tindices_cell{pidx} ; 
        for hidx=1:size(uasses,1)
            curr_I = [] ; % current hypothesis components
            update_hypcmp_tmp= uasses(hidx,:)'; 
            update_hypcmp_idx=update_hypcmp_tmp ;
            update_hypcmp_idx= cpreds.*update_hypcmp_idx+tindices';
            % Grab the stats and indices to calcualte the correct hypothesis weight
            newBirths = sum(update_hypcmp_idx(1:nbirths)>0) ; 
            numSpawns = 0 ; 
            log_w = 0 ; 
            utidx =1 ; 
            while utidx < length(tindices)+1 % go through tracks
                curr_track = tindices(utidx) ; 
                curr_asgn = update_hypcmp_tmp(utidx) ; 
                notBirth = curr_track > nbirths ;
                % check conditions
                birth_flag = false ; 
                dead_flag = false ; 
                survive_flag = false ; 
                spawn_flag = false ;
                if ~notBirth && ~(curr_asgn<0)
                    birth_flag = true ; 
                end
                if notBirth && ~(curr_asgn<0) && update_hypcmp_idx(utidx+1)<0
                    survive_flag = true ; 
                end
                if (notBirth && curr_asgn<0 && update_hypcmp_idx(utidx+1)<0) || ...
                    (~notBirth && curr_asgn<0)
                    dead_flag = true ;
                end
                if (notBirth && curr_asgn<0) && ~(update_hypcmp_idx(utidx+1)<0)
                    spawn_flag = true ;
                end
                    
                % Construct the tt_update_new and compute the hypothesis
                % weights
                if dead_flag
                    log_w = log_w + log(avqs(curr_track)) ;
                elseif survive_flag || birth_flag
                    if curr_asgn>0
                       log_w = log_w + log_avps(curr_track) + log_avpd(curr_track) + ...
                           tt_update{update_hypcmp_idx(utidx)}.logqz ; 
                    else
                        log_w = log_w + log_avps(curr_track) + log_avqd(curr_track) ;
                    end
                    % Use the new index scheme for curr_I
                    curr_I = [curr_I , update_hypcmp_idx(utidx)] ; 
                elseif spawn_flag
                    numSpawns = numSpawns + 1 ;
                    sp1_idx = tindices(utidx+1) ; 
                    sp2_idx = tindices(utidx+2) ;
                    sp1_meas_idx = update_hypcmp_tmp(utidx+1) ; 
                    sp2_meas_idx = update_hypcmp_tmp(utidx+2) ; 
                    
                    if sp1_meas_idx>0 && sp2_meas_idx>0
                        temp_vec = -ones(1,7) ; 
                        temp_vec(1) = runidx ; 
                        temp_vec(2) = tt_idx ; 
                        temp_vec(3) = sp1_idx ; 
                        temp_vec(4) = sp1_meas_idx ; 
                        temp_vec(5) = tt_idx + 1 ; 
                        temp_vec(6) = sp2_idx ; 
                        temp_vec(7) = sp2_meas_idx ; 
                        temp_idx = find_vec(tt_update_new_info(:,[3 4]), temp_vec(:,[3 4])) ;  
                        if ~isempty(temp_idx)
                            curr_I = [curr_I, tt_update_new_info(temp_idx, 2)+L_tt_update , ...
                                tt_update_new_info(temp_idx, 5)+L_tt_update] ; % I will concat tt_update_new to tt_udapte later
                        else
                            tt_update_new_info = [tt_update_new_info; temp_vec] ; 
                            z = zeros(2*model.z_dim, 1) ; 
                            z(1:model.z_dim) = meas.Z{k}(:,sp1_meas_idx) ; 
                            z(model.z_dim+1:end) = meas.Z{k}(:,sp2_meas_idx) ;
                            z_norm_spawn = zeros(2) ; 
                            z_norm_spawn(1,1) = meas.Z_normal{k}(sp1_meas_idx) ; 
                            z_norm_spawn(1,2) = meas.Z_spawn{k}(sp1_meas_idx) ;
                            z_norm_spawn(2,1) = meas.Z_normal{k}(sp2_meas_idx) ; 
                            z_norm_spawn(2,2) = meas.Z_spawn{k}(sp2_meas_idx) ;
                            ah = [sp1_meas_idx, sp2_meas_idx] ; 
                            [tt_sp1 , tt_sp2] = kalman_update_multi_object(model, z, z_norm_spawn, ah, tt_predict{sp1_idx}, tt_predict{sp2_idx}) ; 
                            tt_update_new{tt_idx}= tt_sp1 ; 
                            tt_update_new{tt_idx+1}= tt_sp2 ; 
                            curr_I = [curr_I, L_tt_update+tt_idx, L_tt_update+tt_idx+1];   
                        end
                    else
                        curr_I = [curr_I, update_hypcmp_idx(utidx+1), update_hypcmp_idx(utidx+2)] ; 
                    end

                    % spawn weight
                    log_w = log_w + log_avps(sp1_idx) ; 
                    % spawn 1
                    if sp1_meas_idx>0
                        log_w = log_w + log_avpd(sp1_idx) ; 
                    else
                        log_w = log_w + log_avqd(sp1_idx) ; 
                    end
                    % spawn 2
                    if sp2_meas_idx>0
                        log_w = log_w + log_avpd(sp2_idx) ; 
                    else
                        log_w = log_w + log_avqd(sp2_idx) ; 
                    end
                    % likelihood
                    if sp1_meas_idx>0 && sp2_meas_idx>0 && isempty(temp_idx)
                        log_w = log_w + tt_update_new{tt_idx}.logqz + tt_update_new{tt_idx+1}.logqz ; 
                        
                        tt_idx = tt_idx + 2 ;
                        if tt_idx>curr_max_tt_new
                            curr_max_tt_new = curr_max_tt_new + 2000 ; 
                            tt_temp = cell(2000,1) ; 
                            tt_update_new = [tt_update_new ; tt_temp];
                        end
                    else
                        log_w = log_w + tt_update{update_hypcmp_idx(utidx+1)}.logqz + ...
                            tt_update{update_hypcmp_idx(utidx+2)}.logqz ;
                    end
                end
                if ~notBirth % if birth just increase by 1
                    utidx = utidx + 1 ;
                else
                    utidx = utidx + 3 ; 
                end
            end
            % Solve for the clutter
            num_Z_0 = m - sum(update_hypcmp_tmp>0) ; 
            scidx = N_clutter_in(pidx) +1 ;
            N_top_this = sum(clutter_table{scidx , num_Z_0+1}(1,:)>0) ; 
            varClutterWeight = 1;
            N_top_this = ceil(min(N_top_this , filter.H_clutter) * varClutterWeight) ; 
            hyps_clutter = clutter_table{scidx , num_Z_0+1}(: , 1 : N_top_this) ; 
            meas_ru_ind_temp = false(m,1) ; 
            meas_ru_ind_temp(uasses(hidx,(uasses(hidx,:)>0))) = true ;
            if model.unknown_clutter
                for chidx =1 : size(hyps_clutter , 2)
                    glmb_nextupdate_w{pidx}(runidx) = log_w  + hyps_clutter(2,chidx) + log(w_in(pidx)); 
                    glmb_nextupdate_N_clutter{pidx}(runidx) = hyps_clutter(1,chidx)-1 ; 
                    glmb_nextupdate_N_clutter_count{pidx}(runidx) = m - sum(uasses(hidx,:)>0) ; 
                    glmb_nextupdate_I{pidx}{runidx}= curr_I;                                              
                    glmb_nextupdate_I_old{pidx}(runidx) = pidx ;
                    glmb_nextupdate_n{pidx}(runidx)= sum(update_hypcmp_idx>0);  
                    glmb_nextupdate_total_births{pidx}(runidx) = newBirths ;
                    glmb_nextupdate_total_spawns{pidx}(runidx) = numSpawns ;
                    meas_ru_ind_cell{pidx}(:,runidx) = meas_ru_ind_temp ; 
                    runidx = runidx  + 1 ;  
                end
            else
                glmb_nextupdate_w{pidx}(runidx) = log_w - model.lambda_c + num_Z_0*log(model.lambda_c*model.pdf_c)  + log(w_in(pidx)); 
                glmb_nextupdate_N_clutter{pidx}(runidx) = model.lambda_c ; 
                glmb_nextupdate_N_clutter_count{pidx}(runidx) = num_Z_0 ; 
                glmb_nextupdate_I{pidx}{runidx}= curr_I;                                                                     
                glmb_nextupdate_I_old{pidx}(runidx) = pidx ;
                glmb_nextupdate_n{pidx}(runidx)= sum(update_hypcmp_idx>0);  

                glmb_nextupdate_total_births{pidx}(runidx) = newBirths ;
                glmb_nextupdate_total_spawns{pidx}(runidx) = numSpawns ;
                meas_ru_ind_cell{pidx}(:,runidx) = meas_ru_ind_temp ; 
                runidx = runidx  + 1 ; 
            end
                
        end
        num_new_comp(pidx) = runidx - 1 ;
    end
    if tt_idx == 1 
        tt_update_new = [] ;
    else
        tt_update_new(tt_idx:end) = [] ;
    end
    % re-structuring the glmb_nextupdate
    stidx =1  ; 
    meas_ru_ind = zeros(m, sum(num_new_comp)) ;
    for pidx = 1 : num_comp
        glmb_nextupdate.w(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_w{pidx} ; 
        glmb_nextupdate.N_clutter(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_N_clutter{pidx} ; 
        glmb_nextupdate.N_clutter_count(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_N_clutter_count{pidx} ; 
        glmb_nextupdate.I(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_I{pidx} ; 
        glmb_nextupdate.I_old(stidx:stidx+num_new_comp(pidx)-1) = pidx ; 
        glmb_nextupdate.n(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_n{pidx} ; 
        glmb_nextupdate.total_births(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_total_births{pidx} ; 
        glmb_nextupdate.total_spawns(stidx:stidx+num_new_comp(pidx)-1) = glmb_nextupdate_total_spawns{pidx} ; 
        meas_ru_ind(:, stidx:stidx+num_new_comp(pidx)-1) = meas_ru_ind_cell{pidx} ; 
        stidx = stidx + num_new_comp(pidx) ; 
    end
    glmb_nextupdate.tt = [tt_update; tt_update_new] ; 
    glmb_nextupdate.w= exp(glmb_nextupdate.w-logsumexp(glmb_nextupdate.w,[],2));     
    meas_ru = meas_ru_ind * glmb_nextupdate.w' ; 
    %normalize weights
    for card=0:max(glmb_nextupdate.n)
        glmb_nextupdate.cdn(card+1)= sum(glmb_nextupdate.w(glmb_nextupdate.n==card));                                                                                                       %extract probability of n targets
    end
    for card_spawn = 0:max(glmb_nextupdate.total_spawns)
        glmb_nextupdate.cdn_spawn(card_spawn+1) = sum(glmb_nextupdate.w(glmb_nextupdate.total_spawns==card));  
    end    
    % prune and cap
    glmb_nextupdate = prune(glmb_nextupdate , filter) ; 
    glmb_nextupdate = cap(glmb_nextupdate,filter) ;  
    glmb_nextupdate= clean_update(glmb_nextupdate);
    %% Estimate tracks for each hypothesis
    for hidx = 1 : length(glmb_nextupdate.w)
        est_X_pD = zeros(model.x_dim+1, glmb_nextupdate.n(hidx)) ; 
        est_L = cell(1, glmb_nextupdate.n(hidx)) ; 
        glmb_nextupdate.est_X_pD{hidx} = glmb_in.est_X_pD{glmb_nextupdate.I_old(hidx)} ; 
        glmb_nextupdate.est_L{hidx} = glmb_in.est_L{glmb_nextupdate.I_old(hidx)} ; 
        glmb_nextupdate.est_others{hidx} = glmb_in.est_others{glmb_nextupdate.I_old(hidx)} ;
        for n = 1 : glmb_nextupdate.n(hidx)
            [~,idxcmp] = max(glmb_nextupdate.tt{glmb_nextupdate.I{hidx}(n)}.w) ; 
            est_X_pD(:, n) = [glmb_nextupdate.tt{glmb_nextupdate.I{hidx}(n)}.m(:,idxcmp) ; ...
                glmb_nextupdate.tt{glmb_nextupdate.I{hidx}(n)}.pD]; 
            est_L{n} = glmb_nextupdate.tt{glmb_nextupdate.I{hidx}(n)}.l ; 
        end
        glmb_nextupdate.est_X_pD{hidx}{k} = est_X_pD ;
        glmb_nextupdate.est_L{hidx}{k} = est_L ;
        % other estimates stuff ( [n; n_birth; n_spawn; clutter_rate; clutter_count] )
        glmb_nextupdate.est_others{hidx}(:,k) = [glmb_nextupdate.n(hidx) ; 
                                                 glmb_nextupdate.total_births(hidx) ; 
                                                 glmb_nextupdate.total_spawns(hidx) ;
                                                 glmb_nextupdate.N_clutter(hidx) ;
                                                 glmb_nextupdate.N_clutter_count(hidx)]  ;
    end
end

    
function glmb_clean= clean_update(glmb_temp)
    %flag used tracks
    usedindicator= zeros(length(glmb_temp.tt),1);
    for hidx= 1:length(glmb_temp.w)
        usedindicator(glmb_temp.I{hidx})= usedindicator(glmb_temp.I{hidx})+1;
    end
    trackcount= sum(usedindicator>0);

    %remove unused tracks and reindex existing hypotheses/components
    newindices= zeros(length(glmb_temp.tt),1); newindices(usedindicator>0)= 1:trackcount;
    glmb_clean.tt= glmb_temp.tt(usedindicator>0);
    glmb_clean.w= glmb_temp.w;
    for hidx= 1:length(glmb_temp.w)
        glmb_clean.I{hidx}= newindices(glmb_temp.I{hidx});
    end
    glmb_clean.n= glmb_temp.n;
    glmb_clean.cdn= glmb_temp.cdn;
    glmb_clean.cdn_spawn= glmb_temp.cdn_spawn;
    glmb_clean.total_births= glmb_temp.total_births;
    glmb_clean.total_spawns= glmb_temp.total_spawns;
    glmb_clean.N_clutter= glmb_temp.N_clutter;
    glmb_clean.N_clutter_count= glmb_temp.N_clutter_count;
    glmb_clean.I_old = glmb_temp.I_old ; 
end


function glmb_out= prune(glmb_in,filter)
    %prune components with weights lower than specified threshold
    idxkeep= find(glmb_in.w > filter.hyp_threshold);
    if length(idxkeep)<filter.H_min
        [~,idxkeep] = sort(glmb_in.w, 'descend') ; 
        idxkeep = idxkeep(1:filter.H_min) ; 
    end
    glmb_out.tt= glmb_in.tt;
    glmb_out.w= glmb_in.w(idxkeep);
    glmb_out.I= glmb_in.I(idxkeep);
    glmb_out.n= glmb_in.n(idxkeep);
    glmb_out.total_spawns = glmb_in.total_spawns(idxkeep) ; 
    glmb_out.total_births = glmb_in.total_births(idxkeep) ; 
    glmb_out.N_clutter = glmb_in.N_clutter(idxkeep) ;
    glmb_out.N_clutter_count= glmb_in.N_clutter_count(idxkeep);
    glmb_out.w= glmb_out.w/sum(glmb_out.w);
    for card=0:max(glmb_out.n)
        glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
    end
    for card = 0:max(glmb_out.total_spawns)
        glmb_out.cdn_spawn(card+1) = sum(glmb_out.w(glmb_out.total_spawns==card));
    end
    glmb_out.I_old = glmb_in.I_old(idxkeep) ; 
end



function glmb_out= cap(glmb_in,filter)
%cap total number of components to specified maximum
    if length(glmb_in.w) > filter.H_max
        [~,idxsort]= sort(glmb_in.w,'descend');
        idxkeep=idxsort(1:filter.H_max);
        glmb_out.tt= glmb_in.tt;
        glmb_out.w= glmb_in.w(idxkeep);
        glmb_out.I= glmb_in.I(idxkeep);
        glmb_out.I_old= glmb_in.I_old(idxkeep);
        glmb_out.n= glmb_in.n(idxkeep);
        glmb_out.n= glmb_in.n(idxkeep);
        glmb_out.total_spawns = glmb_in.total_spawns(idxkeep) ; 
        glmb_out.total_births = glmb_in.total_births(idxkeep) ; 
        glmb_out.N_clutter = glmb_in.N_clutter(idxkeep) ;
        glmb_out.N_clutter_count= glmb_in.N_clutter_count(idxkeep);
        glmb_out.w= glmb_out.w/sum(glmb_out.w);
        for card=0:max(glmb_out.n)
            glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
        end
        for card = 0:max(glmb_out.total_spawns)
            glmb_out.cdn_spawn(card+1) = sum(glmb_out.w(glmb_out.total_spawns==card));
        end
    else
        glmb_out= glmb_in;
    end
end

function idx = find_vec(mat, vec)
    if ~isempty(mat)
        rows = size(mat,1) ; 
        vec = repmat(vec, [rows,1]) ; 
        temp = sum( (vec-mat) ~= 0, 2 ) ; 
        idx = find(temp==0) ; 
    else
        idx = [] ; 
    end
end