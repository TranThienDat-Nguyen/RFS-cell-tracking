function [mixture_nextupdate , meas_ru] = exact_filter(model , filter , mixture , tt_birth , ~ , meas , k,clutter_table)
%% Create birth tracks
    numBirths = length(tt_birth) ; 
    r_b = zeros(numBirths, 1) ; 
    for bidx = 1 : numBirths
        r_b(bidx) = tt_birth{bidx}.r_b  ;
    end 
%% Generate children hypotheses ..
        % initilialization for parallel
    m= size(meas.Z{k},2); 
    meas_ru_ind = cell(length(mixture.w),1) ; 
    mixture_para_w = cell(length(mixture.w),1)  ; 
    mixture_para_I = cell(length(mixture.w),1);
    mixture_para_N_clutter =  cell(length(mixture.w),1) ; 
    mixture_para_N_clutter_count = cell(length(mixture.w),1) ;
    mixture_para_n = cell(length(mixture.w),1) ; 
    mixture_para_total_births = cell(length(mixture.w),1) ; 
    mixture_para_total_spawns = cell(length(mixture.w),1) ; 
    mixture_para_I_old = cell(length(mixture.w),1) ; 
    mixture_para_jtt = cell(length(mixture.w),1) ;             
                
    % loop through each prior hypothesis            
    parfor hidx = 1 : length(mixture.w) 
        if mixture.n(hidx) ~= 0 
            tt_in = jtt2tt(mixture.jtt{mixture.I(hidx)}, model) ; % convert joint track tables to independent track tables
            jtt_in = mixture.jtt{mixture.I(hidx)} ;     
        else
            tt_in = cell(0,1); 
        end
        sampling_factor = sqrt(mixture.w(hidx))/sum(sqrt(mixture.w)) ; % control number of components to be generated
        hypo = gen_hypo(model, tt_birth, r_b,  tt_in, meas, mixture.N_clutter(hidx), sampling_factor, filter, k) ;  % sample for components
        runidx = 1 ;
        jtt_idx = 1 ;
        jtt_new = cell(0,1) ; 
        % Run components update using sampled components
        for chidx = 1 : length(hypo)
            % processing birth tracks
            log_qz_birth = 0 ; 
            meas_ass_birth = hypo{chidx}.meas_ass(1:numBirths) ; 
            exist_birth = find(meas_ass_birth>=0) ;   
            this_birth_m = zeros(model.x_dim*length(exist_birth ),1) ; 
            this_birth_P = zeros(model.x_dim*length(exist_birth )) ;
            ibidx = 1 ; 
            for bidx = 1 : length(exist_birth)
                if meas_ass_birth(exist_birth(bidx)) == 0
                    this_birth_m(ibidx:ibidx+model.x_dim-1,1) = tt_birth{exist_birth(bidx)}.m ; 
                    this_birth_P(ibidx:ibidx+model.x_dim-1, ibidx:ibidx+model.x_dim-1) = tt_birth{exist_birth(bidx)}.P ; 
                elseif meas_ass_birth(exist_birth(bidx)) > 0
                    emm = meas_ass_birth(exist_birth(bidx)) ; 
                    [w_temp,m_temp,P_temp] = kalman_update_single_object(meas.Z{k}(:,emm),model, ...
                        tt_birth{exist_birth(bidx)}.w, tt_birth{exist_birth(bidx)}.m, tt_birth{exist_birth(bidx)}.P); 
                    log_qz_birth  = log_qz_birth  + log(w_temp) ;
                    this_birth_m(ibidx:ibidx+model.x_dim-1,1) = m_temp ; 
                    this_birth_P(ibidx:ibidx+model.x_dim-1, ibidx:ibidx+model.x_dim-1) = P_temp ; 
                end
                ibidx = ibidx + model.x_dim ; 
            end
            % processing survivals
            null_flag = true ; % indicate if this hypothesis has 0 cardinality
            log_qz_exist = 0 ; 
            if mixture.n(hidx)>0 % if there are tracks exist in previous time step
                curr_jtt = jtt_in ; 
                dead_indc = zeros(1, (length(hypo{chidx}.meas_ass)-numBirths)/3); % indices of dead tracks
                iidx = 1 ; 
                for ii = numBirths+1:3:length(hypo{chidx}.meas_ass)
                    if hypo{chidx}.meas_ass(ii)<0 && hypo{chidx}.meas_ass(ii+1)<0 
                        dead_indc(iidx) = iidx ; 
                    end
                    iidx = iidx + 1 ;
                end
                dead_indc(dead_indc==0) = [] ; 
                % remove the dead tracks from the joint track table
                ididx = 1 ;
                del_indc = zeros(model.x_dim*length(dead_indc),1) ; 
                for didx = 1 : length(dead_indc)
                    del_indc(ididx:ididx+model.x_dim-1) = (dead_indc(didx)-1)*model.x_dim+1 : dead_indc(didx)*model.x_dim ; 
                    ididx = ididx + model.x_dim ;
                end
                curr_jtt.m(del_indc,:) =[] ; 
                curr_jtt.P(del_indc,:,:) =[] ;
                if ~isempty(curr_jtt.P)
                    curr_jtt.P(:,del_indc,:) =[] ;
                end
            else
                curr_jtt.w = [] ;
                curr_jtt.m = [] ;
                curr_jtt.P = [] ;
            end
            % predict-update existing tracks
            if ~isempty(curr_jtt.m) % if there are still tracks exist
                [jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P] ...
                    = model.kalman_predict_mixture(model, curr_jtt.w, curr_jtt.m, curr_jtt.P, ...
                    hypo{chidx}.pred_model) ;
                 % clean up the Gaussian
                [jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P] ...
                    = gaus_merge(jtt_new{jtt_idx,1}.w', jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P, filter.merge_threshold) ;  
                [jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P] ...
                    = gaus_prune(jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P, filter.elim_threshold) ;  
                [jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P] ...
                    = gaus_cap(jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P, filter.L_max) ;  

                meas_ass_exist = hypo{chidx}.meas_ass(numBirths+1:end) ;
                if any(meas_ass_exist>0)
                    [log_qz_exist_temp, jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, jtt_new{jtt_idx,1}.P] ...
                        =kalman_update_mixture(model, jtt_new{jtt_idx,1}.w, jtt_new{jtt_idx,1}.m, ...
                        jtt_new{jtt_idx,1}.P, meas_ass_exist(meas_ass_exist>=0), meas.Z{k}) ;
                    log_qz_exist = log_qz_exist + log_qz_exist_temp ; 
                end
                % if birth exist then update them
                if ~isempty(this_birth_m) % if there are some new births
                    num_new_comp = length(jtt_new{jtt_idx,1}.w) ; 
                    mat_dim = length(this_birth_m) + size(jtt_new{jtt_idx,1}.m,1) ; 
                    m_temp = zeros(mat_dim  , num_new_comp) ; 
                    m_temp(1:length(this_birth_m),1) = this_birth_m ; 
                    m_temp(length(this_birth_m)+1:end,1) = jtt_new{jtt_idx,1}.m(:,1); 
                    P_temp = zeros(mat_dim, mat_dim, num_new_comp) ; 
                    P_temp(:,:,1) = blkdiag(this_birth_P, jtt_new{jtt_idx,1}.P(:,:,1)); 
                    if num_new_comp>1
                        for cidx = 2 : num_new_comp
                            m_temp(1:length(this_birth_m),cidx) = this_birth_m ; 
                            m_temp(length(this_birth_m)+1:end,cidx) = jtt_new{jtt_idx,1}.m(:,cidx); 
                            P_temp(:,:,cidx) = blkdiag(this_birth_P, jtt_new{jtt_idx,1}.P(:,:,cidx)) ; 

                        end
                    end
                    jtt_new{jtt_idx,1}.m = m_temp ; 
                    jtt_new{jtt_idx,1}.P = P_temp ; 
                end
                null_flag = false;
            else % if there is no exsiting tracks
                if ~isempty(this_birth_m) % if there are some new births
                    jtt_new{jtt_idx,1}.w = 1; 
                    jtt_new{jtt_idx,1}.m = this_birth_m ; 
                    jtt_new{jtt_idx,1}.P = this_birth_P ;
                    null_flag = false;
                end
            end
            % update the mixture 
            if ~null_flag % if there are tracks in the new joint track table, update them
                jtt_new{jtt_idx,1}.l = hypo{chidx}.l ; 
                jtt_new{jtt_idx,1}.st = hypo{chidx}.st ; 
                jtt_new{jtt_idx,1}.mode = hypo{chidx}.mode ;  
                jtt_new{jtt_idx,1}.pD = hypo{chidx}.pD ; 
                jtt_idx = jtt_idx + 1 ;
            end
            % compute the correct the hypothesis weight 
            log_w = log(mixture.w(hidx)) + hypo{chidx}.log_pE_weight + ...
                hypo{chidx}.log_qD_pD_qz_mode_weight + log_qz_birth + log_qz_exist ; 
            % compute clutter hypotheses
            meas_ru_ind_temp = false(1, m) ; 
            meas_ru_ind_temp( 1, hypo{chidx}.meas_ass(hypo{chidx}.meas_ass>0) ) = true ;
            if model.unknown_clutter
                num_Z_0 = m - sum(hypo{chidx}.meas_ass>0) ; 
                scidx = mixture.N_clutter(hidx) +1 ;
                N_top_this = sum(clutter_table{scidx , num_Z_0+1}(1,:)>0) ;  
                varClutterWeight = 1;
                N_top_this = ceil(min(N_top_this , filter.H_clutter) * varClutterWeight) ; 
                hyps_clutter = clutter_table{scidx , num_Z_0+1}(: , 1 : N_top_this) ;
                for clidx = 1 : size(hyps_clutter , 2)
                    mixture_para_w{hidx}(runidx) = log_w + hyps_clutter(2,clidx) ; 
                    if ~null_flag
                        mixture_para_I{hidx}(runidx) = jtt_idx-1 ;
                    else
                        mixture_para_I{hidx}(runidx) = 0;
                    end
                    mixture_para_N_clutter{hidx}(runidx) =  hyps_clutter(1,clidx)-1 ; 
                    mixture_para_N_clutter_count{hidx}(runidx) = num_Z_0 ;
                    mixture_para_n{hidx}(runidx) = hypo{chidx}.n ; 
                    mixture_para_total_births{hidx}(runidx) = hypo{chidx}.nbirths ; 
                    mixture_para_total_spawns{hidx}(runidx) = hypo{chidx}.nspawns ; 
                    mixture_para_I_old{hidx}(runidx) = hidx ;  
                    meas_ru_ind{hidx}(runidx,:) = meas_ru_ind_temp ; 
                    runidx = runidx + 1;
                end
           else
                mixture_para_w{hidx}(runidx) = log_w - model.lambda_c + num_Z_0*log(model.lambda_c*model.pdf_c) ; 
                if ~null_flag
                    mixture_para_I{hidx}(runidx) = jtt_idx-1 ;
                else
                    mixture_para_I{hidx}(runidx) = 0;
                end
                mixture_para_N_clutter{hidx}(runidx) =  model.lambda_c ; 
                mixture_para_N_clutter_count{hidx}(runidx) = num_Z_0 ;
                mixture_para_n{hidx}(runidx) = hypo{chidx}.n ; 
                mixture_para_total_births{hidx}(runidx) = hypo{chidx}.nbirths ; 
                mixture_para_total_spawns{hidx}(runidx) = hypo{chidx}.nspawns ; 
                mixture_para_I_old{hidx}(runidx) = hidx ;  
                meas_ru_ind{hidx}(runidx,:) = meas_ru_ind_temp ; 
                runidx = runidx + 1;
            end
        end % children hypotheses
        mixture_para_jtt{hidx} = jtt_new ;
        if jtt_idx == 1 && null_flag
            mixture_para_jtt{hidx} = cell(0,1); % if there is no tracks, set the joint track table to empty
        else
            mixture_para_jtt{hidx} = jtt_new ; % if there are tracks, get the new joint track table
        end
    end 
    % Bring the para structures to standard format
    runidx = 1 ;
    curr_jtt_idx = 0 ; 
    mixture_nextupdate.jtt = cell(0,1) ; 
    meas_ru_ind_temp = meas_ru_ind ; 
    meas_ru_ind = zeros(0,m) ; 
    for hidx = 1 : length(mixture.w)
        for chidx = 1 : length(mixture_para_w{hidx})
            mixture_nextupdate.w(runidx) = mixture_para_w{hidx}(chidx) ;
            if mixture_para_I{hidx}(chidx) == 0 
                mixture_nextupdate.I(runidx) = 0 ; 
            else
                mixture_nextupdate.I(runidx) = mixture_para_I{hidx}(chidx) + curr_jtt_idx ; 
            end
            mixture_nextupdate.N_clutter(runidx) =  mixture_para_N_clutter{hidx}(chidx) ; 
            mixture_nextupdate.N_clutter_count(runidx) = mixture_para_N_clutter_count{hidx}(chidx) ;
            mixture_nextupdate.n(runidx) = mixture_para_n{hidx}(chidx) ; 
            mixture_nextupdate.total_births(runidx) = mixture_para_total_births{hidx}(chidx) ; 
            mixture_nextupdate.total_spawns(runidx) = mixture_para_total_spawns{hidx}(chidx) ; 
            mixture_nextupdate.I_old(runidx) = hidx ;
            runidx = runidx + 1 ;
        end
        mixture_nextupdate.jtt = [mixture_nextupdate.jtt; mixture_para_jtt{hidx}] ; 
        meas_ru_ind = [meas_ru_ind; meas_ru_ind_temp{hidx}] ; 
        curr_jtt_idx = curr_jtt_idx + length(mixture_para_jtt{hidx}) ; 
    end        
    % Normalize weights and compute cdn distribution
    mixture_nextupdate.w= exp(mixture_nextupdate.w-logsumexp(mixture_nextupdate.w,[],2));   
    meas_ru = meas_ru_ind' * mixture_nextupdate.w' ;
    for card=0:max(mixture_nextupdate.n)
        mixture_nextupdate.cdn(card+1)= sum(mixture_nextupdate.w(mixture_nextupdate.n==card));
    end
    for card = 0:max(mixture_nextupdate.total_spawns)
        mixture_nextupdate.cdn_spawn(card+1) = sum(mixture_nextupdate.w(mixture_nextupdate.total_spawns==card));
    end
    % clean the mixture density
    mixture_nextupdate = prune_mixture(mixture_nextupdate, filter);  
    mixture_nextupdate = cap_mixture(mixture_nextupdate, filter);  
    mixture_nextupdate = clean_mixture(mixture_nextupdate);  
    %% Estimate tracks for each hypothesis
    for hidx = 1 : length(mixture_nextupdate.w)
        est_X_pD = zeros(model.x_dim+1, mixture_nextupdate.n(hidx)) ; 
        est_L = cell(1, mixture_nextupdate.n(hidx)) ; 
        mixture_nextupdate.est_X_pD{hidx} = mixture.est_X_pD{mixture_nextupdate.I_old(hidx)} ; 
        mixture_nextupdate.est_L{hidx} = mixture.est_L{mixture_nextupdate.I_old(hidx)} ; 
        mixture_nextupdate.est_others{hidx} = mixture.est_others{mixture_nextupdate.I_old(hidx)} ;
        if mixture_nextupdate.n(hidx)>0
            % check for the best component
            [~,idxtrk]= max(mixture_nextupdate.jtt{mixture_nextupdate.I(hidx)}.w);
            ii = 1 ; 
            for n = 1 : mixture_nextupdate.n(hidx)
                est_X_pD(:, n) = [mixture_nextupdate.jtt{mixture_nextupdate.I(hidx)}.m(ii:ii+model.x_dim-1, idxtrk); ...
                    mixture_nextupdate.jtt{mixture_nextupdate.I(hidx)}.pD(n)] ; 
                ii = ii + model.x_dim ; 
            end
            est_L = mixture_nextupdate.jtt{mixture_nextupdate.I(hidx)}.l ; 
        end
        mixture_nextupdate.est_X_pD{hidx}{k} = est_X_pD ;
        mixture_nextupdate.est_L{hidx}{k} = est_L ;
        mixture_nextupdate.est_others{hidx}(:,k) = [mixture_nextupdate.n(hidx) ; 
                                                 mixture_nextupdate.total_births(hidx) ; 
                                                 mixture_nextupdate.total_spawns(hidx) ;
                                                 mixture_nextupdate.N_clutter(hidx) ;
                                                 mixture_nextupdate.N_clutter_count(hidx)]  ;
    end
end

    
    
    
function mixture_clean= clean_mixture(mixture_temp)
%flag used tracks
usedindicator= zeros(length(mixture_temp.jtt),1);
for hidx= 1:length(mixture_temp.w)
    if mixture_temp.I(hidx)~=0
        usedindicator(mixture_temp.I(hidx))= usedindicator(mixture_temp.I(hidx))+1;
    end
end
trackcount= sum(usedindicator>0);

%remove unused tracks and reindex existing hypotheses/components
newindices= zeros(length(mixture_temp.jtt),1); newindices(usedindicator>0)= 1:trackcount;
mixture_clean.jtt= mixture_temp.jtt(usedindicator>0);
mixture_clean.w= mixture_temp.w;
for hidx= 1:length(mixture_temp.w)
    if mixture_temp.I(hidx)~=0
        mixture_clean.I(hidx)= newindices(mixture_temp.I(hidx));
    else
        mixture_clean.I(hidx) = 0 ; 
    end
end
mixture_clean.n= mixture_temp.n;
mixture_clean.cdn= mixture_temp.cdn;
mixture_clean.cdn_spawn= mixture_temp.cdn_spawn;
mixture_clean.total_births= mixture_temp.total_births;
mixture_clean.total_spawns= mixture_temp.total_spawns;
mixture_clean.N_clutter= mixture_temp.N_clutter;
mixture_clean.N_clutter_count= mixture_temp.N_clutter_count;
mixture_clean.I_old = mixture_temp.I_old ; 
end


function mixture_out= prune_mixture(mixture_in,filter)
%prune components with weights lower than specified threshold
idxkeep= find(mixture_in.w > filter.hyp_threshold);
if length(idxkeep)<filter.H_min
    [~,idxkeep] = sort(mixture_in.w, 'descend') ; 
    idxkeep = idxkeep(1:filter.H_min) ; 
end
mixture_out.jtt= mixture_in.jtt;
mixture_out.w= mixture_in.w(idxkeep);
mixture_out.I= mixture_in.I(idxkeep);
mixture_out.n= mixture_in.n(idxkeep);
mixture_out.total_spawns = mixture_in.total_spawns(idxkeep) ; 
mixture_out.total_births = mixture_in.total_births(idxkeep) ; 
mixture_out.N_clutter = mixture_in.N_clutter(idxkeep) ;
mixture_out.N_clutter_count= mixture_in.N_clutter_count(idxkeep);
mixture_out.w= mixture_out.w/sum(mixture_out.w);
for card=0:max(mixture_out.n)
    mixture_out.cdn(card+1)= sum(mixture_out.w(mixture_out.n==card));
end
for card = 0:max(mixture_out.total_spawns)
    mixture_out.cdn_spawn(card+1) = sum(mixture_out.w(mixture_out.total_spawns==card));
end
mixture_out.I_old = mixture_in.I_old(idxkeep) ; 
end



function mixture_out= cap_mixture(mixture_in,filter)
%cap total number of components to specified maximum
if length(mixture_in.w) > filter.H_max
    [~,idxsort]= sort(mixture_in.w,'descend');
    idxkeep=idxsort(1:filter.H_max);
    mixture_out.jtt= mixture_in.jtt;
    mixture_out.w= mixture_in.w(idxkeep);
    mixture_out.I= mixture_in.I(idxkeep);
    mixture_out.I_old= mixture_in.I_old(idxkeep);
    mixture_out.n= mixture_in.n(idxkeep);
    mixture_out.n= mixture_in.n(idxkeep);
    mixture_out.total_spawns = mixture_in.total_spawns(idxkeep) ; 
    mixture_out.total_births = mixture_in.total_births(idxkeep) ; 
    mixture_out.N_clutter = mixture_in.N_clutter(idxkeep) ;
    mixture_out.N_clutter_count= mixture_in.N_clutter_count(idxkeep);
    mixture_out.w= mixture_out.w/sum(mixture_out.w);
    for card=0:max(mixture_out.n)
        mixture_out.cdn(card+1)= sum(mixture_out.w(mixture_out.n==card));
    end
    for card = 0:max(mixture_out.total_spawns)
        mixture_out.cdn_spawn(card+1) = sum(mixture_out.w(mixture_out.total_spawns==card));
    end
else
    mixture_out= mixture_in;
end
end

        
    