function est = extract_tracks(glmb, model)
    K = length(glmb.est_L{1}) ; 
    est(1).X= cell(K,1);
    est(1).pD= cell(K,1) ; 
    [~,mode] = max(glmb.cdn);
    N = mode-1;
    [~,idxcmp]= max(glmb.w.*(glmb.n==N));
    est(1).L = glmb.est_L{idxcmp} ; 
    est(1).N = glmb.est_others{idxcmp}(1,:) ; 
    est(1).nbirths = glmb.est_others{idxcmp}(2,:) ; 
    est(1).nspawns = glmb.est_others{idxcmp}(3,:) ; 
    est(1).nclutter = glmb.est_others{idxcmp}(4,:) ; 
    est(1).nclutter_count = glmb.est_others{idxcmp}(5,:) ; 
    
    for k = 1 : K
        if ~isempty(glmb.est_X_pD{idxcmp}{k})
            est(1).X{k} = glmb.est_X_pD{idxcmp}{k}(1:model.x_dim,:) ;
            est(1).pD{k} = glmb.est_X_pD{idxcmp}{k}(end,:) ;
        end
    end
end
