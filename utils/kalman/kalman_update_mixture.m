function [logqz_all, w_out, m_out, P_out] = kalman_update_mixture(model, w, m, P, meas_ass, Z)
% kalman_update_mixture performs update the joint objects density

    pos_indc = find(meas_ass>0) ; 
    num_pos_indc = length(pos_indc) ;
    if num_pos_indc>0
        % select portion of matrix corresponding to detected objects
        rem_indc = find(meas_ass==0) ;
        jj = 1 ;
        del_indc = zeros(length(rem_indc)*model.x_dim,1) ; 
        for ii = 1 : length(rem_indc)
            del_indc(jj:jj+model.x_dim-1,1) = (rem_indc(ii)-1)*model.x_dim+1 : rem_indc(ii)*model.x_dim ; 
            jj = jj + model.x_dim ; 
        end
        m_in = m ; 
        m(del_indc,:) = [] ; 
        P_in = P ; 
        P(del_indc,:,:) = [] ; 
        if ~isempty(P)
            P(:,del_indc,:) = [] ; 
        end

        % run Kalman update on the selected portion
        rep_indc = zeros(num_pos_indc*model.x_dim,1) ; 
        jj = 1 ;
        for ii = 1 :num_pos_indc
            rep_indc(jj:jj+model.x_dim-1,1) = (pos_indc(ii)-1)*model.x_dim+1 : pos_indc(ii)*model.x_dim ; 
            jj = jj + model.x_dim ; 
        end
        for cidx = 1 : size(P,3)
            P(:,:,cidx) = (P(:,:,cidx) + P(:,:,cidx)')/2; 
            [~,flag] = chol(P(:,:,cidx)) ; 
            if flag~=0
                P(:,:,cidx) = nearestSPD(P(:,:,cidx));
            end
        end
            % fabricate the H, R matrices
        R = kron(eye(num_pos_indc),model.R); % observation noise
        H = kron(eye(num_pos_indc),model.H); % observation matrix
            % fabricate the joint measurements
        meas_indc = meas_ass(meas_ass>0) ; 
        z = zeros(model.z_dim * length(meas_indc), 1) ; 
        ii = 1 ; 
        for midx =1 : length(meas_indc)
            z(ii : ii + model.z_dim - 1,1) = Z(:,meas_indc(midx)) ; 
            ii = ii + model.z_dim ; 
        end
        w_out = w ; 
        for pidx = 1 : length(w)
            mu = H * m(:,pidx) ; 
            S  =  R + H * P(:,:,pidx) * H';
            Vs = chol(S); 
            det_S = prod(diag(Vs))^2; 
            inv_sqrt_S = inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
            K  = P(:,:,pidx) * H' * iS;
            logqz = -0.5 * length(z) * log(2*pi) - 0.5 * log(det_S) - 0.5 *(z - mu)' * iS * (z - mu);
            w_out(pidx) = w_out(pidx) * exp(logqz) ; 
            m(:,pidx) = m(:,pidx) + K*(z-mu) ; 
            P(:,:,pidx) = (eye(size(P(:,:,pidx))) - K * H ) * P(:,:,pidx);
        end

        % put the measurement-updated portion of the matrix back to the original matrix
        logqz_all = log(sum(w_out)) ; 
        w_out = w_out/sum(w_out) ; 
        m_out = m_in ; 
        m_out(rep_indc,:) = m ; 
        P_out = P_in ; 
        P_out(rep_indc,rep_indc,:) = P ; 
    else % if all objects are miss-detected, no need to run joint update
        logqz_all = 0 ; 
        w_out = w ; 
        m_out = m ; 
        P_out = P ;
    end

end
    
    



    