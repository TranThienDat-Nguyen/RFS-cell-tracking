function [qz_update,m_update,P_update] = kalman_update_single_object(z,model,w,m,P)
% kalman_update_single_object performs single object Kalman update

    plength= size(m,2);
    zlength= size(z,2);
    qz_update= zeros(zlength , plength);
    m_update = zeros(model.x_dim,plength,zlength);
    P_update = zeros(model.x_dim,model.x_dim,plength);
    H = model.H ; 
    R = model.R ; 
    for idxp=1:plength
        eta = H * m(:,idxp) ;
        S= R + H * P(:,:,idxp) * H' ; S= (S+ S')/2;   % addition step to avoid numerical problem
        Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
        K  = P(:,:,idxp)*H'*iS;
        qz_temp = exp(-0.5*size(z,1)*log(2*pi) - 0.5*log(det_S) - 0.5*dot(z-repmat(eta,[1 size(z,2)]),iS*(z-repmat(eta,[1 size(z,2)]))))';
        m_temp = repmat(m(:,idxp),[1 size(z,2)]) + K*(z-repmat(eta,[1 size(z,2)]));
        P_temp = (eye(size(P(:,:,idxp)))-K*H)*P(:,:,idxp);
        qz_update(:,idxp)   = w(idxp) * qz_temp + eps;
        m_update(:,idxp,:) = m_temp;
        P_update(:,:,idxp) = P_temp;
    end
end




 



