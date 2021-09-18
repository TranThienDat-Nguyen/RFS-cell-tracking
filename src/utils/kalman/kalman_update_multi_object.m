function [tt_sp1 , tt_sp2] = kalman_update_multi_object(model, z, z_norm_spawn, ah, tt_predict_sp1, tt_predict_sp2)
% kalman_update_multi_object performs joint update of tracks spawned from
% the same parent. It returns 2 independent track tables.
% z: joint kinematic measurements, z_norm_spawn: 2 separate modes measurements.
% Note: consider the input prediction is independent (since the prediction model
% induces no correlation terms anyway).
% Assuming single Gaussian likelihood

% grab the prediction info
w_1 = tt_predict_sp1.w ; 
m_1 = tt_predict_sp1.m ; 
P_1 = tt_predict_sp1.P ; 

w_2 = tt_predict_sp2.w ; 
m_2 = tt_predict_sp2.m ; 
P_2 = tt_predict_sp2.P ; 

w = (w_1.*w_2)/sum(w_1.*w_2) ; 

plength = length(w) ; 
P = zeros(2*model.x_dim,2*model.x_dim,plength) ; 
m = zeros(2*model.x_dim, plength) ; 
R = kron(eye(2),model.R); % observation noise
H = kron(eye(2),model.H); % observation matrix
logqz_all = zeros(1, plength) ; 

% perform joint update
for pidx = 1 : plength
    m(:,pidx) = [m_1(:,pidx) ; m_2(:,pidx)] ; 
    P(1:model.x_dim, 1:model.x_dim, pidx) = P_1(:,:,pidx) ; 
    P(model.x_dim+1:end, model.x_dim+1:end, pidx) = P_2(:,:,pidx) ; 
    
    mu = H * m(:,pidx) ; 
    S  =  R + H * P(:,:,pidx) * H';
    Vs = chol(S); 
    det_S = prod(diag(Vs))^2; 
    inv_sqrt_S = inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
    K  = P(:,:,pidx) * H' * iS;
    logqz = -0.5 * length(z) * log(2*pi) - 0.5 * log(det_S) - 0.5 *(z - mu)' * iS * (z - mu);
    logqz_all(pidx) = log(w(pidx)) + logqz ; 
    m(:,pidx) = m(:,pidx) + K*(z-mu) ; 
    P(:,:,pidx) = (eye(size(P(:,:,pidx))) - K * H ) * P(:,:,pidx);
end
logqz_out = logsumexp(logqz_all,[],2) ; 
w_out = exp(logqz_all - logqz_out) ;

% marginalize and form independent track tables
pD_1 = beta(tt_predict_sp1.st(1)+1,tt_predict_sp1.st(2))/ ...
    beta(tt_predict_sp1.st(1),tt_predict_sp1.st(2)) ; 

pD_2 = beta(tt_predict_sp2.st(1)+1,tt_predict_sp2.st(2))/ ...
    beta(tt_predict_sp2.st(1),tt_predict_sp2.st(2)) ; 

tt_sp1 = tt_predict_sp1 ;
tt_sp1.w = w_out ; 
tt_sp1.m = m(1:model.x_dim, :) ; 
tt_sp1.P = P(1:model.x_dim, 1:model.x_dim, :) ; 
tt_sp1.mode = tt_predict_sp1.mode .* z_norm_spawn(1,:) ; 
tt_sp1.logqz = logqz_out + log(sum(tt_sp1.mode)) ; 
tt_sp1.mode = tt_sp1.mode/sum(tt_sp1.mode) ; 
tt_sp1.st(1) = tt_sp1.st(1) + 1 ;

tt_sp1.ah = [tt_predict_sp1.ah, ah(1)] ; 
tt_sp1.pD = pD_1 ;

tt_sp2 = tt_predict_sp2 ;
tt_sp2.w = w_out ; 
tt_sp2.m = m(model.x_dim+1 : end, :) ; 
tt_sp2.P = P(model.x_dim+1 : end, model.x_dim+1 : end, :) ; 
tt_sp2.mode = tt_predict_sp2.mode .* z_norm_spawn(2,:) ; 
tt_sp2.logqz = log(sum(tt_sp2.mode)) ; % the joint kinematic qz is in the first tt
tt_sp2.mode = tt_sp2.mode/sum(tt_sp2.mode) ; 
tt_sp2.st(1) = tt_sp2.st(1) + 1 ;

tt_sp2.ah = [tt_predict_sp2.ah, ah(2)] ; 
tt_sp2.pD = pD_2 ;
    
    
    