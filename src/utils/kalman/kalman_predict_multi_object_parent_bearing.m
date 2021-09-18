function [w_pred , m_out , P_out] = kalman_predict_multi_object_parent_bearing(model , w_kin_in , m_kin_in , P_kin_in , spawn_idx , marg_flag)
% kalman_predict_multi_object performs Kalman prediction from multiple
% input objects to multiple output objects.
% marg_flag: marginalizing the joint objects density to form independent single-object density.
% Note: current usage only considers 1 input object.

    if ~isempty(spawn_idx)
        [w_trans , F , Q] = gen_F_Q (model , 1);
        [w_pred , m_pred , P_pred] = kalman_predict(model,w_trans,F,Q,w_kin_in,m_kin_in,P_kin_in,spawn_idx)  ;               
    else
        [w_trans , F , Q] = gen_F_Q (model , 2);
        [w_pred , m_pred , P_pred] = kalman_predict(model,w_trans,F,Q,w_kin_in,m_kin_in,P_kin_in,spawn_idx)  ;  
    end
    numobj_out = size(m_pred ,1) / model.x_dim ;
    
    % marginalizing the tracks
    if marg_flag && numobj_out > 1
        m_out = cell(numobj_out , 1) ; 
        P_out = cell(numobj_out , 1) ;   
        for oidx = 1 : numobj_out
            span = model.x_dim*(oidx-1) + 1 : model.x_dim * oidx ; 
            m_out{oidx} = m_pred(span,:) ; 
            P_out{oidx} = P_pred(span , span , :) ; 
        end
    else 
        m_out = m_pred ; 
        P_out = P_pred ;
    end
end

function [ w_kin_pred , m_kin_pred , P_kin_pred ] = kalman_predict(model,w_trans,F,Q, w_kin_in , m_kin_in , P_kin_in , spawn_idx)
    numObjout = 1 + length(spawn_idx) ;
    numCompout = length(w_trans) * length(w_kin_in) ; 
    if ~isempty(spawn_idx)
        numCompout = numCompout * length(model.wsp) ;
    end

    w_kin_pred = zeros(1 , numCompout) ;
    m_kin_pred = zeros(numObjout * model.x_dim , numCompout) ; 
    P_kin_pred = zeros(numObjout * model.x_dim , numObjout * model.x_dim , numCompout) ; 
    ocidx = 1 ; 
    
    % loop through previous components
    for cidx = 1 : length(w_kin_in)
        % check if current object spawns
        if ~isempty(spawn_idx) % if yes, stacking matrices on top of each other
            d = zeros (numObjout*model.x_dim , 1) ;
            w_sets = 1 ; 
            d_s = model.spawn_dist ; 
            delta = rad2deg(atan2(m_kin_in(4,cidx),m_kin_in(3,cidx)))+90 ; 
            dx = d_s * cosd(delta) ; 
            dy = d_s * sind(delta) ; 
            d(model.xy_pos) = [dx ; dy] ; 
            d(model.xy_pos+model.x_dim) = [-dx ; -dy] ; 
        else
            d = zeros (numObjout*model.x_dim , 1) ;
            w_sets = 1 ; 
        end
        
        % loop through all transitions
        for fidx = 1 : size(F,3)
            for sgidx = 1 : size(d,2)
                m_kin_pred(:,ocidx) = F(:,:,fidx) * m_kin_in(:,cidx) + d(:,sgidx);
                P_kin_pred(:,:,ocidx) = Q(:,:,fidx) + F(:,:,fidx) * P_kin_in(:,:,cidx) * F(:,:,fidx)' ; 
                P_kin_pred(:,:,ocidx) = (P_kin_pred(:,:,ocidx) + P_kin_pred(:,:,ocidx)') / 2 ; % just to avoid numerical issue
                w_kin_pred(ocidx) = w_kin_in(cidx) * w_trans(fidx) *w_sets (sgidx); 
                ocidx = ocidx + 1 ; 
            end
        end
    end
end

function [w , F , Q] = gen_F_Q (model , type) 
    if type == 1 % single object with spawn
        F = zeros(2*model.x_dim , model.x_dim , 1) ; % double the rows as stack on top 
        Q = zeros(2*model.x_dim , 2*model.x_dim , 1) ; % covariance for  2 objects
        w = model.wtsp ; 
        for fidx = 1 : length(model.wtsp)
            F(:,:,fidx) = [ model.S(:,:,fidx) ; model.S(:,:,fidx) ] ; 
            Q(:,:,fidx) = blkdiag(model.Q_s , model.Q_s) ;
        end
    elseif type == 2  %single object without spawn
        F = zeros(model.x_dim , model.x_dim , length(model.ws)) ; 
        Q = F ; % same size as F 
        w = model.ws ; 
        for fidx = 1 : length(model.ws)
            F(:,:,fidx) = model.F(:,:,fidx) ;
            Q(:,:,fidx) = model.Q(:,:,fidx) ; 
        end
    end          
end
            
            
             

        