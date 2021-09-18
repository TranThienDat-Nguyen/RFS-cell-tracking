function [w_out , m_out , P_out] = kalman_predict_mixture_parent_bearing(model , w_kin_in , m_kin_in , P_kin_in , pred_model)
% kalman_predict_mixture performs kalman prediction for joint object density
% Note: the model to update is given by pred_model, which is
%       length of pred_models is number of models
%       first row of each pred_models cell: transition model (negative denotes spawning object)
%       second row of each pred_models cell: shifting model

    % prepare the models
    num_models = length(pred_model) ; 
    num_objects = size(m_kin_in,1)/model.x_dim ; 
    all_F = cell(num_models,1) ; 
    all_Q = cell(num_models,1) ; 
    w_model = ones(num_models,1) ; 
    all_d = cell(num_models,1) ; 
    for midx = 1 : num_models
        temp_F = cell(num_objects,1) ; 
        temp_Q = cell(num_objects,1) ; 
        num_surv = sum(pred_model{midx}(1,:)>0) ; 
        num_spawn = sum(pred_model{midx}(1,:)<0) ; 
        temp_d = zeros(model.x_dim*num_surv + 2*model.x_dim*num_spawn ,1) ; 
        ii = 1 ; 
        for n = 1 : num_objects
            if pred_model{midx}(1,n) > 0 % this is a survival
                temp_F{n} = model.F(:,:,pred_model{midx}(1,n)) ; 
                temp_Q{n} = model.Q(:,:,pred_model{midx}(1,n)) ;
                w_model(midx) = w_model(midx)*model.ws(pred_model{midx}(1,n)) ; 
                ii = ii + model.x_dim ; 
            elseif pred_model{midx}(1,n)<0 % this is a spawning object
                % stack matrices on top of each other
                temp_F{n} = repmat(model.S(:,:,-pred_model{midx}(1,n)),2,1) ; 
                temp_Q{n} = blkdiag(model.Q_s(:,:,-pred_model{midx}(1,n)),...
                    model.Q_s(:,:,-pred_model{midx}(1,n))) ;
                temp = zeros(model.x_dim,1) ; 
                parent_state = m_kin_in((n-1)*model.x_dim+1 : n*model.x_dim) ; 
                d_s = model.spawn_dist ; 
                delta = rad2deg(atan2(parent_state(4),parent_state(3)))+90 ; 
                dx = d_s * cosd(delta) ; 
                dy = d_s * sind(delta) ; 
                temp(model.xy_pos,1) = [dx ; dy] ;
                temp_d(ii:ii+model.x_dim-1,:) = temp; 
                ii = ii + model.x_dim ; 
                temp_d(ii:ii+model.x_dim-1,:) = -temp ; 
                ii = ii + model.x_dim ; 
                w_model(midx) = w_model(midx) * model.wsp(pred_model{midx}(2,n)) ;      
            end
        end
        all_F{midx} = blkdiag(temp_F{:}) ;
        all_Q{midx} = blkdiag(temp_Q{:}) ;
        all_d{midx} = temp_d ; 
    end
    w_model = w_model/sum(w_model) ; 
    mat_dim = size(all_F{1},1) ; 
    num_init = length(w_kin_in) ; 
    w_out = zeros(num_models*num_init, 1) ; 
    m_out = zeros(mat_dim, num_models*num_init) ; 
    P_out = zeros(mat_dim, mat_dim, num_models*num_init) ; 
    cidx = 1 ; 
    
    % perform the prediction
    for pidx = 1 : num_init
        for midx =1 : num_models
            F = all_F{midx}; 
            Q = all_Q{midx}; 
            w_out(cidx) = w_kin_in(pidx) * w_model(midx) ; 
            m_out(:,cidx) = F*m_kin_in(:,pidx) + all_d{midx} ; 
            P_out(:,:,cidx) = F*P_kin_in(:,:,pidx)*F' +Q; 
            P_out(:,:,cidx)  = (P_out(:,:,cidx) + P_out(:,:,cidx)') / 2 ;
            cidx = cidx + 1 ;
        end
    end
end

            
             

        