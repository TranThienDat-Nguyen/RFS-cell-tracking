function tt_birth = gen_meas_driven_birth(model , meas , meas_ru , k , bindc, offset)
% gen_meas_driven_birth generate adaptive births from the measurement
    if nargin == 5
        scale_flag =  false ;
    else
        scale_flag = true ;
        max_scale = 1 ; 
        min_scale = 0.5 ; 
    end
    Z = meas.Z{k} ; 
    meas_rb = (1 - meas_ru) / sum(1-meas_ru) ;
    meas_rb = meas_rb + eps ; % do this to avoid zero-weight
    tt_birth = cell(size(Z , 2) , 1) ;
    for bidx = 1 : length(tt_birth)
        m_temp = zeros(model.x_dim ,1);
        m_temp(model.xy_pos,:) = Z(:,bidx);
        tt_birth{bidx}.m = m_temp ; 
        tt_birth{bidx}.P = model.Q_b ; 
        tt_birth{bidx}.w = 1 ; 
        if k == 1
            tt_birth{bidx}.r_b = model.pB_init ; 
        else
            if scale_flag
                scale = max_scale - (max_scale - min_scale) * ...
                    gbellmf(tt_birth{bidx}.m(model.xy_pos(1)),[model.range_c(1,2)/2-offset 40 model.range_c(1,2)/2]);
                scale = scale * (max_scale - (max_scale - min_scale) * ...
                    gbellmf(tt_birth{bidx}.m(model.xy_pos(2)),[model.range_c(2,2)/2-offset 40 model.range_c(2,2)/2]));
            else
                scale =1 ;
            end
            tt_birth{bidx}.r_b = scale * min([model.pB_max , model.lambda_b*meas_rb(bidx)]); 
        end
        tt_birth{bidx}.ah = [] ;
        tt_birth{bidx}.l= [k; bindc(bidx)] ;
        tt_birth{bidx}.mode = [(1-model.pSp) , model.pSp] ;
        tt_birth{bidx}.st = model.init_st ; 
    end
end
    