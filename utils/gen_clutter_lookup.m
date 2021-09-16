function [optm_costm_clutter_first , optm_costm_clutter] = gen_clutter_lookup(model, N_max_clt_births_first, N_max_cl, N_top_cl, N_max_z, N_min_z)
% I grab this from Yuthika codes
% Unknown clutter related parameters ----- --------------------------------

%clutter target birth modeal
N_max_clt_births_first = 100;

clutter_bar_q_first = model.clutter_P_B;
log_cl_bar_q_first = log(clutter_bar_q_first);
log_cl_neg_bar_q_first = log(1 - clutter_bar_q_first);

N_max_clt_births = model.clutter_gen;%

clutter_bar_q = model.clutter_P_B;

clutter_Ps = model.clutter_P_S;
clutter_Qs = 1 - clutter_Ps;
clutter_Pd = model.clutter_P_D;
clutter_Qd = 1 - clutter_Pd;



log_cl_bar_q = log(clutter_bar_q);
log_cl_neg_bar_q = log(1 - clutter_bar_q);
log_cl_Ps = log(clutter_Ps);
log_cl_Qs = log(clutter_Qs);
clutter_Pd = clutter_Pd * model.pdf_c;
log_cl_Pd = log(clutter_Pd);
log_cl_Qd = log(clutter_Qd);



% clval_first = clutterpdf;
% clval = clutterpdf;

N_max_cl = 100; 
N_top_cl = 100;
N_max_z = 100;
N_min_z = 0;

%Unknown clutter related pre-calculations --------------------------------
%time step 1
costm_clutter_first = cell( 1, N_max_z-N_min_z+1); % no existing clutter targets
optm_costm_clutter_first = cell( 1, N_max_z-N_min_z+1);
zidx = 0 ; 
for z_a = N_min_z:N_max_z
    zidx = zidx + 1; % number of unassigned measurements
    for b_a = 0:N_max_clt_births_first
        bidx = b_a + 1; % number of births
        
        if (b_a >= z_a )
            costm_clutter_first{zidx}(bidx,1) = lognck(N_max_clt_births_first, b_a) + (log_cl_bar_q_first * b_a) + (log_cl_neg_bar_q_first * (N_max_clt_births_first - b_a)) + ...
                                                sum(log(1:z_a)) +  lognck( b_a, z_a) + (log_cl_Pd * z_a) + (log_cl_Qd * ( b_a - z_a)); 
                                                 
        else
            costm_clutter_first{zidx}(bidx,1) = -inf;            
        end
    end
    [sorted_costm_clutter_total, idxs] = sort( -costm_clutter_first{zidx});    
    optm_costm_clutter_first{zidx}(1,1:N_top_cl) = idxs(1:N_top_cl);    
    optm_costm_clutter_first{zidx}(2,1:N_top_cl) = -sorted_costm_clutter_total(1:N_top_cl);    
end
clear sorted_costm_clutter_total; clear idxs;
%for after 1st time step
%initialization
costm_clutter = cell( N_max_cl+1, N_max_z-N_min_z+1);
costm_clutter_total = cell( N_max_cl+1, N_max_z-N_min_z+1);
optm_costm_clutter = cell( N_max_cl+1, N_max_z-N_min_z+1);
for cidx = 1:N_max_cl+1
    for zidx = 1 : N_max_z-N_min_z+1
        costm_clutter{cidx,zidx} = zeros( N_max_clt_births+1, cidx);
        costm_clutter_total{cidx,zidx} = zeros( 1, N_max_cl+1);
        optm_costm_clutter{cidx,zidx} = zeros( 1, min(cidx,N_top_cl));
    end
end
%calculations
for cidx = 1:N_max_cl+1
    c_a = cidx - 1;
    zidx = 0 ; 
    for z_a = N_min_z:N_max_z
        zidx = zidx + 1;        
        for b_a = 0:N_max_clt_births
            bidx = b_a + 1;
            for s_a = 0:c_a
                sidx = s_a + 1;
                t_a = b_a + s_a;
                if ( t_a >= z_a )
                    costm_clutter{cidx,zidx}(bidx,sidx) = lognck(N_max_clt_births, b_a) + (log_cl_bar_q * b_a) + (log_cl_neg_bar_q * (N_max_clt_births - b_a)) + ...
                                        lognck(c_a, s_a) + (log_cl_Ps * s_a) + (log_cl_Qs * (c_a - s_a)) + ...
                                        sum(log(1:z_a)) +  lognck( t_a, z_a) + (log_cl_Pd * z_a) + (log_cl_Qd * (t_a - z_a)); 
                else
                    costm_clutter{cidx,zidx}(bidx,sidx) = -inf;
                end                
            end
        end
      
        costm_clutter_total{cidx,zidx}(1,1:N_max_cl+1) = -inf;
        for b_a = 0:N_max_clt_births
            for s_a = 0:c_a
                tidx = b_a+s_a+1;
                if ( tidx <= N_max_cl+1)                    
                    costm_clutter_total{cidx,zidx}(1,tidx) = logsumexp([costm_clutter_total{cidx,zidx}(1,tidx), ...
                                                                costm_clutter{cidx,zidx}(b_a+1,s_a+1)],[],2);
                end
            end
        end              
        [sorted_costm_clutter_total, idxs] = sort( -costm_clutter_total{cidx,zidx});
%         optm_costm_clutter{cidx,zidx}(1,1:N_top_cl) = idxs(1:N_top_cl);
%         optm_costm_clutter{cidx,zidx}(2,1:N_top_cl) = -sorted_costm_clutter_total(1:N_top_cl);
        optm_costm_clutter{cidx,zidx}(1,1:N_top_cl) = idxs(1:N_top_cl);
        optm_costm_clutter{cidx,zidx}(2,1:N_top_cl) = -sorted_costm_clutter_total(1:N_top_cl);
    end
end