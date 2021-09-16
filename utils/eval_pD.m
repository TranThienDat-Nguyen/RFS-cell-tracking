function pD = eval_pD(glmb)
pD = 0 ; 
for hidx = 1 : length(glmb.I)
    if ~isempty(glmb.I{hidx})
        temp_pD = 0 ; 
        for tidx = 1 : length(glmb.I{hidx})
            temp_pD = temp_pD + glmb.tt{glmb.I{hidx}(tidx)}.pD ; 
        end
        temp_pD = temp_pD/length(glmb.I{hidx}) ; 
        pD = pD + temp_pD*glmb.w(hidx) ; 
    else
        pD = pD + glmb.w(hidx) ; 
    end
end