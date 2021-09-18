function log_nckval = lognck( n, c )
%nchoosek function in log domain
%
    log_nckval = 0;
    fac_n = 0;
    fac_c = 0;
    fac_n_c = 0;
    
    for xx = 1:n
        fac_n = fac_n + log(xx);
        if ( xx <= c)
            fac_c = fac_c + log(xx);
        end
        if ( xx <= (n-c))
            fac_n_c = fac_n_c + log(xx);
        end        
    end
    log_nckval = fac_n - fac_c - fac_n_c;
    %disp(exp(log_nckval));
end



