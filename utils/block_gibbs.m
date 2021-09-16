function [assignments]= block_gibbs(P0,m,nbirths)
n1 = size(P0,1);
burn_in = 0 ; 
assignments= zeros(m,n1);

% Use no births + miss detection of survivals as initial solution
currsoln = zeros(1,n1) ; 
currsoln(1:nbirths) = (1:nbirths) ; % not exist
for ii = nbirths+1 : n1
    if mod(ii-nbirths,3) == 1 % miss detection
        currsoln(ii) = ii + n1 ; 
    else
        currsoln(ii) = ii ; % not existing
    end
end
% main sampling
assignments(1,:)= currsoln;
for sol= 2:m+burn_in
    var = 1 ;
    while var<=n1
        notBirth = var-nbirths>0 ; 
        if ~notBirth % do the normal Gibbs sampling
            tempsamp= exp(-P0(var,:)); %grab row of costs for current association variable
            tempsamp(currsoln([1:var-1,var+1:end]))= 0; %lock out current and previous iteration step assignments except for the one in question
            idxold= find(tempsamp>0); tempsamp= tempsamp(idxold);
            [~,currsoln(var)]= histc(rand(1,1),[0;cumsum(tempsamp(:))/sum(tempsamp)]); 
            currsoln(var)= idxold(currsoln(var));
            var = var + 1 ; 
        else % do the join Gibbs sampling
             % need to form the distribution to sample from
            tempsamp= exp(-P0(var,:)); % survival cost
            tempsamp(currsoln([1:var-1,var+3:end]))= 0;
            tempsamp1 = exp(-P0(var+1,:)); % spawn 1 cost
            tempsamp1(currsoln([1:var-1,var+3:end]))= 0;
            tempsamp2 = exp(-P0(var+2,:)); % spawn 2 cost
            tempsamp2(currsoln([1:var-1,var+3:end]))= 0;
            % form the overall distribution
            idxold_survival = find(tempsamp>0) ; 
            idxold_spawn1 = find(tempsamp1>0) ; % the first non-zero cost is the miss detection one
            idxold_spawn2 = find(tempsamp2>0) ; % the first non-zero cost is the miss detection one
            tempsamp1 = tempsamp1(idxold_spawn1) ; 
            tempsamp2 = tempsamp2(idxold_spawn2) ; 
            LTP = length(idxold_survival) ; 
            LTP1 = length(idxold_spawn1) ; 
            LTP2 = length(idxold_spawn2) ;
            
            %vectorrization
            mapOld = zeros(2 , LTP1*LTP2) ; 
            mapOld(1,:) = repmat(idxold_spawn1 , [1 , LTP2]) ; 
            mapOld(2,:)= reshape(repmat(idxold_spawn2,[LTP1 , 1]) , [1 LTP1*LTP2]) ; 
            tempsamp1 = repmat(tempsamp1 , [1 , LTP2]) ; 
            tempsamp2 = reshape(repmat(tempsamp2,[LTP1 , 1]) , [1 LTP1*LTP2]) ; 
            invalidIndc = (mapOld(1,:)-mapOld(2,:) == 0) ;
            tempsampspawn = prod([tempsamp1 ; tempsamp2],1) ;
            tempsampspawn(invalidIndc) = 0 ; 
            all_tempsamp =  [tempsamp(idxold_survival) , tempsampspawn] ;       
            idxold= find(all_tempsamp>0); all_tempsamp= all_tempsamp(idxold);
            [~,tempsoln]= histc(rand(1,1),[0;cumsum(all_tempsamp(:))/sum(all_tempsamp)]); 
            % convert the solution to the correct format
            tempsoln = idxold(tempsoln) ; 
            if tempsoln>LTP % spawns
                oldIndc = mapOld(:,tempsoln-LTP) ; 
                currsoln(var) = var ; 
                currsoln(var+1) = oldIndc(1) ; 
                currsoln(var+2) = oldIndc(2) ; 
            else % no spawns
                currsoln(var) = idxold_survival(tempsoln) ; 
                currsoln(var+1) = var+1 ; % not existing flag
                currsoln(var+2) = var+2 ; % not existing flag
            end
            var = var + 3 ; 
        end
    end
    if sol>burn_in
        assignments(sol-burn_in,:)= currsoln;
    end
end
[C,~,~]= unique(assignments,'rows');
assignments= C;