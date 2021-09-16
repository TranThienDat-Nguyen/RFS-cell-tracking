% Def: run modified k-shortest paths for to find repeated meas. I call the
% output uasses to make it consistent with Gibbs sampler.
function [uasses,costs]= kshortestwrap_meas_res(init_mat_kshort,k)

if k==0
    uasses= [];
    costs= [];
    return;
end
  CMPad = gen_kshortmat ( init_mat_kshort ) ; 
  numtracks = size(init_mat_kshort , 1) ;
  numassign = size(init_mat_kshort , 2) ;
  destination = size(CMPad , 1); 
  [paths,costs]= kShortestPath_any(CMPad,1,destination,k); %do k-shortest path
  uasses = zeros(k , numtracks) ; 
  uidx = 1 ; 
  for p=1:length(paths)
      if isequal(paths{p},[1 destination])
          paths{p}= [];  %0 indicates no nodes selected
      else
          paths{p}= paths{p}(2:end-1)-1; %strip dummy entry and finish nodes
          % Convert into assignment format (similar to Gibbs sampling
          % output)
          mod_ = mod(paths{p} , numassign) ; 
          death_idx = mod_ == 1 ; 
          mis_idx = mod_ == 2 ; 
          last_meas_idx = mod_ == 0 ;
          uasses(uidx , :) = mod_ - 2 ; % Minus 2 because I want to take the meas indices (take the death & mis detec out)
          uasses(uidx , death_idx) = -Inf ; 
          uasses(uidx , mis_idx) = 0 ; 
          uasses(uidx , last_meas_idx) = numassign - 2 ;
          uidx = uidx + 1;
      end
  end
  % Prune the uasses
  uasses(uidx : end , :) = [] ; 
 
 