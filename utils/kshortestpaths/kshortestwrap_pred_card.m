function [paths , cost] = kshortestwrap_pred_card ( init_mat_kshort , k ) 
    if k == 0 
        paths  = [] ;
        return ; 
    end
  numtracks = size(init_mat_kshort , 1) ;
  numassign = size(init_mat_kshort , 2) ;
  CMPad = gen_kshortmat ( init_mat_kshort ) ; 
  destination = size(CMPad , 1); 
  [paths , cost] = kShortestPath_any(CMPad,1,destination,k); %do k-shortest path
  uidx = 1 ; 
  for p=1:length(paths)
      if isequal(paths{p},[1 destination])
          paths{p}= [];  %0 indicates no nodes selected
      else
          paths{p}= paths{p}(2:end-1)-1; %strip dummy entry and finish nodes
          % Convert into assignment format (similar to Gibbs sampling
          % output)
          mod_ = mod(paths{p} , numassign) ; 
          paths{p}(mod_ == 0) = numassign ; 
          paths{p}(mod_ ~= 0) = mod_(mod_~=0) ; 
      end
  end
end