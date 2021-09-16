function [paths , cost] = kshortestwrap_pred_model ( init_mat_kshort , k ) 
    if k == 0 
        paths  = [] ;
        return ; 
    end
  CMPad = gen_kshortmat ( init_mat_kshort ) ; 
  destination = size(CMPad , 1); 
  numassign = size(init_mat_kshort , 2) ;
  [paths , cost] = kShortestPath_any(CMPad,1,destination,k); %do k-shortest path
  for p=1:length(paths)
      if isequal(paths{p},[1 destination])
          paths{p}= [];  %0 indicates no nodes selected
      else
          paths{p}= paths{p}(2:end-1)-1; %strip dummy entry and finish nodes
          mod_ = mod(paths{p} , numassign) ; 
          paths{p}(mod_ == 0) = numassign ; 
          paths{p}(mod_ ~= 0) = mod_(mod_~=0) ; 
      end
  end
end