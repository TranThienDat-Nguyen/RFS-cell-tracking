% def: this function creates cost matrix for k-shortest path from standard
% cost matrix (for Gibbs)
function cost_mat_kshort = gen_kshortmat ( init_mat_kshort )
if ~isempty(init_mat_kshort)
    numtracks = size(init_mat_kshort , 1) ;
    cost_mat_kshort = zeros ( 2 + numel(init_mat_kshort) ) ; % extra 2 for start, end nodes
    % Start putting elements into cost_mat_kshort
    numcols = size( init_mat_kshort , 2 ) ; 
        % From start node to first track
    cost_mat_kshort(1 , 2 : 2 + numcols-1) = init_mat_kshort(1,:) ; 
        % From track 2 to last track
    spanidx = 2 + numcols  ;
    for tidx = 2 : numtracks
        temp_mat = repmat(init_mat_kshort(tidx , :) , numcols , 1) ;
        span_col = spanidx : spanidx + numcols - 1 ; 
        span_row = span_col - numcols ;
        cost_mat_kshort(span_row , span_col) = temp_mat ; 
        spanidx = spanidx + numcols ; 
    end
        % From last track to end node
    cost_mat_kshort(end-numcols  : end-1 , spanidx) = eps ;
else
    cost_mat_kshort = [0 eps ; 0 0] ; 
end
        
    
