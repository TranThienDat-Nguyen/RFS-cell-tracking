function plot_lineage(est, id)
    % generate a tree for tree plot
    % find fami
    fam_ances = zeros(2,0) ; 
    fam = cell(1,0) ;
    % init_ances = cell(1,0) ; 
    for k = 1 : length(est.L)
        for t = 1 : length(est.L{k})
            lab = est.L{k}{t} ; 
            ances_lab = lab([1:2],1) ; 
            curr_fam_idx = find_vector_in_matrix(ances_lab, fam_ances) ; 
            if curr_fam_idx == 0
                fam_ances = [fam_ances, ances_lab] ; 
                fam{length(fam)+1} = cell(1,0) ; 
            else
                temp_idx = find_vector_in_cell(lab, fam{curr_fam_idx}) ; 
                if temp_idx == 0
                    fam{curr_fam_idx}{length(fam{curr_fam_idx})+1} = lab ; 
                end
            end
        end
    end
    nodes_fam = cell(1,length(fam)) ; 
    for fidx = 1 : length(fam)
        fam_temp = fam{fidx} ; 
        for n = 1 : length(fam{fidx})
            fam_temp{n} = sprintf('%.0f-' , fam_temp{n});
        end
        nodes_fam{fidx} = form_nodes(fam_temp) ; 
    end
    my_treeplot(nodes_fam, id) ; 
end


function nodes = form_nodes(unique_labs)
nodes = zeros(1, length(unique_labs)) ;
if length(unique_labs)>1
    for lidx = 1 : length(unique_labs)
        % split the strings
        temp = unique_labs{lidx} ; 
        temp(end) = [] ; 
        temp = split(temp,'-') ; 
        if length(temp) == 2
            nodes(lidx) = 0 ;
        else
            L = length(temp) - 2 ;
            target_str =[] ; 
            for ii = 1 : L
                target_str = [target_str, temp{ii}, '-'] ; 
            end
            cnt = 1 ; 
            not_found = 1 ;
            while cnt<=length(unique_labs) && not_found
                if strcmp(target_str, unique_labs(cnt))
                    not_found = 0 ; 
                    nodes(lidx) = cnt ; 
                else
                    cnt = cnt + 1 ;
                end
            end
        end
    end
end
end

function idx = find_vector_in_matrix(v,m)
idx = 0 ;
L = size(m,2) ;
ii = 1 ; 
while ii<=L && idx==0
    if all(v-m(:,ii) == 0)
        idx = ii ; 
    end
    ii = ii + 1 ;
end
end

function idx = find_vector_in_cell(v,c)
idx = 0 ; 
L = length(c) ; 
ii = 1 ; 
while ii<=L && idx==0
    if length(v) == length(c{ii})
        if all(v-c{ii}==0)
            idx = ii ; 
        end
    end
    ii = ii + 1 ; 
end
end