function store_glmb(setting, glmb, output_path, k)
    glmb = prune(glmb,  setting.store_glmb_prune_thres) ; 
    glmb = cap(glmb, setting.store_glmb_cap) ; 
    if ~isfolder(fullfile(output_path, 'glmb'))
        mkdir(fullfile(output_path, 'glmb')) ; 
    end 
    filename = fullfile(output_path, 'glmb', [sprintf('%03d', k), '.mat']) ; 
    save(filename, 'glmb') ; 
end

function glmb_out= prune(glmb_in,hyp_threshold)
%prune components with weights lower than specified threshold
idxkeep= find(glmb_in.w > hyp_threshold);
if length(idxkeep)<1
    [~,idxkeep] = sort(glmb_in.w, 'descend') ; 
    idxkeep = idxkeep(1) ; 
end
glmb_out.tt= glmb_in.tt;
glmb_out.w= glmb_in.w(idxkeep);
glmb_out.I= glmb_in.I(idxkeep);
glmb_out.n= glmb_in.n(idxkeep);
glmb_out.total_spawns = glmb_in.total_spawns(idxkeep) ; 
glmb_out.total_births = glmb_in.total_births(idxkeep) ; 
glmb_out.N_clutter = glmb_in.N_clutter(idxkeep) ;
glmb_out.N_clutter_count= glmb_in.N_clutter_count(idxkeep);
glmb_out.w= glmb_out.w/sum(glmb_out.w);
for card=0:max(glmb_out.n)
    glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
end
for card = 0:max(glmb_out.total_spawns)
    glmb_out.cdn_spawn(card+1) = sum(glmb_out.w(glmb_out.total_spawns==card));
end
glmb_out.I_old = glmb_in.I_old(idxkeep) ; 
end



function glmb_out= cap(glmb_in,H_max)
%cap total number of components to specified maximum
if length(glmb_in.w) > H_max
    [~,idxsort]= sort(glmb_in.w,'descend');
    idxkeep=idxsort(1:H_max);
    glmb_out.tt= glmb_in.tt;
    glmb_out.w= glmb_in.w(idxkeep);
    glmb_out.I= glmb_in.I(idxkeep);
    glmb_out.I_old= glmb_in.I_old(idxkeep);
    glmb_out.n= glmb_in.n(idxkeep);
    glmb_out.n= glmb_in.n(idxkeep);
    glmb_out.total_spawns = glmb_in.total_spawns(idxkeep) ; 
    glmb_out.total_births = glmb_in.total_births(idxkeep) ; 
    glmb_out.N_clutter = glmb_in.N_clutter(idxkeep) ;
    glmb_out.N_clutter_count= glmb_in.N_clutter_count(idxkeep);
    glmb_out.w= glmb_out.w/sum(glmb_out.w);
    for card=0:max(glmb_out.n)
        glmb_out.cdn(card+1)= sum(glmb_out.w(glmb_out.n==card));
    end
    for card = 0:max(glmb_out.total_spawns)
        glmb_out.cdn_spawn(card+1) = sum(glmb_out.w(glmb_out.total_spawns==card));
    end
else
    glmb_out= glmb_in;
end
end