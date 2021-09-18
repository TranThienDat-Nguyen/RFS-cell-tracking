function tt = jtt2tt(jtt, model)
mat_dim = size(jtt.m,1) ; 
num_objects = mat_dim/model.x_dim ; 
tt = cell(num_objects, 1) ; 
blk_indc = 1 : model.x_dim ; 
for n = 1 : num_objects
    tt{n}.w = jtt.w ; 
    tt{n}.m = jtt.m(blk_indc,:) ; 
    tt{n}.P = jtt.P(blk_indc, blk_indc, :) ; 
    tt{n}.l = jtt.l{n} ; 
    tt{n}.st = jtt.st(n,:) ; 
    tt{n}.mode = jtt.mode(n,:) ; 
    blk_indc = blk_indc + model.x_dim ; 
end
    