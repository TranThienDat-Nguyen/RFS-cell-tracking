function gen_synthesized_image(model, ccd_noise, path, truth)
img_size =  [model.range_c(1,2) model.range_c(2,2)] ;
nucli_radius = 3 ;
shape_param = [0.1 0.1] ; 
texture_param_normal = [0.5 2 5 0.2] ; 
texture_param_spawn = [5 5 5 5] ; 
illumscale =1 ; 
autofluorscale = 0.05 ;
O_S = 5;
O_V = 0.5;
ccd = ccd_noise ;
for k = 1 : 100
    img = zeros(img_size) ; 
    bw = img ; 
    for ind = 1 : truth.N(k)
        Y = truth.X{k}(3,ind) ; 
        X = truth.X{k}(1,ind) ;
        if ismember(truth.track_list{k}(ind), truth.spawnNextTimeTrack{k})
            n = nucleus([Y,X] , ind , nucli_radius , shape_param , texture_param_spawn) ; 
        else
            n = nucleus([Y,X] , ind , nucli_radius , shape_param , texture_param_normal) ; 
        end
        [img , bw] = object2image(n , img , bw) ; 
        % compensate for overlapping cells
        tmp = bw;
        tmp(tmp == 0) = Inf;
        img = img./tmp;
    end
    % Generate ideal image

    bkg_autofluor = autofluor(img_size,1,1,2,0);
    ba = sum(bkg_autofluor(:).^2);

    %Uneven illumination
    bkg_ill = illumination(fliplr(img_size),[0 , 0]);
    bi = sum(bkg_ill(:).^2);
    ie = sum(img(:).^2);

 
    ill_scale = sqrt(illumscale*ie/bi);
    autofluor_scale = sqrt(autofluorscale*ie/ba);

    blurred = my_optics(O_S,O_V,img,bw);
    
    img = imnoise(blurred + ill_scale*bkg_ill + autofluor_scale*bkg_autofluor,....
            'gaussian',0,ccd);
    imwrite(img , fullfile(path , [num2str(k,'%03.f'),'.tif'])) ;    
end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTICS Function for generating optical aberrations
% Input:   (1) scale of one object
%          (2) variance which is used for defining different focus levels
%          (3) struct containing images for all objects
%          (4) struct containing binary images for all objects
% Output:  (1) struct for blurred images
%
% (C) Antti Lehmussola, 22.2.2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[blurred] = my_optics(R,v,image, bw)

blurred = [];


%Define the overall area of cell
if ~isempty(bw)
	objects = bw;
end


L = bwlabel(objects);
D = zeros(size(objects));

%Define size of the gaussian kernel, here we use half of object scale
wl = round(R/2);

%Make size odd
if mod(wl,2) == 0
    wl = wl + 1;
end

    
%Generate depth information

%Depth bias
db = 0.0;

for ind = 1:max(L(:))
    r = db+randn*v; 
    D(L == ind) = abs(r);
end

%Blur images
if v > 0
	if ~isempty(image)
		blurred = blurimage(image,D,wl,.5,5);
	end
	
else
	blurred = image;
end
end
