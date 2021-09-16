%% provided by Seyed Hamid Rezatofighi (The Australian National University)
% based on the the following paper with some modifications
% S. H. Rezatofighi, R. Hartley, and W. E. Hughes, “A new approach for 
% spot detection in total internal reflection fluorescence microscopy,”
% ISBI 2012, pp. 860--863.
 
%% MPHD enhancement for a single 2D/3D image.
% The output of the MPHD_Enhancer function is the enhanced spots and the estimated background
 
% The main inputs of this function are:
% 1- I = a single 2D/3D image which includes bright spots and a non-homogeneous background
% 2- Zigma = The standard deviation of the Gaussian kernel. To avoid eliminating small objects,
% Zigma may be set by the size of the smallest spot.
% 3- Rmx = The radius Rmx is set by the user to limit the search area around
% each regional maximum and should be bigger than the size of the biggest
% spot. 
 
% The optional inputs are:
% 4- Hnd = This value makes the regional maxima with height< = Hnd flat. By default, Hnd=0. 
% To discard these regional maxima, this value can be changed.  
% 5- Hp =  Using this value, the regional maxima with height> = Hp is only 
% considered as candidate maxima. By default, Hp=1. 
% 6- Rmn = The radius Rmn is used to search for optimal base point in the area
% between Rmn and Rmx for a regional maximum. By default, Rmn=0.
 
% [Spots,Background] = MPHD_Enhancer(I,Zigma,Rmx,Hnd,Hp,Rmn);
[Spots,Background] = MPHD_Enhancer(I,Zigma,Rmx);
 
%% MPHD enhancement for 2D+t/3D+t sequences.
% Sicne this function uses temporal information for improving the results,
% this function is recommended for the 2D+t/3D+t sequences.
 
% The outputs of the MPMD_Enhancement_Sequences function are:
% 1- Temp_Spots = The enhanced spots using temporal averaging (the main
% output).
% 2- Spots = The enhanced spots without using temporal information. This
% output is exactly equal to "Spots" in the MPHD_Enhancer function for each frame.
% 3- Temp_Background = The estimated background using temporal averaging (the main
% output). Temp_Background + Temp_Spots = Is (The smoothed version of
% original sequenses).
% 4- Smoothed_Background (optional output) = The smoothed version of 
% Temp_Background using a Gaussian kernel with standard deviation Zigma_B.
 
% The main inputs of this function are:
% 1- I = the main 2D+t/3D+t sequences which includes bright spots and a non-homogeneous background
% 2- MPHD_P = MPHD parameters including [Zigma,Rmx,Hnd(optional),Hp(optional),Rmn(optional)];
% 3- AV_rate = Supposing the background’s characteristics change gradually
% and the background is constant over short periods of time = AV_rate, this value is used for
% temporal averaging. 
 
% The optional input is:
% 4- Zigma_B = for smoothing Temp_Background using a Gaussian kernel with
% standard deviation Zigma_B. By default, Zigma_B=2;
 
 
% [Temp_Spots,Spots,Temp_Background,Smoothed_Background]=MPMD_Enhancement_Sequences(I,MPHD_P,AV_rate,Zigma_B);
[Temp_Spots,Spots,Temp_Background,~]=MPMD_Enhancement_Sequences(I,MPHD_P,AV_rate);
 
%% Spot detection using the enhanced image.
% This function uses the enhanced images and the estimated background and 
% different thresholds to detect the objects of interest. This function works 
% on a single 2D/3D frame. For the sequences, it should be called for each
% frame separately
 
% The outputs of the Spot_Detector_MPHD function are:
% 1- centroids = The center of each detected spot.
% 2- Iz = The intensity of the each spot that can be either its maximum or
% its mean intensity. 
 
% The main inputs of this function are:
% 1- Spots = The "Spots" in the MPHD_Enhancer function or
% one of the frame of the "Spots" in the MPMD_Enhancement_Sequences function. 
% 2- Is = The smoothed 2D/3D image using the Gaussian kernel with the 
% standard deviation = Zigma. For the MPHD_Enhancer function, Is = Spots + Background
% and for the MPMD_Enhancement_Sequences function, Is = one frame of Temp_Spots + Temp_Background
% 3- Background = The estimated background. For the MPHD_Enhancer function,
% Background = Background. For the MPMD_Enhancement_Sequences function, use
% a single frame of Temp_Background.
% treI = A threshold (single value) on intensity of the enhanced spots that 
% is "Spots" (for the MPHD_Enhancer function) and "Temp_Spots" 
% (for the MPMD_Enhancement_Sequences function)
 
% The optional inputs are:
% 4- Intensity_mode = The intensity each spot can be calculated using the  
% maximum ('Max') or the mean ('Mean') intensity of each spot's region. By 
% default, Intensity_mode_B='Max';
% 5- treIb = A threshold (single value) on intensity of the enhanced spots
% without temporal averaging that is "Spots" (for both the MPHD_Enhancer and 
% MPMD_Enhancement_Sequences functions). By default, threIb=0;
% 6- treG = A threshold (single value) on the averaged gradient of the each 
% enhanced spot that is "Spots" (for the MPHD_Enhancer function) and "Temp_Spots" 
% (for the MPMD_Enhancement_Sequences function). By default, treG=0;
 
% For a single 2D/3D image
% [centroids,Iz]=Spot_Detector_MPHD(Spots,Spots+Background,Background,treI);
% [centroids,Iz]=Spot_Detector_MPHD(Spots,Spots+Background,Background,treI,'Mean',treI,treG);
 
% For a 2D+t movie
for f=1:Frame
    % [centroids{1,f},Iz{1,f}]=Spot_Detector_MPHD(Spots(:,:,f),Temp_Spots(:,:,f)+Temp_Background(:,:,f),Temp_Background(:,:,f),treI,'Mean',treI,treG);
    [centroids{1,f},Iz{1,f}]=Spot_Detector_MPHD(Spots(:,:,f),Temp_Spots(:,:,f)+Temp_Background(:,:,f),Temp_Background(:,:,f));
end
 


