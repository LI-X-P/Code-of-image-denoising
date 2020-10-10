close all;
clear;clc;
tic;
%% Read image and Initialization
dir = pwd;
input_image = imread(strcat(dir,'\image\pillar.jpg'));
[~,~,d] = size(input_image);
if d == 3
    im_gray=rgb2gray(input_image);
    im_gray_1=im2double(im_gray);
else
    im_gray_1=im2double(input_image);
end

Noise_density = 0.80;%amount of noise
im_noised = imnoise(im_gray_1,'salt & pepper',Noise_density);

[p,q]=size(im_noised);
im_denoised=0.63*ones(p+20,q+20);% Padding for edges
im_denoised(10:p+9,10:q+9)=im_noised; %Image at center
im_Omega=0.63*ones(p+20,q+20);% Padding for edges
im_Omega(10:p+9,10:q+9)=im_noised; %Image at center

%% Type-2 Fuzzy Identifier
im_denoised_pixels = zeros(p+20,q+20);%denoised image
im_noised_pixels = zeros(p+20,q+20);%mask showing location of SAP
for rIdx=10:q+9 % To scan rows
    for cIdx=10:p+9 % To scan col.
        M = 1;  % Window size
        N_init = 6; % initial number of good pixels
        S_max = 2; % upper bound of 'M'
        n_gp = N_init; % number of good pixels
        
        while (im_denoised(cIdx,rIdx)==0)||(im_denoised(cIdx,rIdx)==1)
            subImg = windowImg(im_denoised,cIdx,rIdx,M); % Vector containing pixel around possible currupted pixel
            subImg_onlyGP = zeros(size(subImg));
            [T_Threshold,ave_PI] = type2_MF(subImg);
            windowLen = length(subImg);
            idx = 1;
            GP=zeros(1,n_gp);%for storing good pixels
            for x=1:windowLen
                if  ave_PI(x)>=T_Threshold %px. is good if > threshold
                    GP(idx)=subImg(x); 
                    idx=idx+1;
                    subImg_onlyGP(x) = subImg(x);
                elseif  (subImg(x)~=0)&&(subImg(x)~=1)
                    GP(idx)=subImg(x);
                    idx=idx+1;
                    subImg_onlyGP(x) = subImg(x);
                end
                if idx==n_gp+1
                    break
                end
            end
            
            im_Omega((cIdx-M):(cIdx+M),(rIdx-M):(rIdx+M))= reshape(subImg_onlyGP,[2*M+1,2*M+1 ]);
            neta=length(find(GP));%number of good pixels found
            if (neta<n_gp) && (M< S_max)% if we find less good pixels, inc. window size
                M=M+1;
                continue
            elseif (neta<n_gp) && (M== S_max)% if we find less good pixels&window size is max
                n_gp=n_gp-1;% try with less no. of good pixels
                if n_gp<1%if still unable to find, then try with even larger max window size
                    S_max=S_max+1;
                    n_gp=1;
                end
                continue
            end
            break
        end
    end
end
%% Matrix Completion Denoiser
im_Omega=im_Omega(10:p+9,10:q+9);
array_Omega = im_Omega;
array_Omega(array_Omega~=0) = 1;
[Y]  = Denoiser(im_gray_1, im_Omega, 40, array_Omega);


% %% Plotting
figure;
subplot(1,3,1);
imshow(im_gray_1);
title('Input image');
subplot(1,3,2);
imshow(im_noised);
title('Input noisy image');
subplot(1,3,3);
imshow(Y);
title('Recovered Image');


