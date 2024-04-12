close all
clear
clc
run E:\matconvnet-1.0-beta25\matlab\vl_setupnn.m
addpath(genpath('Hysure'))
addpath(genpath('superpixel'))
 addpath(genpath('CNN_subspace'))
 addpath(genpath('FSTRD'))
addpath(genpath('GLORIA'))
downsampling_scale=4;
 
 aa=load('C:\Users\admin\Desktop\LTMR\Indian_pines_corrected.mat');
 aa=aa.indian_pines_corrected;
aa=aa(1:128,1:128,:);
cc=double(aa);
cc=max(cc,0);
x1= max(max(max(cc)));
dd=cc/x1;
S=double(dd);
  [M, N, L]=size(S);      


spec0 = load( 'E:\我的代码\information fusion 2017\comparision\AVIRIS_spec.ascii');
spec = spec0(:,1); % Center wavelength of AVIRIS
wavelength = [450 520; 520 600; 630 690; 760 900; 1550 1750; 2080 2350]; % ETM/Landsat
m_band = size(wavelength,1); % Number of MSI bands
F = zeros(m_band,224);

for b=1:6
    b_i = find(spec>wavelength(b,1),1);
    b_e = find(spec<wavelength(b,2),1,'last');
    F(b,b_i:b_e) = 1/(b_e+1-b_i);
end
F=F(:,1:L);

for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
end
 


    
[M,N,L] = size(S);

%  simulate LR-HSI
S_bar = hyperConvert2D(S);

sag=3;
psf        =    fspecial('gaussian',7,sag);
par.fft_B      =    psf2otf(psf,[M N]);
par.fft_BT     =    conj(par.fft_B);
s0=1;
par.H          =    @(z)H_z(z, par.fft_B, downsampling_scale, [M N],s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, downsampling_scale,  [M N],s0);
Y_h_bar=par.H(S_bar);

  
SNRh=25;
sigma = sqrt(sum(Y_h_bar(:).^2)/(10^(SNRh/10))/numel(Y_h_bar));
rng(10,'twister')
   Y_h_bar = Y_h_bar+ sigma*randn(size(Y_h_bar));
HSI=hyperConvert3D(Y_h_bar,M/downsampling_scale, N/downsampling_scale );
[G,~,~] = Construct_Toeplitz_G(size(S,1),size(S,2),7,sag,downsampling_scale);





  %  simulate HR-MSI
rng(10,'twister')
Y = F*S_bar;
SNRm=30;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
MSI=hyperConvert3D(Y,M,N);

para.K=400;
para.eta=[1e-2, 1e-4, 1e-2];
para.p=10;
t0=clock;
 Z7 = Gernalized_TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
  t7=etime(clock,t0)
 %[metirc7] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);
 [psnr7,rmse7, ergas7, sam7, uiqi7,ssim7,DD7,CC7]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);