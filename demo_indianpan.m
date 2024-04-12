close all
clear
clc
run E:\matconvnet-1.0-beta25\matlab\vl_setupnn.m
addpath(genpath('NSSR'))
addpath(genpath('IR_TenSR-main'))
addpath(genpath('TSVD'))
addpath(genpath('CSTF'))
 addpath(genpath('NSSR'))
addpath(genpath('HySure'))
 addpath(genpath('FSTRD'))
 addpath(genpath('CNNMF'))
  addpath(genpath('functions'))
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



spec0 = load( 'E:\ÎÒµÄ´úÂë\information fusion 2017\comparision\AVIRIS_spec.ascii');
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
  sz=[M N];
%  simulate LR-HSI
S_bar = hyperConvert2D(S);

sag=2;
sf=4;
psf        =    fspecial('gaussian',7,sag);
par.fft_B      =    psf2otf(psf,[M N]);
par.fft_BT     =    conj(par.fft_B);
s0=1;
par.H          =    @(z)H_z(z, par.fft_B, downsampling_scale, [M N],s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, downsampling_scale,  [M N],s0);
Y_h_bar=par.H(S_bar);

  
SNRh=30;
sigma = sqrt(sum(Y_h_bar(:).^2)/(10^(SNRh/10))/numel(Y_h_bar));
rng(10,'twister')
   Y_h_bar = Y_h_bar+ sigma*randn(size(Y_h_bar));
HSI=hyperConvert3D(Y_h_bar,M/downsampling_scale, N/downsampling_scale );



%  simulate HR-MSI
rng(10,'twister')
Y = F*S_bar;
SNRm=35;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
MSI=hyperConvert3D(Y,M,N);
%% CNNMF
 t0=clock;
 Z5 = CNMF_fusion(HSI,MSI,F,par);
  t5=etime(clock,t0)
 [psnr5,rmse5, ergas5, sam5, uiqi5,ssim5,DD5,CC5] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/downsampling_scale);

%% IR_TenSR
para   =  Parameters_setting1( downsampling_scale, psf ,'pavia',sz,s0 );
MSI_h={};HSI_h={};X={}; MSI2=MSI;HSI2=HSI;MSI_clean=MSI;HSI_clean=HSI; 
t0=clock;
Z7 = IR_TenSR(HSI,MSI,F,para);
% for i=1:para.IR   
%      MSI_h{i}=MSI_clean;  
%      HSI_h{i}=HSI_clean;
%      t0=clock;
%      X{i} = IR_TenSR(HSI2,MSI2,F,para);
%      S_bar1=hyperConvert2D(X{i});
%      HSI1=hyperConvert3D(para.H(S_bar1),M/downsampling_scale,N/downsampling_scale);
%      MSI1 = hyperConvert3D(F*S_bar1, M, N);
%      HSI2=HSI_h{i}-HSI1;
%      MSI2=MSI_h{i}-MSI1;
%      MSI_clean=MSI2;
%      HSI_clean=HSI2; 
%  end
%     Z7=0;
%  for ii=1:para.IR
%       Z7=Z7+X{ii};
%  end 
 t7=etime(clock,t0);
  [psnr7,rmse7, ergas7, sam7, uiqi7,ssim7,DD7,CC7]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);
 



 %% GTSVD
para.K=400;
para.eta=[1e-2, 1e-6, 1e-1];
para.p=10;
t0=clock;
 Z10 = Gernalized_TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
  t10=etime(clock,t0)
 %[metirc7] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);
 [psnr10,rmse10, ergas10, sam10, uiqi10,ssim10,DD10,CC10]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z10)), 0, 1.0/sf);
 

  
%% GSA
t0=clock;
Z1=GSA_wrapper(HSI,MSI,downsampling_scale);
t1=etime(clock,t0);
 [psnr1,rmse1, ergas1, sam1, uiqi1,ssim1,DD1,CC1] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z1)), 0, 1.0/sf);
%% NSSR
par.P=F;
par.w=size(S,1);
par.h=size(S,2);
par.eta2       =  1e-4;    % 0.03
    par.eta1       =   1e-4;
    par.mu         =  2e-4;   % 0.004
    par.ro         =   1.1; 
    par.Iter       =   26;
par.K          =    80;
par.lambda     =    0.001;
par.s0=s0;
t0=clock;
Z2     =    NSSR_HSI_SR1( HSI,MSI,S_bar, sf,par,sz,s0 );
Z2=hyperConvert3D(Z2,sz(1),sz(2));
t2=etime(clock,t0)
 [psnr2,rmse2, ergas2, sam2, uiqi2,ssim2,DD2,CC2] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z2)), 0, 1.0/sf);

%% CSTF
par1.W=260; par1.H=260;  par1.S=15; par1.lambda=1e-5;
BW       =    fspecial('gaussian',[7 1],2);
 BW1=psf2otf(BW,[M 1]);
 BH       =    fspecial('gaussian',[7 1],2);
 BH1=psf2otf(BH,[N 1]);
t0=clock;
 Z3 = CSTF_FUS(HSI,MSI,F,BW1,BH1,sf,par1,s0,S);
  t3=etime(clock,t0)
 [psnr3,rmse3, ergas3, sam3, uiqi3,ssim3,DD3,CC3] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z3)), 0, 1.0/sf);

%% TSVD
para.K=200;
para.eta=1e-3;
para.p=10;
t0=clock;
 Z4 = TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
  t4=etime(clock,t0)
 [psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z4)), 0, 1.0/sf);

%  %% GLORIA
% SNR = 100;
% t0=clock;
% [G,~,~] = Construct_Toeplitz_G(size(S,1),size(S,2),7,sag,sf);
% Z5 = GLORIA_simplified(hyperConvert2D(HSI),hyperConvert2D(MSI),F,G,size(MSI,1),size(MSI,2),par,'Mu',30/SNR);
% Z5=hyperConvert3D(Z5,M,N);
% t5=etime(clock,t0);
%   %[metirc5] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/downsampling_scale);
% [psnr5,rmse5, ergas5, sam5, uiqi5,ssim5,DD5,CC5]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/sf);
%  %a5=[metirc5([1,4,6,5,3]),t6]

 


 %% FSTRD
par.lambda = 2;% MSI term
par.tau = 1e-6;     % TV term
par.rho = 1;     % PAM term
par.beta = 0.1;  % penalty term
par.r = [4,100,4]; %[2,150,2]; %TR rank
BW       =    fspecial('gaussian',[7 1],sag);
BW1=psf2otf(BW,[M 1]);
BH       =    fspecial('gaussian',[7 1],sag);
BH1=psf2otf(BH,[N 1]);
t0=clock;
[Z6] = FSTRD(HSI,MSI,F,BW1,BH1,sf,par,s0,S);
t6=etime(clock,t0);
%[metirc6] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z6)), 0, 1.0/downsampling_scale);
 %a6=[metirc6([1,4,6,5,3]),t6]
 [psnr6,rmse6, ergas6, sam6, uiqi6,ssim6,DD6,CC6]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z6)), 0, 1.0/sf);

%   %% CNN_FUS
% para.gama=1.1;
% para.p=10;
% para.sig=8e-4;
% t0=clock;
%  [Z8]= CNN_Subpace_FUS( HSI, MSI,F,par.fft_B,downsampling_scale,S,para,1);
%  t8=etime(clock,t0);
%  metric8 = quality_assessment(double(im2uint8(S)), double(im2uint8(Z8)), 0, 1.0/downsampling_scale);
  
 
save('indan.mat','Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z10','HSI','MSI')