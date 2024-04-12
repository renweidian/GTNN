
clear
clc
addpath(genpath('NSSR'))
addpath(genpath('IR_TenSR-main'))
addpath(genpath('TSVD'))
addpath(genpath('CSTF'))
 addpath(genpath('NSSR'))
addpath(genpath('HySure'))
 addpath(genpath('FSTRD'))
 addpath(genpath('CNNMF'))
  addpath(genpath('functions'))
S=imread('E:\ÎÒµÄ´úÂë\LTMR\CSTF\CSTF\original_rosis.tif');
F=load('C:\Users\admin\Desktop\LTMR\CSTF\CSTF\R.mat');
S=S(1:256,1:256,11:end);
S=double(S);
S=S/max(S(:));
[M N L]=size(S);

sag=2;
sf =4;
downsampling_scale=sf;
sz=[M N];
s0=1;
 psf        =    fspecial('gaussian',7,sag);
  par.fft_B      =    psf2otf(psf,sz);
  par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);

% R=load('C:\Users\Utilizador\Desktop\drw\learningcode\HSI-superresolution\TGRS-2015\HySure-master1\HySure-master\data\ikonos_spec_resp.mat');
% R=R.ikonos_sp;
% [~, valid_ik_bands] = intersect(R(:,1), 430:860);
% no_wa = length(valid_ik_bands);
% xx  = linspace(1, no_wa, L);
% x = 1:no_wa;
% F = zeros(5, L);
% for i = 1:5 % 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
%     F(i,:) = spline(x, R(valid_ik_bands,i+1), xx);
% end
% F = F(2:4,:);

F=F.R;
F=F(:,1:end-10);
for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
end
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
MSI = hyperConvert3D((F*S_bar), M, N);
SNRh=35;
sigma = sqrt(sum(MSI(:).^2)/(10^(SNRh/10))/numel(MSI));
rng(10,'twister')
 MSI = MSI+ sigma*randn(size(MSI));

  HSI =hyperConvert3D(hyper,M/sf, N/sf );
SNRh=30;
sigma = sqrt(sum(HSI(:).^2)/(10^(SNRh/10))/numel(HSI));
rng(10,'twister')
 HSI = HSI+ sigma*randn(size(HSI));
 
  
%% IR_TenSR
para   =  Parameters_setting1( downsampling_scale, psf ,'pavia',sz,s0 );
MSI_h={};HSI_h={};X={}; MSI2=MSI;HSI2=HSI;MSI_clean=MSI;HSI_clean=HSI; 
t0=clock;
for i=1:para.IR   
     MSI_h{i}=MSI_clean;  
     HSI_h{i}=HSI_clean;
     t0=clock;
     X{i} = IR_TenSR(HSI2,MSI2,F,para);
     S_bar1=hyperConvert2D(X{i});
     HSI1=hyperConvert3D(para.H(S_bar1),M/downsampling_scale,N/downsampling_scale);
     MSI1 = hyperConvert3D(F*S_bar1, M, N);
     HSI2=HSI_h{i}-HSI1;
     MSI2=MSI_h{i}-MSI1;
     MSI_clean=MSI2;
     HSI_clean=HSI2; 
 end
    Z7=0;
 for ii=1:para.IR
      Z7=Z7+X{ii};
 end 
 t7=etime(clock,t0);
  [psnr7,rmse7, ergas7, sam7, uiqi7,ssim7,DD7,CC7]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);
 %% CNNMF
 t0=clock;
 Z5 = CNMF_fusion(HSI,MSI,F,par);
  t5=etime(clock,t0)
 [psnr5,rmse5, ergas5, sam5, uiqi5,ssim5,DD5,CC5] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/downsampling_scale);

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
  
  
  
 %% GTSVD
para.K=400;
para.eta=[1e-2, 1e-4, 1e-2];
para.p=10;
t0=clock;
 Z10 = Gernalized_TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
  t10=etime(clock,t0)
 %[metirc7] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z7)), 0, 1.0/sf);
 [psnr10,rmse10, ergas10, sam10, uiqi10,ssim10,DD10,CC10]=quality_assessment(double(im2uint8(S)), double(im2uint8(Z10)), 0, 1.0/sf);
 

 
save('pavia.mat','Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z10','HSI','MSI')
