
clear
clc
addpath(genpath('CSU1'))
addpath(genpath('TSVD'))
addpath(genpath('CSTF'))
 addpath(genpath('NSSR'))
 addpath(genpath('HySure'))
F=create_F();
sf =32;
pathstr = fileparts('E:\调参数据集\高光谱集数据\CAVE\');
dirname = fullfile(pathstr,'*.mat');
imglist = dir(dirname);

sz=[512 512];
s0=1;
 psf        =    fspecial('gaussian',33,4);
  par.fft_B      =    psf2otf(psf,sz);
  par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);
F=F(:,3:31);
 for band = 1:size(F,1)
        div = sum(F(band,:));
        for i = 1:size(F,2)
            F(band,i) = F(band,i)/div;
        end
    end
par.P=F;

dd=30;

% for yy=dd
% im_structure =load(fullfile(pathstr, 'shuju1', imglist(yy).name));
% S = im_structure.b; 
% S=S(:,:,3:31);
% [M,N,L] = size(S);
% S_bar = hyperConvert2D(S);
% hyper= par.H(S_bar);
% multi=F*S_bar;
%  HSI= hyperConvert3D(hyper, M/sf, N/sf);
%  MSI = hyperConvert3D((F*S_bar), M, N);
%  basis_type = 'VCA';
% lambda_phi = 5e-5;
% lambda_m = 1;
% p=10;
% B=ifft2( par.fft_B );
% t0=clock;
%  Z1=data_fusion( HSI, MSI, sf, F, B, p, basis_type, lambda_phi, lambda_m,s0 );
% t1(yy)=etime(clock,t0)
%  [psnr1(yy),rmse1(yy), ergas1(yy), sam1(yy), uiqi1(yy),ssim1(yy),DD1(yy),CC1(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z1)), 0, 1.0/sf);
%  end



%% CSU
% for yy=dd
% im_structure =load(fullfile(pathstr, 'shuju1', imglist(yy).name));
% S = im_structure.b; 
% S=S(:,:,3:31);
% [M,N,L] = size(S);
% S_bar = hyperConvert2D(S);
% hyper= par.H(S_bar);
% multi=F*S_bar;
% par.w=size(S,1);
% par.h=size(S,2);
% p=10;
% t0=clock;
% % [E,A] = SupResPALM(hyper, multi, S_bar, F,p,par);
%  Z1 = hyperConvert3d(E*A);
% t1(yy)=etime(clock,t0)
%  [psnr1(yy),rmse1(yy), ergas1(yy), sam1(yy), uiqi1(yy),ssim1(yy),DD1(yy),CC1(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z1)), 0, 1.0/sf);
%  end

%% NSSR
% for yy=dd
% im_structure =load(fullfile(pathstr, 'shuju1', imglist(yy).name));
% S = im_structure.b;
% S=S(:,:,3:31);
% [M,N,L] = size(S);
% S_bar = hyperConvert2D(S);
% hyper= par.H(S_bar);
% Y_h = hyperConvert3D(hyper, M/sf, N/sf);
% Y = hyperConvert3D((F*S_bar), M, N);
% par.P=F;
% par.w=size(S,1);
% par.h=size(S,2);
% par.eta2       =  1e-4;    % 0.03
%     par.eta1       =   1e-2;
%     par.mu         =  2e-4;   % 0.004
%     par.ro         =   1.1; 
%     par.Iter       =   26;
% par.K          =    80;
% par.lambda     =    0.001;
% par.s0=s0;
% t0=clock;
% Z2     =    NSSR_HSI_SR1( Y_h,Y,S_bar, sf,par,sz,s0 );
% Z2=hyperConvert3D(Z2,sz(1),sz(2));
% t2(yy)=etime(clock,t0)
%  [psnr2(yy),rmse2(yy), ergas2(yy), sam2(yy), uiqi2(yy),ssim2(yy),DD2(yy),CC2(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z2)), 0, 1.0/sf);
%  end



%% CSTF
%  for yy=dd
%     yy
% im_structure =load(fullfile(pathstr, 'shuju1', imglist(yy).name));
% S = im_structure.b;  
% S=S(:,:,3:31);
% [M,N,L] = size(S);
% NN=M/sf;
% S_bar = hyperConvert2D(S);
% hyper= par.H(S_bar);
%   HSI =hyperConvert3D(hyper,NN, NN );
% MSI = hyperConvert3D((F*S_bar), M, N);
% par1.W=260; par1.H=260;  par1.S=15; par1.lambda=1e-5;
% BW       =    fspecial('gaussian',[7 1],2);
%  BW1=psf2otf(BW,[M 1]);
%  BH       =    fspecial('gaussian',[7 1],2);
%  BH1=psf2otf(BH,[N 1]);
% t0=clock;
%  Z3 = CSTF_FUS(HSI,MSI,F,BW1,BH1,sf,par1,s0,S);
%   t3(yy)=etime(clock,t0)
%  [psnr3(yy),rmse3(yy), ergas3(yy), sam3(yy), uiqi3(yy),ssim3(yy),DD3(yy),CC3(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z3)), 0, 1.0/sf);
%  end






%% TSVD
 for yy=1:32
    yy
im_structure =load(fullfile(pathstr,  imglist(yy).name));
S = im_structure.b; 
S=S(:,:,3:31);
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
Y_h = hyperConvert3D(hyper, M/sf, N/sf);
Y = hyperConvert3D((F*S_bar), M, N);
para.K=300;
para.eta=3e-2;
para.p=10;
t0=clock;
 Z4 = TSVD_Subpace_FUS(Y_h,Y,F, par.fft_B,sf,S,para);
  t4(yy)=etime(clock,t0)
 [psnr4(yy),rmse4(yy), ergas4(yy), sam4(yy), uiqi4(yy),ssim4(yy),DD4(yy),CC4(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z4)), 0, 1.0/sf);
 end

for yy=1:32
    yy
im_structure =load(fullfile(pathstr,  imglist(yy).name));
S = im_structure.b; 
S=S(:,:,3:31);
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
HSI = hyperConvert3D(hyper, M/sf, N/sf);
MSI = hyperConvert3D((F*S_bar), M, N);
para.K=300;
para.eta=[3e-2, 1e-3, 2e-2];
para.p=10;
t0=clock;
 Z5 =  Gernalized_TSVD_Subpace_FUS(HSI,MSI,F, par.fft_B,sf,S,para);
  t5(yy)=etime(clock,t0)
 [psnr5(yy),rmse5(yy), ergas5(yy), sam5(yy), uiqi5(yy),ssim5(yy),DD5(yy),CC5(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/sf);
end
 
% QM.Z1=Z1;
% QM.Z2=Z2;
% QM.Z3=Z3;
% QM.Z4=Z4;
% QM.S=S;

% QM.psnr=[mean(psnr1), mean(psnr2),mean(psnr3),mean(psnr4)];
% QM.rmse=[mean(rmse1), mean(rmse2),mean(rmse3),mean(rmse4)];
% QM.t=[mean(t1), mean(t2),mean(t3),mean(t4)];
% QM.sam=[mean(sam1), mean(sam2),mean(sam3),mean(sam4)];
% QM.ergas=[mean(ergas1), mean(ergas2),mean(ergas3),mean(ergas4)];
% QM.uiqi=[mean(uiqi1), mean(uiqi2),mean(uiqi3),mean(uiqi4)];


% QM.psnr=[(uiqi1); (uiqi2);(uiqi3);(uiqi4)];
% QM.psnr=[(psnr1); (psnr2);(psnr3);(psnr4)];

save('cave_30.mat','QM')
