clear
clc
addpath(genpath('CSU1'))
addpath(genpath('TSVD'))
addpath(genpath('CSTF'))
 addpath(genpath('NSSR'))
F=create_F();
sf = 16;
pathstr = fileparts('E:\STC\CAVE\');
dirname = fullfile(pathstr, '*.mat');
imglist = dir(dirname);

sz=[512 512];
s0=1;
 psf        =    fspecial('gaussian',17,3);
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

dd=1;
im_structure =load(fullfile(pathstr,  imglist(dd).name));
S = im_structure.b; 
S=S(:,:,3:31);
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
Y_h = hyperConvert3D(hyper, M/sf, N/sf);
Y = hyperConvert3D((F*S_bar), M, N);

%% p
para.K=160;
para.eta=1e-2;
para.p=10;
for yy=1:8
    para.p=2*yy+2;
 Z = TSVD_Subpace_FUS(Y_h,Y,F, par.fft_B,sf,S,para);
 [psnr1(yy),rmse1(yy), ergas1(yy), sam1(yy), uiqi1(yy),ssim1(yy),DD1(yy),CC1(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z)), 0, 1.0/sf);
end


para.K=160;
para.eta=1e-2;
para.p=10;
for yy=1:7
    para.eta=10^(yy-7);
 Z = TSVD_Subpace_FUS(Y_h,Y,F, par.fft_B,sf,S,para);
 [psnr2(yy),rmse2(yy), ergas2(yy), sam2(yy), uiqi2(yy),ssim2(yy),DD2(yy),CC2(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z)), 0, 1.0/sf);
end


para.K=160;
para.eta=1e-2;
para.p=10;
for yy=1:8
    para.K=1+50*(yy-1);
 Z = TSVD_Subpace_FUS(Y_h,Y,F, par.fft_B,sf,S,para);
 [psnr3(yy),rmse2(yy), ergas2(yy), sam2(yy), uiqi2(yy),ssim2(yy),DD2(yy),CC2(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z)), 0, 1.0/sf);
end

save('para_cave_32.mat','psnr1','psnr2','psnr3')