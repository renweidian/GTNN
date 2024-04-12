function  [Z] =  Gernalized_TSVD_FUS( HSI, MSI,R,FBm,sf,S,para)
K=para.K;
eta=para.eta;

p=para.p;
mu=1e-3;
patchsize=8;
overlap=4;

HSI3=Unfold(HSI,size(HSI),3);
[D,~]=svds(HSI3,p);

L1=size(R,2);




 bparams.block_sz = [patchsize, patchsize];
 bparams.overlap_sz=[overlap overlap];
[nr, nc,~]=size(MSI);
L=size(HSI,3);
num1=(nr-patchsize)/(patchsize-overlap)+1;
num2=(nc-patchsize)/(patchsize-overlap)+1;
 bparams.block_num=[num1 num2]; 

predenoised_blocks = ExtractBlocks(MSI, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
if K==1
    aa=ones(num1*num2,1);
else
  [aa ]=fkmeans(Y2,K);
end

 HSI_int=zeros(nr,nc,L);
    HSI_int(1:sf:end,1:sf:end,:)=HSI;
    FBmC  = conj(FBm);
    FBs  = repmat(FBm,[1 1 L]);
       FBs1  = repmat(FBm,[1 1 L1]);
       FBCs=repmat(FBmC,[1 1 L]);
FBCs1=repmat(FBmC,[1 1 L1]);
HHH=ifft2((fft2(HSI_int).*FBCs));
  HHH1=hyperConvert2D(HHH);




%% iteration

MSI3=Unfold(MSI,size(MSI),3);

n_dr=nr/sf;
n_dc=nc/sf;

HR_load1=imresize(HSI, sf,'bicubic');

V1=hyperConvert2D(HR_load1);
V2=V1;
V3=V1;



G1=zeros(size(V2));
G2=G1;
G3=G1;

CCC=R'*MSI3+HHH1;
 C1=R'*R+3*mu*eye(size(R,2)); 
 [Q,Lambda]=eig(C1);
Lambda=reshape(diag(Lambda),[1 1 L1]);
InvLbd=1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
B2Sum=PPlus(abs(FBs1).^2./( sf^2),n_dr,n_dc);
InvDI=1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));



for i=1:100
    HR_HSI3=mu*(V1+G1/(2*mu)+V2+G2/(2*mu)+V3+G3/(2*mu));
C3=CCC+HR_HSI3;
C30=fft2(reshape((Q\C3)',[nr nc L1   ])).*InvLbd;
temp  = PPlus_s(C30/( sf^2).*FBs1,n_dr,n_dc);
invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1; % The operation: C5bar- temp*(\lambda_j d Im+\Sum_i=1^d Di^2)^{-1}Dv^H)
VXF    = Q*reshape(invQUF,[nc*nc L1])';
A = reshape(real(ifft2(reshape(VXF',[nr nc L1   ]))),[nc*nc L1])'; 
%   [ZE1] = Sylvester(C1,psfY.B, sf,n_dr,n_dc,C3);  
    
Zt=hyperConvert3D(A,nr, nc );
% aa1=C1*ZE-C3+ creat_HSI_T(creat_HSI(ZE ,psfY),psfY);
% norm(aa1(:))

  rmse2(i)=getrmse(double(im2uint8(S)),double(im2uint8(Zt))) 

%% spatial similarties

  

 B1=hyperConvert3D(A-G1/(2*mu),nr, nc );
 predenoised_blocks1 = ExtractBlocks(B1, bparams);
 
  B2=hyperConvert3D(A-G2/(2*mu),nr, nc );
 predenoised_blocks2 = ExtractBlocks(B2, bparams);
 
  B3=hyperConvert3D(A-G3/(2*mu),nr, nc );
 predenoised_blocks3 = ExtractBlocks(B3, bparams);
 
   Z_block1=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));
  Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));
  Z_block3=zeros( bparams.block_sz(1), bparams.block_sz(2),L1, bparams.block_num(1)* bparams.block_num(2));

  eta2=1e-3;
  eta3=1e-3;
for mn=1:max(aa)
    gg=find(aa==mn);
 XES=predenoised_blocks1(:,:,:,gg);
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
     V1=Log_prox_tnn( XES, eta/2/mu );
V1=reshape(V1,[a b c d]); 
  Z_block1(:,:,:,gg)=V1;
  
   XES=predenoised_blocks2(:,:,:,gg);
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
    XES=permute(XES,[1,3,2]);
       V2=Log_prox_tnn( XES, eta2/2/mu );
       V2=permute(V2,[1,3,2]);
V2=reshape(V2,[a b c d]); 
  Z_block2(:,:,:,gg)=V2;
 
    XES=predenoised_blocks3(:,:,:,gg);
   [a, b, c, d ]=size(XES);
    XES = reshape(XES,[a*b c d]);
    XES=permute(XES,[2,3,1]);
       V3=Log_prox_tnn( XES, eta3/2/mu );
       V3=permute(V3,[3,1,2]);
V3=reshape(V3,[a b c d]); 
  Z_block3(:,:,:,gg)=V3;
end
    
V1= JointBlocks(Z_block1, bparams);
V1=hyperConvert2D(V1);
G1=G1+2*mu*(V1-A);

V2= JointBlocks(Z_block2, bparams);
V2=hyperConvert2D(V2);
G2=G2+2*mu*(V2-A);

V3= JointBlocks(Z_block3, bparams);
V3=hyperConvert2D(V3);
G3=G3+2*mu*(V3-A);
end
Z=A;