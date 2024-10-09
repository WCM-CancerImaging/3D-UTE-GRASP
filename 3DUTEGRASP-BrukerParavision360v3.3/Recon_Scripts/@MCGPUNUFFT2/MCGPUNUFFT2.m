function [res] = MCGPUNUFFT2(k,w,sens);
% function m = gpuNUFFT(k,w,osf,wg,sw,imageDim,sens,varargin)
%
%     k -- k-trajectory, scaled -0.5 to 0.5
%          dims: 3 ... x, y and z
%                N ... # sample points
%                nCh ... # channels / coils
%     w -- k-space weighting, density compensation
%     osf -- oversampling factor (usually between 1 and 2)
%     wg -- kernel width (usually 3 to 7)
%     sw -- sector width to use
%     imageDim -- image dimensions [n n n] 
%     sens -- coil sensitivity data
%  res -- gpuNUFFT operator
%
%  A. Schwarzl, Graz University of Technology
%  F. Knoll, NYU School of Medicine
%  Modified by Li Feng, NYU
%%%%%%  mod(imageDim*osf,sw))  has to be 0   Li Feng

osf = 1.5; wg = 3; sw = 8;
atomic = true;
use_textures = true;    
balance_workload = true;
res.adjoint = 0;
imageDim=size(sens(:,:,:,1));
res.imageDim = imageDim;
res.normsqr=10;

[nn,nd,nt]=size(k);

for tt=1:nt
    tt;
    res.op{tt}.params.img_dims = uint32(imageDim);
    res.op{tt}.params.osr = single(osf);
    res.op{tt}.params.kernel_width = uint32(wg);
    res.op{tt}.params.sector_width = uint32(sw);
    res.op{tt}.params.trajectory_length = uint32(length(k));
    res.op{tt}.params.use_textures = use_textures;
    res.op{tt}.params.balance_workload = balance_workload;
    
    [res.op{tt}.dataIndices,res.op{tt}.sectorDataCount,res.op{tt}.densSorted,res.op{tt}.coords,res.op{tt}.sectorCenters,res.op{tt}.sectorProcessingOrder] = mex_gpuNUFFT_precomp_f(single(k(:,:,tt)),single(ones(1,nn)),res.op{tt}.params);
    res.op{tt}.atomic = atomic;
    res.op{tt}.verbose = false;
end
tmp = [real(sens(:))'; imag(sens(:))'];
res.b1 = reshape(tmp,[2 imageDim(1)*imageDim(2)*imageDim(3) size(sens,4)]);
res.w=permute(sqrt(w),[1,3,2]);
res = class(res,'MCGPUNUFFT2');
