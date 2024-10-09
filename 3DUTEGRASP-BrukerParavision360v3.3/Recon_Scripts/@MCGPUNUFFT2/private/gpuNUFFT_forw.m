function ress = gpuNUFFT_forw(a,bb)
% ress = gpuNUFFT_forw(a,bb)
% Performs forward gpuNUFFT 
% from image to k-space
%
% supports multi-channel data
%
% a  ... GpuNUFFT Operator
% bb ... image data
%        W x H x D x (nChn) for 3d 
%        W x H x (nChn)     for 2d  
%
nChn = size(bb,4);
bb = [real(bb(:))'; imag(bb(:))'];
bb = reshape(bb,[2 a.params.img_dims(1)*a.params.img_dims(2)*a.params.img_dims(3) nChn]);
sens = a.sens;
nChn = a.sensChn;
m = mex_gpuNUFFT_forw_atomic_f(single(bb),(a.dataIndices),single(a.coords),(a.sectorDataCount),(a.sectorProcessingOrder),(a.sectorCenters(:)),single(sens),a.params);
ress(:,:) = squeeze(m(1,:,:) + 1i*(m(2,:,:)));


