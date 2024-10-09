function ress = gpuNUFFT_adj(a,bb)
% ress = gpuNUFFT_adj(a,bb)
% Performs adjoint gpuNUFFT 
% from k-space to image space 
%
% supports multi-channel data
%
% a  ... GpuNUFFT Operator
% bb ... k-space data
%        k x nChn
%
nChn = size(bb,2);
sens = a.sens;
kspace = bb(:,:);
kspace = [real(kspace(:))'; imag(kspace(:))'];
kspace = reshape(kspace,[2 a.params.trajectory_length nChn]);
ress = mex_gpuNUFFT_adj_atomic_f(single(kspace),(a.dataIndices),single(a.coords),(a.sectorDataCount),(a.sectorProcessingOrder),(a.sectorCenters(:)),single(a.densSorted),single(sens),a.params);
ress = squeeze(ress(1,:,:,:,:) + 1i*(ress(2,:,:,:,:)));

