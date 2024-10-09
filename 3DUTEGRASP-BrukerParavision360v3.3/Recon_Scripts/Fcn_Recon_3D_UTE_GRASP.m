function [ ] = Fcn_Recon_3D_UTE_GRASP( Mouse, ScanPara, ReconPara, iEcho)
% Fcn_Recon_3D_UTE_GRASP Iterative reconstruction using 3D-UTE-GRASP
%
%   Input: 
%       Mouse: Mouse property structure 
%       ScanPara: Scan parameter structure
%       ReconPaara: Reconstruction parameter structure
% 
%   Output: 
%       reconstruction image and analyze files are saved to local directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of coil sensitivite maps

load([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_Rawdata]);

traj=reshape(k,[ScanPara.nTrajSamples, ScanPara.nTotalProj, 3]);
nFrame=floor(ScanPara.nTotalProj/(ReconPara.T/ScanPara.TR));

load([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_DCF]);
DensityCompen3D=[];
for i=1:1:nFrame
    DensityCompen3D=[DensityCompen3D squeeze(DCF.frame(i,:,:))];
end

clear kdata
kdata=squeeze(rawdata(:,:,:,iEcho));
clear rawdata

nt=1;%Here,define the temporal resolution, for calcuation of b1 maps, use nt=1
nline=floor(size(DensityCompen3D,2)/nt);
for ii=1:nt
    Traj3D1(:,:,:,ii)=traj(:,(ii-1)*nline+1:ii*nline,:);
    DensityCompen3D1(:,:,ii)=DensityCompen3D(:,(ii-1)*nline+1:ii*nline);
    kdata1(:,:,:,ii)=kdata(:,(ii-1)*nline+1:ii*nline,:);
end
[nx,ntviews,nc,nt]=size(kdata1);

Traj3D1=single(reshape(Traj3D1,[nx*ntviews,3,nt]));
DensityCompen3D1=(single(reshape(DensityCompen3D1,[nx*ntviews,nt])));
kdata1=single(reshape(kdata1,[nx*ntviews,nc,nt]));
y=kdata1.*permute(repmat(sqrt(DensityCompen3D1),[1,1,nc]),[1,3,2]);

E=MCGPUNUFFT2(single(Traj3D1),single(DensityCompen3D1),ones(128,128,128));
for ii=1:size(y,2)
    ref(:,:,:,ii)=E'*y(:,ii);
end
%size(ref)
b1=adapt_array_3d(ref);
b1=b1/max(abs(b1(:)));

MouseEchoFolder=fullfile(Mouse.PathBase, Mouse.PathRecon, Mouse.ID, ['echo' num2str(iEcho, '%1d')]);
if ~exist(MouseEchoFolder, 'dir')
    mkdir(MouseEchoFolder);
end

save(fullfile(MouseEchoFolder, Mouse.SavedFile_CoilSens), 'b1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D UTE GRASP Reconstruction
clearvars -except nFrame Mouse ScanPara ReconPara MouseEchoFolder iEcho;
load([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_Rawdata]);
traj=reshape(k,[ScanPara.nTrajSamples, ScanPara.nTotalProj, 3]);

load(fullfile(MouseEchoFolder, Mouse.SavedFile_CoilSens));

load([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_DCF]);
DensityCompen3D=[];
for i=1:1:nFrame
    DensityCompen3D=[DensityCompen3D squeeze(DCF.frame(i,:,:))];
end

clear kdata
kdata=squeeze(rawdata(:,:,:,iEcho));
clear rawdata

nt=nFrame;%Here,define the temporal resolution
nline=floor(ReconPara.T/ScanPara.TR);
clear Traj3D1 DensityCompen3D1 kdata1
for ii=1:nt
    Traj3D1(:,:,:,ii)=traj(:,(ii-1)*nline+1:ii*nline,:);
    DensityCompen3D1(:,:,ii)=DensityCompen3D(:,(ii-1)*nline+1:ii*nline);
    kdata1(:,:,:,ii)=kdata(:,(ii-1)*nline+1:ii*nline,:);
end
[nx,ntviews,nc,nt]=size(kdata1);

Traj3D1=single(reshape(Traj3D1,[nx*ntviews,3,nt]));
DensityCompen3D1=(single(reshape(DensityCompen3D1,[nx*ntviews,nt])));
kdata1=single(reshape(kdata1,[nx*ntviews,nc,nt]));
param.y=kdata1.*permute(repmat(sqrt(DensityCompen3D1),[1,1,nc]),[1,3,2]);
param.E=MCGPUNUFFT2(single(Traj3D1),single(DensityCompen3D1),b1);
recon_cs=param.E'*param.y;

Max_NUFFT=recon_cs;

save(fullfile(MouseEchoFolder, Mouse.SavedFile_NUFFT), '-v7.3', 'recon_cs');

param.TV=TV_Temp3D;

param.TVWeight=max(abs(recon_cs(:)))*ReconPara.lambda/100;

param.L1Weight=0;
%param.L1Weight=max(abs(recon_cs(:)))*ReconPara.lambda/100;

param.TVWeightRes=0;
param.nite =5;
param.display = 1;

param.initial=recon_cs;

n=1;
bStop=0;
error=1;
while ~bStop
    [recon_cs, recon_fval]= CSL1NlCg_4DRadial(recon_cs, param);
    cost_f_val(n)=recon_fval;
    save(fullfile(MouseEchoFolder, [Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_nIter_' num2str(n,'%.3d') '.mat']), '-v7.3', 'recon_cs', 'recon_fval', 'cost_f_val', 'Mouse', 'ScanPara', 'ReconPara');
    if n>=2
        error=(cost_f_val(n)-cost_f_val(n-1))/cost_f_val(n-1);
    end
    if abs(error)<2.5/100
        bStop=1;
    end    
    n=n+1;
end

recon_cs=abs(recon_cs);
save(fullfile(MouseEchoFolder, Mouse.SavedFile_GRASP), '-v7.3', 'recon_cs', 'recon_fval', 'cost_f_val', 'Mouse', 'ScanPara', 'ReconPara');
delete(fullfile(MouseEchoFolder, [Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_nIter*.mat']));

% save the 3D image
pixdim=ScanPara.FOV./ScanPara.dim;
ImgWhole=squeeze(abs(recon_cs(:,:,:,end))); % using the last frame
desPath=MouseEchoFolder;
desFile=[Mouse.SavedFile_Analyze];
dtype='int32';
origin=[0,0,0];
Fcn_write_analyze(ImgWhole, ScanPara.dim, pixdim, desPath, desFile, dtype, origin);


DynamicFrame.D1=squeeze(recon_cs(:,round(ScanPara.dim(1)/2),:,:));
DynamicFrame.D2=squeeze(recon_cs(round(ScanPara.dim(2)/2),:,:,:));
DynamicFrame.D3=squeeze(recon_cs(:,:,round(ScanPara.dim(3)/2),:));

dim=[size(recon_cs,1), size(recon_cs,3), size(recon_cs,4)];
desFile=[Mouse.SavedFile_Analyze_D1];
Fcn_write_analyze(DynamicFrame.D1, dim, pixdim, desPath, desFile, dtype, origin);

dim=[size(recon_cs,2), size(recon_cs,3), size(recon_cs,4)];
desFile=[Mouse.SavedFile_Analyze_D2];
Fcn_write_analyze(DynamicFrame.D2, dim, pixdim, desPath, desFile, dtype, origin);

dim=[size(recon_cs,1), size(recon_cs,2), size(recon_cs,4)];
desFile=[Mouse.SavedFile_Analyze_D3];
Fcn_write_analyze(DynamicFrame.D3, dim, pixdim, desPath, desFile, dtype, origin);

end
