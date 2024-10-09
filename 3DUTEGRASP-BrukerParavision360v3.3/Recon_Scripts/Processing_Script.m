% 3D-UTE-GRASP reconstruction demo code for Bruker PV360 (v3.3)

% Reference: 
% Magn Reson Med. 2019 Jan; 81(1):140-152. doi: 10.1002/mrm.27357. Epub 2018 Jul 29.
% <<Rapid dynamic contrast-enhanced MRI for small animals at 7T using 3D
% ultra-short echo time and golden-angle radial sparse parallel MRI>>
% Author: Zhang J, Feng L, Otazo R, Kim SG.

% Author: Jin Zhang, Li Feng, Ricardo Otazo, Sungheon Gene Kim 
% New Your University Langone Health
% Date: 08/01/2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to reconsturcte 3D-UTE-GRASP dynamic series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mouse.PathBase='/athena/kimlab/scratch/jzh4009/Recon_3DUTEGRASP_PV360/'; % Base Path 
Mouse.PathRecon='Recon_Data/'; % Recon data Path
% Mouse.ID='P230731_SingleEcho'; % single echo acquisition 
Mouse.ID='P230713_DualEcho'; % dual echo acquisition

% Reconstruction Parameter
ReconPara.T=05; % Temporal resolution 
ReconPara.lambda=1; % Reconstruction regularization parameter in % of max NUFFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saved file names
Mouse.SavedFile_Rawdata=[Mouse.ID '_rawdata.mat']; % k trajectory and rawdata from scan
Mouse.SavedFile_DCF=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_DCF.mat']; % density compenstation function
Mouse.SavedFile_CoilSens=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_CoilSens.mat']; % coil sensitivity
Mouse.SavedFile_NUFFT=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconNUFFT.mat']; % NUFFT recon
Mouse.SavedFile_GRASP=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconGRASP.mat']; % Grasp recon

Mouse.SavedFile_Analyze=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconGRASP_3D']; % 3D last frame analyze image

Mouse.SavedFile_Analyze_D1=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconGRASP_Coronal']; % Dynamic frames of center slice in Coronal
Mouse.SavedFile_Analyze_D2=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconGRASP_Sagittal']; % Dynamic frames of center slice in Sagittal
Mouse.SavedFile_Analyze_D3=[Mouse.ID '_T' num2str(ReconPara.T, '%.3d') '_ReconGRASP_Axial']; % Dynamic frames of center slice in Axial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan parameters
ScanPara.nReceiver=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'method', 'PVM_EncNReceivers');
ScanPara.nTotalProj=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'method', 'NPro');
ScanPara.nEchoes=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'acqp', 'NECHOES');
ScanPara.dim=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'method', 'PVM_Matrix');

% compute the trajectory points for PV360
ScanPara.nTrajSamples=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'method', 'PVM_TrajSamples');
ScanPara.nPostPoints=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'method', 'PostPoints');
ScanPara.nTrajSamples=ScanPara.nTrajSamples-ScanPara.nPostPoints;

ScanPara.TR=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'acqp', 'ACQ_repetition_time'); % in ms
ScanPara.TR=ScanPara.TR/1000; % in s

ScanPara.FOV=Fcn_Bruker_Method([Mouse.PathBase Mouse.PathRecon Mouse.ID], 'acqp', 'ACQ_fov'); % in cm
ScanPara.FOV=ScanPara.FOV*10; % in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Extract raw data and save as a matlab dat file
if ~isfile([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_Rawdata])
    Fcn_Extract_From_FID(Mouse, ScanPara); % extract data and scan parameters
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute voronoi density compensation function
if ~isfile([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_DCF])
    Fcn_Compute_Voronoi_DCF(Mouse, ScanPara, ReconPara);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 3D UTE GRASP reconstruction %%%%%%%
for iEcho=1:1:ScanPara.nEchoes
    Fcn_Recon_3D_UTE_GRASP(Mouse, ScanPara, ReconPara, iEcho);
end
