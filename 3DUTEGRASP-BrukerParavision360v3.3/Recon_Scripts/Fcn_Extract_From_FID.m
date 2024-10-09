function [ ] = Fcn_Extract_From_FID(Mouse, ScanPara)
% Fcn_Extract_From_FID Extract fid and traj information from the raw data and save to local direcory
%
%   Input: 
%       Mouse: Mouse property structure 
%       ScanPara: Scan parameter structure
% 
%   Output: 
%       fid and traj varaiables are saved to local directory

pname=[Mouse.PathBase Mouse.PathRecon Mouse.ID '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% read trajectory file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileName=[pname '/' 'traj'];
fi=fopen(FileName, 'rb', 'ieee-le');
[data, cnt]=fread(fi, inf, 'double');
fclose(fi);
temp=reshape(data,3,cnt/3);
temp=temp';
Km=temp(1:ScanPara.nTotalProj*ScanPara.nTrajSamples,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read fid file
FileName=[pname 'rawdata.job0'];
acqpFile='acqp';

strFileFormat=Fcn_Bruker_Method(pname, acqpFile, 'GO_raw_data_format');
if strcmp(strFileFormat, 'GO_32BIT_SGN_INT')
    FileFormat='int32';
    nBit=32;
elseif strcmp(strFileFormat, 'GO_16_BIT_SGN_INT')
    FileFormat='int16';
    nBit=16;
elseif strcmp(strFileFormat, 'GO_32_BIT_FLOAT')
end

nACQsize=Fcn_Bruker_Method(pname, acqpFile, 'ACQ_size');

fi=fopen(FileName, 'rb', 'ieee-le');
[data, cnt]=fread(fi, inf, 'int32');
fclose(fi);

data=reshape(data, 2, cnt/2);
cdata=complex(data(1,:), data(2,:));

rawdata=reshape(cdata, [ScanPara.nTrajSamples, ScanPara.nReceiver, ScanPara.nEchoes, ScanPara.nTotalProj]);
rawdata=permute(rawdata,[1,4,2,3]); % [nTrajSamples nTotalProj nReceiver nEchoes]

k=Km;
    
save([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_Rawdata], 'k', 'rawdata');

end