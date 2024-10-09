function [ ] = Fcn_Compute_Voronoi_DCF( Mouse, ScanPara, ReconPara )
% Fcn_Compute_Voronoi_DCF Compute the voronoi based density compensation
% funciotn and save to local directory
%
%   Input: 
%       Mouse: Mouse property structure 
%       ScanPara: Scan parameter structure
%       ReconPaara: Reconstruction parameter structure
% 
%   Output: 
%       density compensation function is saved as a matlab .mat file

load([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_Rawdata]);

traj=reshape(k,[ScanPara.nTrajSamples,ScanPara.nTotalProj,3]);
nNT=floor(ScanPara.nTotalProj/(ReconPara.T/ScanPara.TR));
nSpokes=floor(ReconPara.T/ScanPara.TR);

for ii=1:1:nNT
    %ii 
    x=traj(:,(ii-1)*nSpokes+1:ii*nSpokes,1);
    x=x(:);
    y=traj(:,(ii-1)*nSpokes+1:ii*nSpokes,2);
    y=y(:);
    z=traj(:,(ii-1)*nSpokes+1:ii*nSpokes,3);
    z=z(:);
    X=[x y z];
    
    [V,C]=voronoin(X);
    
    bound=max(X(:));
    V(1,:)=bound;
           
    nHulls=length(C);
    for jj=1:1:nHulls
        temp=C{jj};
        for kk=1:1:length(temp)
            XYZ(kk,:)=V(temp(kk),:);
        end
        % [K,v] = convhulln(XYZ,{'Qt'});
        [K,v] = convhulln(XYZ,{'Qt', 'Pp'});
        volume(jj)=v;
        XYZ=[];
    end
    
    % extropolate the last 6 points
    tempW=reshape(volume, [ScanPara.nTrajSamples, nSpokes]);
    Weight=tempW;
    x=1:6;
    for jj=1:1:nSpokes
        y=tempW(ScanPara.nTrajSamples-6:ScanPara.nTrajSamples-1,jj);
        p=polyfit(x,y',1);
        y_end=p(1)*7+p(2);
        Weight(ScanPara.nTrajSamples,jj)=y_end;
    end

    DCF.raw(ii,:,:)=Weight;
    DCF.frame(ii,:,:)=Weight./max(Weight(:));    
end
temp=DCF.raw;
tempmax=max(temp(:));
DCF.whole=DCF.raw./tempmax;

% save the density compensation function
save([Mouse.PathBase Mouse.PathRecon Mouse.ID '/' Mouse.SavedFile_DCF], 'DCF');
end
