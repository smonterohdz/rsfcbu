function [HbO_brain,HbR_brain] = ImageReconstruction_FC(snirfObj,dodObj,dcObj,fwFolder,flags)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Load Brain Adot
load([fwFolder 'Adot.mat'],'Adot');%
index_select = log10(sum(Adot(:,:,1),1))>=-2;
n_vox = sum(index_select);
Adot_select = Adot(:,index_select,:);
Adot = Adot_select;


% Load Scalp Adot
load([fwFolder 'Adot_scalp.mat'],'Adot_scalp');%
index_select_scalp = log10(sum(Adot_scalp(:,:,1),1))>=-2;
Adot_select = Adot_scalp(:,index_select_scalp,:);
Adot_scalp = Adot_select;


% Get long separation channels and mark inacive all sd pairs that have short sep optodes
SD = snirfObj.GetSDG();
SrcPos = SD.SrcPos();
DetPos = SD.DetPos();

ml_dod  = dodObj.GetMeasurementList('matrix');
if strcmp(flags.gsr,'channel') %if departs from Conc
    SD_ = snirfObj.Get_SD;
    [dc,time_,ml_] = dcObj.GetDataTimeSeries('reshape');    
    dod = hmrConc2OD( dc, SD_, [1 1] );
else
    dod = dodObj.dataTimeSeries;
end

for iML = 1:length(ml_dod)
    rho = sum( (SrcPos(ml_dod(iML,1),:) - DetPos(ml_dod(iML,2),:)) .^ 2) ^ 0.5;
    if rho < flags.rhoSD_ssThresh        
        k1 = ml_dod(:,1)==ml_dod(iML,1) & ml_dod(:,2)==ml_dod(iML,2);        
        ml_dod(k1,3) = 0;        
    end
end
k = ml_dod(:,4)==1;
activeChLst_SDpairs_singWL = find(ml_dod(k,3) == 1);
activeChLst_SDpairs_bothWL = find(ml_dod(:,3) == 1);

dod = dod(:,activeChLst_SDpairs_bothWL);

Adot = Adot(activeChLst_SDpairs_singWL,:,:);
Adot_scalp = Adot_scalp(activeChLst_SDpairs_singWL,:,:);

% Preparing AdotBrain,AdotScalp
%snirf.probe.wavelengths = [760,   850]
E = GetExtinctions([SD.Lambda(1) SD.Lambda(2)]);
E = E/10; %convert from /cm to /mm  E raws: wavelength, columns 1:HbO, 2:HbR

if strcmp(flags.imagerecon,'brain')
% Amatrix = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
%         squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)];
%function [HbO, HbR, err] = hmrImageReconConc(dodavgimg, dodresid, alpha, Adot)
Amatrix = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
    squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)];
alpha = 0.1;
[HbO, HbR, err] = hmrImageReconConc(dod', [], alpha, Amatrix);
if isa(HbO,'double')
    HbO = single(HbO);
end
if isa(HbR,'double')
    HbR = single(HbR);
end
end

if strcmp(flags.imagerecon,'brain+scalp')
    Amatrix = [[squeeze(Adot(:,:,1)) squeeze(Adot_scalp(:,:,1))]*E(1,1) [squeeze(Adot(:,:,1)) squeeze(Adot_scalp(:,:,1))]*E(1,2);
        [squeeze(Adot(:,:,2)) squeeze(Adot_scalp(:,:,2))]*E(2,1) [squeeze(Adot(:,:,2)) squeeze(Adot_scalp(:,:,2))]*E(2,2)];


    iA = Tikhonov_invert_Amat(Amatrix, 0.01, 0.1);

    % remove bad channels
    %Bad channels pruning
    %badChannels = Proc_data.BadChannels;
    %badChannels = ~logical(output.misc.mlActAuto{1}(:,3));
    %badChannels = [];
    iA_pruned = iA;
    %iA_pruned(:,badChannels) = [];
    %dod(:,badChannels) = [];

    Hb = iA_pruned*dod';
    HbO = Hb(1:length(iA_pruned)/2,:);
    HbR = Hb((length(iA_pruned)/2)+1:end,:);
end
HbO_brain = [HbO(1:n_vox,:)]';
HbR_brain = [HbR(1:n_vox,:)]';

% Global signal regression in image space?
if strcmp(flags.gsr,'image')
    HbO_brain = GlobalRegression(HbO_brain);
    HbR_brain = GlobalRegression(HbR_brain);
end
end


function Hb_brain = GlobalRegression(Hb_brain)
A = mean(Hb_brain,1);
A_inv = A'/(A*A');
for i = 1:size(Hb_brain,1)
    y = Hb_brain(i,:);
    b = y*A_inv;
    Hb_brain(i,:) = y - A*b;
end
end