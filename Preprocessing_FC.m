function [snirfObj,dc,dod,dcNP,dodNP] = Preprocessing_FC(fullPathSnirf,flags)
%Preprocessing_FC Performs fnirs preprocessing on channel space
%   Detailed explanation goes here


% load snirf file
snirfObj = SnirfClass(fullPathSnirf);


% Prune Channels
%mlActAuto = hmrR_PruneChannels(data, probe, mlActMan, tInc, dRange, SNRthresh, SDrange)
mlActAuto = [];
% mlActMan = [];
% tInc = [];
% dRange = [1e-4 1e1];
% SNRthresh = 5;
% SDrange = [0 45];
% mlActAuto = hmrR_PruneChannels(snirfObj.data, snirfObj.probe, mlActMan, tInc, dRange, SNRthresh, SDrange);
% % if mlActAuto=0 for either wavelength of a channel, then set all
% % wavelengths to 0
% lst1=find(mlActAuto{1}(:,4)==1);
% lst2=find(mlActAuto{1}(:,4)==2);
% mlActAuto{1}(lst1,3) = min(mlActAuto{1}(lst1,3),mlActAuto{1}(lst2,3));


% Convert to dOD
dod = hmrR_Intensity2OD( snirfObj.data );
dodNP = dod;
if strcmp(flags.macorrect,'spline')
    % Find Motion Artifacts
    % tIncAuto = hmrR_MotionArtifact(data, probe, mlActMan, mlActAuto, tIncMan, tMotion, tMask, STDEVthresh, AMPthresh)
    tMotion = 0.5;
    tMask = 0.5;
    STDEVthresh = 40;
    AMPthresh = 100;
    %tIncAuto = hmrR_MotionArtifact(dod, snirfObj.probe, [], mlActAuto, [], tMotion, tMask, STDEVthresh, AMPthresh);
    % By channel
    % [tInc, tIncCh] = hmrR_MotionArtifactByChannel(data, probe, mlActMan, mlActAuto, tIncMan, tMotion, tMask, std_thresh, amp_thresh)
    [tInc, tIncCh] = hmrR_MotionArtifactByChannel(dod, snirfObj.probe, [], mlActAuto, [], tMotion, tMask, STDEVthresh, AMPthresh);

    % motion artifacts correction
    p =0.99;
    turnon = 1;
    FrameSize_sec = 10;
    %USAGE data_dod = hmrR_MotionCorrectSpline(data_dod, mlAct, tIncCh, p, turnon)
    dod = hmrR_MotionCorrectSpline(dod, mlActAuto, tIncCh, p, turnon);
    %USAGE function data_d = hmrR_MotionCorrectSplineSG(data_d, mlActAuto, p, FrameSize_sec, turnon)
    %dod = hmrR_MotionCorrectSplineSG(dod, mlActAuto, p, FrameSize_sec, turnon);
end

if strcmp(flags.bpfilt,'channel')
    % Band Pass Filter
    dod = hmrR_BandpassFilt( dod, 0.009, 0.080);
end

% Convert to Conc
dc = hmrR_OD2Conc( dod, snirfObj.probe, [1 1]);
dcNP = hmrR_OD2Conc(dodNP,snirfObj.probe,[1 1]);

% Global signal regression in image space?
if strcmp(flags.gsr,'channel')
    d = dc.GetDataTimeSeries('reshape');
    HbO = GlobalRegression(squeeze(d(:,1,:)));
    HbR = GlobalRegression(squeeze(d(:,2,:)));
    d(:,1,:) = HbO;
    d(:,2,:) = HbR;
    d(:,3,:) = HbO+HbR;
    dc.SetDataTimeSeries(d);
end

end
