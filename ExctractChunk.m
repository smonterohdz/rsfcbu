function [HbXBrain_chunk] = ExctractChunk(HbXBrain,snirfObj,twindow)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(twindow,'init_sec') || isempty(twindow.init_sec) 
    twindow.init_sec = 90;
end
if ~isfield(twindow,'dur_sec') || isempty(twindow.dur_sec)
    twindow.dur_sec = 900;
end
if ~isfield(twindow,'offset_sec') || isempty(twindow.offset_sec)||twindow.offset_sec<1
    twindow.offset_sec = 0;
end
if twindow.init_sec == -1
    if ~(isempty(snirfObj.stim))
        twindow.init_sec = snirfObj.stim.data(1);
        twindow.dur_sec = snirfObj.stim.data(2);
    else
        twindow.init_sec = 60;
        twindow.dur_sec = 900;
        warning('No stimuli information.\n');
    end
end
twindow.init_sec = twindow.init_sec + twindow.offset_sec;
twindow.dur_sec = twindow.dur_sec - twindow.offset_sec;
twindow.end_sec = twindow.init_sec + 2;
fs = 1/mean(diff(snirfObj.data.time));

tinit = floor(twindow.init_sec*fs);% floor(twindow.shift_sec*fs):floor(twindow.end_sec*fs)
tend = tinit + floor(twindow.dur_sec*fs);
if size(HbXBrain,1)<tend
    tend = size(HbXBrain,1) - ceil(fs);
    twindow.dur_sec = floor((tend-tinit)/fs);    
end
fprintf('Using onset:%.1fs and duration:%.1fs\n',twindow.init_sec,twindow.dur_sec)
HbXBrain_chunk = HbXBrain(tinit:tend,:);
end