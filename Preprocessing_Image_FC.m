function [outputArg1,outputArg2] = Preprocessing_Image_FC(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
% Band pass filtering in image space
if strcmp(flags.bpfilt,'image')
    %[y2] = image_BandpassFilt(y,hpf,lpf,fs)
    fs = mean(1./diff(snirfObjr1.data.time));
    HbO_brain_chunkr1 = image_BandpassFilt(HbO_brain_chunkr1, 0.009, 0.080,fs);
    HbR_brain_chunkr1 = image_BandpassFilt(HbR_brain_chunkr1, 0.009, 0.080,fs);
    HbO_brain_chunkr2 = image_BandpassFilt(HbO_brain_chunkr2, 0.009, 0.080,fs);
    HbR_brain_chunkr2 = image_BandpassFilt(HbR_brain_chunkr2, 0.009, 0.080,fs);
end
% Global signal regression in image space?
if strcmp(flags.gsr,'image') && ~strcmp(flags.bpfilt,'none')
    HbO_brain_chunkr1 = GlobalRegression(HbO_brain_chunkr1);
    HbR_brain_chunkr1 = GlobalRegression(HbR_brain_chunkr1);
    HbO_brain_chunkr2 = GlobalRegression(HbO_brain_chunkr2);
    HbR_brain_chunkr2 = GlobalRegression(HbR_brain_chunkr2);
end

end