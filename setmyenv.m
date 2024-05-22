function [fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv(flags)
%Sets my environment (paths, etc.)
%   Detailed explanation goes here
if ismac()
    addpath('./cbrewer/');
    if isempty(which('AtlasViewerGUI'))
        mypwd=pwd;cd('/Volumes/rdsprojects/m/montersa-fnirs-connectivity/sw_repos/atlasviewer_smh');setpaths;cd(mypwd);
    end 
    if isempty(which('Homer3'))
        mypwd=pwd;cd('/Volumes/rdsprojects/m/montersa-fnirs-connectivity/sw_repos/Homer3_smh/');setpaths;cd(mypwd);
    end
    if contains(flags.parcel_scheme,'aal')
        fwFolder = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/fw/';
        anatomFolder = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/probe_10MPhotons/anatomical/';
    elseif contains(flags.parcel_scheme,'schaefer')
        fwFolder = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/colinredo/fw/';
        anatomFolder = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/colinredo/anatomical/';
    end
    derivFolder = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/derivatives/rsfc/';
    dataDir = '/Volumes/rdsprojects/m/montersa-fnirs-connectivity/RSFC_BU/Data/DataRSFC_Analysis/';
elseif isunix()
    maxNumCompThreads(32);
    addpath('./cbrewer/');
    if isempty(which('AtlasViewerGUI'))
        mypwd=pwd;cd('../atlasviewer_repo/');setpaths;cd(mypwd);
    end
    if isempty(which('Homer3'))
        %mypwd=pwd;cd('../homer3_repo/');setpaths;cd(mypwd);
        mypwd=pwd;cd('../Homer3_smh_repo/');setpaths;cd(mypwd);
    end
    if contains(flags.parcel_scheme,'aal')
        fwFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/fw/';
        anatomFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/probe_10MPhotons/anatomical/';
    elseif contains(flags.parcel_scheme,'schaefer')
        fwFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/colinredo/fw/';
        anatomFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/colinredo/anatomical/';
    end
    derivFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/derivatives/rsfc/';
    dataDir = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/';
end
if ispc()
    HD_ = 'X';
    usrname=getenv('USERNAME');
    if strcmp(usrname,'samue')
        HD_ = 'C';
    end
    if isempty(which('cbrewer'))
        addpath('./cbrewer/');
    end
    if isempty(which('Homer3'))
        mypwd = pwd;cd([HD_,':\Homer3_smh_repo\']);setpaths;cd(mypwd);
    end
    if isempty(which('AtlasViewerGUI'))
        mypwd = pwd;cd([HD_,':\atlasViewer-repo\']);setpaths;cd(mypwd);
    end
    if contains(flags.parcel_scheme,'aal')
        fwFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\fw\'];
        anatomFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_10MPhotons\anatomical\'];
    elseif contains(flags.parcel_scheme,'schaefer')
        fwFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_schaefer\fw\'];
        anatomFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_schaefer\anatomical\'];
    end
  
    derivFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\derivatives\rsfc\'];
    dataDir = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\'];
end

end