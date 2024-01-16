function [fwFolder,anatomFolder,derivFolder,dataDir] = setmyenv()
%Sets my environment (paths, etc.)
%   Detailed explanation goes here
if isunix()
    maxNumCompThreads(32);
    addpath('./cbrewer/');
    if isempty(which('AtlasViewerGUI'))
        mypwd=pwd;cd('../atlasviewer_repo/');setpaths;cd(mypwd);
    end
    if isempty(which('Homer3'))
        mypwd=pwd;cd('../homer3_repo/');setpaths;cd(mypwd);
    end
    fwFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/fw/';
    anatomFolder = '/projectnb/nphfnirs/s/DATA_BU/2022/Rest_Movie_WorkingMemory/DataRSFC_Analysis/probe_10MPhotons/anatomical/';
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
    fwFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\fw\'];
    anatomFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\probe_10MPhotons\anatomical\'];
    derivFolder = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\derivatives\rsfc\'];
    dataDir = ['C:\Users\',usrname,'\OneDrive - Boston University\RS_MovieWatching\Rest_Movie_WorkingMemory\'];
end

end