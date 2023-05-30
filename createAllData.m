% Make sure to import metadata so that 
% the dates and genotypes are text. 
% otherwise you will get weird values.... 

AllData_auditoryTNT = struct;
%% 

% check which OS currently using
if ismac
    prefixdir = ['/Volumes/fileset-mmurthy'];
else
    prefixdir = ['M:'];
end

% define directories for later use
homedir = [prefixdir, filesep, 'Kyle/code/notebooks/auditoryTNT/processExpts'];
resultsdir = [homedir, filesep, 'features'];
dataDir = [prefixdir, filesep, 'Christa/auditoryTNT'];

% import metadata... I import by hand just to make sure datatypes are
% correct


%% 

% Update Info column of allData with metadata info.
for c = 1:length(metadata.exptName)
    AllData_auditoryTNT(c).Info.Sheet = 1;
    AllData_auditoryTNT(c).Info.SheetName = {metadata.femaleGT(c)};
    AllData_auditoryTNT(c).Info.Folder = [dataDir char(metadata.exptName(c))];
    AllData_auditoryTNT(c).Info.ExperimentDate = char(metadata.exptDate(c));
    AllData_auditoryTNT(c).Info.MaleGenotype = metadata.maleGT(c);
    AllData_auditoryTNT(c).Info.MaleDOB = metadata.maleDOB(c);
    if isnan(metadata.copFrame(c))
       AllData_auditoryTNT(c).Info.TimeToCopulation = char('NC');
    else
        AllData_auditoryTNT(c).Info.TimeToCopulation = char(metadata.copFrame(c));
    end
    AllData_auditoryTNT(c).Info.FemaleGenotype = metadata.femaleGT(c);
    AllData_auditoryTNT(c).Info.FemaleDOB = metadata.femaleDOB(c);
    AllData_auditoryTNT(c).Info.MaleAgeAtExperiment = metadata.maleAge(c);
    AllData_auditoryTNT(c).Info.FemaleAgeAtExperiment = metadata.femaleAge(c);
    AllData_auditoryTNT(c).Info.IsReadyForAnalysis = 1;

end

%% 

% find session files.
fName = '_features.h5';
inDir = strcat(resultsdir, filesep);

fprintf('finding files in %s...\n', inDir);
dirList = regexpdir(inDir, fName); 

numSessions = length(dirList);
for sessNum = 1:numSessions 
    fprintf('Loading video data %i... \n', sessNum);
    
%     flopped sexes switching and testing difference
    ftrfile = dirList{sessNum};
    fFV = h5read(ftrfile, '/fFV').';
    mFV = h5read(ftrfile, '/mFV').';
    fLV = h5read(ftrfile, '/fLV').';
    mLV = h5read(ftrfile, '/mLV').';
    fRV = h5read(ftrfile, '/fRS').';
    mRV = h5read(ftrfile, '/mRS').';
    mfDist = h5read(ftrfile, '/mfDist').';
    
    AllData_auditoryTNT(sessNum).Tracking.mFV = mFV;
    AllData_auditoryTNT(sessNum).Tracking.fFV = fFV;
    AllData_auditoryTNT(sessNum).Tracking.mLV = mLV;
    AllData_auditoryTNT(sessNum).Tracking.fLV = fLV;
    AllData_auditoryTNT(sessNum).Tracking.mRV = mRV;
    AllData_auditoryTNT(sessNum).Tracking.fRV = fRV;
    AllData_auditoryTNT(sessNum).Tracking.mfDist = mfDist;

end

%% 
for a = 1:length(metadata.exptName)
    fprintf('Loading song data in %i... \n', a)
    
    audioname = 'daq_segmentation_new.mat';
    songPath = strcat(dataDir, filesep, char(metadata.exptName(a)), filesep, audioname);
    load(songPath, 'bInf', 'pInf')
    MASK = bInf.Mask;
    All_Pulses = pInf.wc;
    stEn = bInf.stEn;
    Type = bInf.Type;
    
    AllData_auditoryTNT(a).Audio.MASK = MASK;
    AllData_auditoryTNT(a).Audio.All_Pulses = All_Pulses;
    AllData_auditoryTNT(a).Audio.stEn = stEn;
    AllData_auditoryTNT(a).Audio.Type = Type;

end
%% 

for sessNum = 1:numSessions
    fprintf('Loading sync data in %i... \n', sessNum);
    ftrfile = dirList{sessNum};
    fprintf('Session... %i \n', sessNum);
    vSampleNumber_at_Frame = h5read(ftrfile, '/sample_at_frame');
    AllData_auditoryTNT(sessNum).Sync.vSampleNumber_at_Frame = vSampleNumber_at_Frame;

end
%% 

savePath = strcat(homedir, filesep, 'auditoryTNT_allData_2021_11_24.mat');
save(savePath, 'AllData_auditoryTNT', '-v7.3');