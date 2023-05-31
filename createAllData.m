% Make sure to import metadata so that 
% the dates and genotypes are text. 
% otherwise you will get weird values.... 

allData = struct;
%% 

% check which OS currently using
prefixdir = '/run/user/1000/gvfs/smb-share:server=cup.pni.princeton.edu,share=murthy';

% define directories for later use
homedir = [prefixdir, filesep, 'Kyle/flies/pc2_tnt'];


% import metadata... I import by hand just to make sure datatypes are
% correct


%% 

% Update Info column of allData with metadata info.
for c = 1:length(metadata.exptName)
    allData(c).Info.Sheet = 1;
    allData(c).Info.SheetName = {metadata.femaleGT(c)};
    allData(c).Info.Folder = strcat(homedir, filesep, metadata.expt(c), filesep, metadata.rig(c));
    allData(c).Info.ExperimentDate = char(metadata.exptDate(c));
    allData(c).Info.MaleGenotype = metadata.maleGT(c);
    allData(c).Info.MaleDOB = metadata.maleDOB(c);
    if isnan(metadata.copFrame(c))
       allData(c).Info.TimeToCopulation = char('NC');
    else
        allData(c).Info.TimeToCopulation = metadata.copFrame(c)/ 60 ;
    end
    allData(c).Info.FemaleGenotype = metadata.femaleGT(c);
    allData(c).Info.FemaleDOB = metadata.femaleDOB(c);
    allData(c).Info.MaleAgeAtExperiment = metadata.maleAge(c);
    allData(c).Info.FemaleAgeAtExperiment = metadata.femaleAge(c);
    allData(c).Info.IsReadyForAnalysis = 1;
    

    ftrfile = metadata.exptpath(c);

    fFV = h5read(ftrfile, '/fFV').';
    mFV = h5read(ftrfile, '/mFV').';
    fLV = h5read(ftrfile, '/fLV').';
    mLV = h5read(ftrfile, '/mLV').';
    fRV = h5read(ftrfile, '/fRS').';
    mRV = h5read(ftrfile, '/mRS').';
    mfDist = h5read(ftrfile, '/mfDist').';
    
    allData(c).Tracking.mFV = mFV;
    allData(c).Tracking.fFV = fFV;
    allData(c).Tracking.mLV = mLV;
    allData(c).Tracking.fLV = fLV;
    allData(c).Tracking.mRV = mRV;
    allData(c).Tracking.fRV = fRV;
    allData(c).Tracking.mfDist = mfDist;

    songpath = metadata.songpath(c);
    load(songpath, 'bInf', 'pInf');
    MASK = bInf.Mask;
    All_Pulses = pInf.wc;
    stEn = bInf.stEn;
    Type = bInf.Type;

    allData(c).Audio.MASK = MASK;
    allData(c).Audio.All_Pulses = All_Pulses;
    allData(c).Audio.stEn = stEn;
    allData(c).Audio.Type = Type;

    datapath = metadata.datapath(c);
    
    vSampleNumber_at_Frame = h5read(datapath, '/vSampleNumber_at_Frame').';
    allData(c).Sync.vSampleNumber_at_Frame = vSampleNumber_at_Frame;
end

%% 

save('/home/kyle/Desktop/allData_2.mat', 'allData', '-v7.3');
