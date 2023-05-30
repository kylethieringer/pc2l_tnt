%The function Analyze_Behavior is comparing male and female behavior (song,
%response to song) for a number of groups (2 groups or more).
%The function can run on a windows or on unix

%function Analyze_Behavior_Jans(varargin)

%parameters
DataFile = 'allData.mat';%here the structure Alldata is saved
SaveTo = 'output';%Here some outputs from this function will be saved
param.nTakeOnly_NC = 0;% if param.nTakeOnly_NC = 1 - the function will consider only experiments with NoCopulation
IsReload = 0;
DataFolder = '/run/user/1000/gvfs/smb-share:server=cup.pni.princeton.edu,share=murthy/Kyle/code/pc2_tnt/';
DataFile = [DataFolder DataFile];
Sep = filesep();
if ~exist('allData.mat','var') || IsReload == 1, load(DataFile),end

Groups = {'pc2_tnt','pc2_control'};
vGroupsToAnalyze = 1:2;%Which tabs in the excell table (=which groups) to use for analysis

Folders = Groups(vGroupsToAnalyze);
N_Groups = length(Folders);

cd(DataFolder)

param.nMin_Frames = 120;
Param.CurrTime = datestr(datenum(clock));%Time at start of function
param.nfps = 60; %Frames per second. Used only to translate time to copulatin from seconds (in the table..) to frames
param.nSamplesPerSeconds = 10000;%In the audio channels
%parameters for female slowing in bins
param.BinWidth_FemaleSongResponse = 60;%In seconds
param.PercentOverlap = 50;%  % overlap between bins
param.IgnoreSpeedBeforeSong = 0.3;%sec - when looking for the female speed - look outside song bouts and param.IgnoreSpeedAroundSong seconds around them
param.IgnoreSpeedAfterSong = 0.3;
%parameters for the statistics of response after start of bout
param.IgnoreSeconds = 0.2;%
param.TransientSeconds = 0.6;
param.SteadyStateSeconds = 0.6;

%Count how many complete files
vSheets = zeros(1,size(allData,2));
for nTrial = 1:size(allData,2)
    vSheets(nTrial) =  allData(nTrial).Info.Sheet;
end
vSheets = unique(vSheets);

sCompleteFolders = struct;
for ii=1:max(vSheets), sCompleteFolders(ii).N_Folders = 0;end
for ii=1:max(vSheets), sCompleteFolders(ii).N_CompleteFolders = 0;end
for nTrial = 1:size(allData,2)
    nSheet = allData(nTrial).Info.Sheet;
    sCompleteFolders(nSheet).Name = allData(nTrial).Info.SheetName;
    sCompleteFolders(nSheet).N_Folders = sCompleteFolders(nSheet).N_Folders + 1;
    if allData(nTrial).Info.IsReadyForAnalysis
        sCompleteFolders(nSheet).N_CompleteFolders = sCompleteFolders(nSheet).N_CompleteFolders +1;
    end
end

for nSheet = vSheets
    disp([sCompleteFolders(nSheet).Name{1},' : ',num2str(sCompleteFolders(nSheet).N_CompleteFolders),'/',...
        num2str(sCompleteFolders(nSheet).N_Folders),' complete folders'])
end

% Folders = cell(1,length(vGroupsToAnalyze));
% for nGroup = 1:length(vGroupsToAnalyze)
%     Folder_to_Analyze = Groups{vGroupsToAnalyze};
%     for nTab = 1:length(sCompleteFolders)
%         Folder_in_Tab = sCompleteFolders(nTab).Name;
%         if strcmp(Folder_to_Analyze,Folder_in_Tab)
%             Folders{nGroup} = sCompleteFolders(vGroupsToAnalyze(nGroup)).Name{1};
%             break
%         end
%     end
% end



%Initialize Female speed Vs male song
MaleSong_FliesKinematics = struct('MaleSong_FemaleSpeed_Binned',[]);
for nGroup = 1:size(Folders,2)
    MaleSong_FliesKinematics(nGroup) = struct('MaleSong_FemaleSpeed_Binned',[]);
end


%% Run over the trial, align female speed to male pulse song

for ii = 1:size(allData,2)
    Folder = allData(ii).Info.FemaleGenotype;
    
    nGroup=[];
    for jj = 1:size(Folders,2)
        
        Folder_Strain = Folder;
        if strcmp(Folder_Strain,Folders{jj});nGroup = jj;break;end%Find to which group this trial belongs
    end
    
    
    %folders to skip - no tracking, too short. Trials with copulation when param.nTakeOnly_NC = 1
    if isempty(nGroup) || ~allData(ii).Info.IsReadyForAnalysis,continue,end
    %load Synchronization file
    
    vSampleNumber_at_Frame = int64(allData(ii).Sync.vSampleNumber_at_Frame); %synchronizing frames and samples
    vfaas = allData(ii).Sync.vSampleNumber_at_Frame;
    Is_Copulation = ~strcmp(allData(ii).Info.TimeToCopulation,'NC');
    if param.nTakeOnly_NC && Is_Copulation,continue,end
    
    
    %Start copulation
    if ~isfield(allData(ii).Info,'ReachStartFrame') || isempty(allData(ii).Info.ReachStartFrame) ||...
            isnan(allData(ii).Info.ReachStartFrame)
    START_COURTSHIP = 1;
    else START_COURTSHIP =  allData(ii).Info.ReachStartFrame;
    end
    
    
    %ignore the frames after copulation, if there was copulation
    if Is_Copulation%copulation happened
        nFrame_StartCopulation = START_COURTSHIP + allData(ii).Info.TimeToCopulation*param.nfps - 1;
    else%no copulation
        nFrame_StartCopulation = inf;
    end
    if nFrame_StartCopulation<param.nMin_Frames*param.nfps,continue,end%movie is too short
    %end of - folders to skip
    
    disp(['Trial ',num2str(ii),', Folder ',Folder,' is being analyzed'])
    
    cd(DataFolder);
    
    %Tracking data from Alldata
    mFV = allData(ii).Tracking.mFV'; fFV = allData(ii).Tracking.fFV';
    mLV = allData(ii).Tracking.mLV'; fLV = allData(ii).Tracking.fLV';
    mRV = allData(ii).Tracking.mRV'; fRV = allData(ii).Tracking.fRV';
    fDist = allData(ii).Tracking.mfDist';
    
    %Song info
    All_Pulses = allData(ii).Audio.All_Pulses;
    if ~~isempty(All_Pulses), continue, end
    
    %BINF = allData(ii).Audio.BINF; %BINF.stEn is N*2 for N chunks, where a chunk j has one song type (NoSong/Pulse/Sine)
    %from BINF.stEn(j,1) to BINF.stEn(j,2) and BINF.Mask is 0/1/2 for each sample
    
    %Last frame to consider - a minute before copulation (if copulation), tracking and synchronization OK.
    nEnd_Tracking = find(mFV~=0); nEnd_Tracking =  nEnd_Tracking(end);%No tracking after nEnd
    nFirst_NonZero = find(vSampleNumber_at_Frame>0,1);
    nEnd_Synch = find(vSampleNumber_at_Frame(nFirst_NonZero:end)==0,1)+nFirst_NonZero-2;if isempty(nEnd_Synch),nEnd_Synch = length(vSampleNumber_at_Frame);end
%     nEndGenderID = find(allData(ii).IsMaleFly(2:end) == -1 & allData(ii).IsMaleFly(1:end-1) ~= -1);

    %Ignore the frames after last pulse
    nFrame_LastPulse = find(vSampleNumber_at_Frame>All_Pulses(end),1);
    if isempty(nFrame_LastPulse), nFrame_LastPulse = inf;end
    
    %LAST FRAME TO CONSIDER IN THIS ANALYSIS.
    %Note - 1 minute before copulation is ignored!
    nEnd = min([nEnd_Tracking nEnd_Synch nFrame_StartCopulation-param.nfps nFrame_LastPulse]); %nEndGenderID]);
    
    vSampleNumber_at_Frame = vSampleNumber_at_Frame(1:nEnd);
    %extract male and female behavior
    fFV = medfilt1(fFV(1:nEnd))';mFV = medfilt1(mFV(1:nEnd))';
    fLV = medfilt1(fLV(1:nEnd))';mLV = medfilt1(mLV(1:nEnd))';
    fRV = medfilt1(fRV(1:nEnd))';mRV = medfilt1(mRV(1:nEnd))';
    dis = medfilt1(fDist(1:nEnd));
    
    
    %Filter
    % fFV = medfilt1(fFV); fLV = medfilt1(fLV); fRV = medfilt1(fRV);
    % mFV = medfilt1(mFV); mLV = medfilt1(mLV); mRV = medfilt1(mRV);
    
    %Absolute velocity
    fAbsoluteVelocity = sqrt(fFV.^2+fLV.^2);
    mAbsoluteVelocity = sqrt(mFV.^2+mLV.^2);
    
    
    %Limit the song - last sample at the same time as last frame for the tracking data
    %     Mask = Alldata(ii).MASK(1:min([floor(vSampleNumber_at_Frame(nEnd)) length(Alldata(ii).MASK)]));
    %     stEn = Alldata(ii).stEn(Alldata(ii).stEn(:,2)<floor(nEnd/param.nfps*param.nSamplesPerSeconds),:);
    
    %Define the moving time window for the analysis
    nBinWidth = param.nfps*param.BinWidth_FemaleSongResponse;
    vBins_femaleresponse = START_COURTSHIP:nBinWidth*(param.PercentOverlap/100):length(fAbsoluteVelocity);
    vBins_femaleresponse=vBins_femaleresponse(2:end-1);
    
    %% Build the matrix for calculating the correlation between 11 binned song
    %features and female speed (based on Clemens 2015 Neuron paper)
    
    %Mark song and around song time - so that female speed will be
    %calculated only when no-song, and away with some distance from songs.
    vAnySongWrapped = zeros(size(allData(ii).Audio.MASK));
    for nSongBin = 1:size(allData(ii).Audio.stEn,1)
        nStartSong_Bin = max(1,allData(ii).Audio.stEn(nSongBin,1)-param.IgnoreSpeedBeforeSong*param.nSamplesPerSeconds);
        nEndSong_Bin = min(length(allData(ii).Audio.MASK),allData(ii).Audio.stEn(nSongBin,2)+param.IgnoreSpeedAfterSong*param.nSamplesPerSeconds);
        vAnySongWrapped(nStartSong_Bin:nEndSong_Bin) = 1;%contains the song periods + some margin around them
    end
    vAnySongWrapped(1:param.nSamplesPerSeconds) = 1;%First minute is ignored
    %Also - last minute before copulation (nEnd =... see above)
    
    
    
    %1-11 song features, based on Clemans et al 2015 (same order as in the x-axis in figure 1)
    %12-15 female speeds (FWD, lateral, rotational, absolute), 16-19 males speeds (same as female), 20 male-female distance
    %21-23 are the index of in AllData and the sample at the beginning/end of each bin
    mMaleSong_FliesKinematics_Binned = zeros(length(vBins_femaleresponse),23);
    for nBin_Index = 1:length(vBins_femaleresponse)
        %Bin start/end - look at male song and female response in this bin
        nStartBinFrame = vBins_femaleresponse(nBin_Index)-nBinWidth*(param.PercentOverlap/100);
        nEndBinFrame = nStartBinFrame + nBinWidth;
        if nStartBinFrame<1 || nEndBinFrame>length(fAbsoluteVelocity), continue,end
        nStartBinSample = floor(vSampleNumber_at_Frame(nStartBinFrame));
        nEndBinSample = floor(vSampleNumber_at_Frame(nEndBinFrame));
        if nStartBinSample<1 || nEndBinSample>length(allData(ii).Audio.MASK); continue;end
        
        %Male song - amount
        Songs = allData(ii).Audio.MASK(nStartBinSample:nEndBinSample);
        mMaleSong_FliesKinematics_Binned(nBin_Index,1) = length(find(Songs==1))/length(Songs);%Pulse amount
        mMaleSong_FliesKinematics_Binned(nBin_Index,2) = length(find(Songs==2))/length(Songs);%Sine amount
        mMaleSong_FliesKinematics_Binned(nBin_Index,3) = length(find(Songs>0))/length(Songs);%Song amount
        
        %for the number (onsets) and duration, only complete bouts are used -
        %Take only full bouts for the analysis of bout duration
        Transitions = find(diff(Songs)~=0);
        if length(Transitions)>=2
            nStartSong = Transitions(1) + 1;
            nEndSong = Transitions(end);
            nLastSong = Songs(nStartSong-1);
            Songs = Songs(nStartSong:nEndSong);Songs_shifted = [nLastSong Songs(1:end-1)];
            bFullBoutsExist = 1;
        else
            bFullBoutsExist = 0;
        end
        
        
        %Number of onsets (only full bouts, see above)
        if bFullBoutsExist
            mMaleSong_FliesKinematics_Binned(nBin_Index,4) = length(find(Songs==1 & Songs_shifted~=1));%Pulse number (=number of pulse bouts)
            mMaleSong_FliesKinematics_Binned(nBin_Index,5) = length(find(Songs==2 & Songs_shifted~=2));%Sine number (=number of sine bouts)
            mMaleSong_FliesKinematics_Binned(nBin_Index,6) = length(find(Songs~=0 & Songs_shifted==0));%Song number (=number of song bouts)
        else
            mMaleSong_FliesKinematics_Binned(nBin_Index,4:6) = -1;
        end
        
        %IPI
        PulsesInBin = All_Pulses(All_Pulses>nStartBinSample & All_Pulses<=nEndBinSample);
        if length(PulsesInBin)>=2
            IPI = diff(PulsesInBin); IPI = IPI(IPI<=2000);%IPI>200ms - new bout
            mMaleSong_FliesKinematics_Binned(nBin_Index,7) = mean(IPI);%mean IPI
        else mMaleSong_FliesKinematics_Binned(nBin_Index,7) = -1;
        end
        
        %Mean bout duration (only full bouts, see above)
        if bFullBoutsExist
            if ~(mMaleSong_FliesKinematics_Binned(nBin_Index,4)==0)%Full pulse bouts exist, Pulse (mean) duration:
                mMaleSong_FliesKinematics_Binned(nBin_Index,8) = length(find(Songs==1))/mMaleSong_FliesKinematics_Binned(nBin_Index,4);
            else mMaleSong_FliesKinematics_Binned(nBin_Index,8) = -1;
            end
            if ~(mMaleSong_FliesKinematics_Binned(nBin_Index,5)==0)%Full sine bouts exist, Sine (mean) duration:
                mMaleSong_FliesKinematics_Binned(nBin_Index,9) = length(find(Songs==2))/mMaleSong_FliesKinematics_Binned(nBin_Index,5);
            else mMaleSong_FliesKinematics_Binned(nBin_Index,9) = -1;
            end
            if ~(mMaleSong_FliesKinematics_Binned(nBin_Index,6)==0)%Full song bouts exist, Song (mean) duration:
                mMaleSong_FliesKinematics_Binned(nBin_Index,10) = length(find(Songs>0))/mMaleSong_FliesKinematics_Binned(nBin_Index,6);
            else mMaleSong_FliesKinematics_Binned(nBin_Index,10) = -1;
            end
        else mMaleSong_FliesKinematics_Binned(nBin_Index,8:10) = -1;
        end
        
        %Pauses
        Songs = allData(ii).Audio.MASK(nStartBinSample:nEndBinSample);
        Songs_shifted = [0 Songs(1:end-1)];
        nStartPauses = find(Songs == 0 & Songs_shifted>0);
        nEndPauses = find(Songs > 0 & Songs_shifted==0);
        if isempty(nStartPauses) || isempty(nEndPauses)
            mMaleSong_FliesKinematics_Binned(nBin_Index,11) = -1;
        else
            nStartPauses = nStartPauses(1); nEndPauses = nEndPauses(end);
            Songs = Songs(nStartPauses:nEndPauses); Songs_shifted = [0 Songs(1:end-1)];
            NumberOfFullPauses = length(find(Songs == 0 & Songs_shifted > 0))+1;
            mMaleSong_FliesKinematics_Binned(nBin_Index,11) = length(Songs(Songs == 0))/NumberOfFullPauses;
        end
        
        
        
        %Female speed
        vFramesInBin = nStartBinFrame:nEndBinFrame;
        vSamplesInBin =round(vSampleNumber_at_Frame(vFramesInBin));%The sample number at the frames in the bin
        vNoSong_Indexes = find(vAnySongWrapped(int32(vSamplesInBin))==0);%Which frames are song-free
        %vNoSong_Indexes is the indexes withis the bin that are far enough
        %from any song
        vFemaleSpeed = fFV(nStartBinFrame:nEndBinFrame);vFemaleSpeed = vFemaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,12) = mean(vFemaleSpeed);%mean female FWD velocity
        vFemaleSpeed = fLV(nStartBinFrame:nEndBinFrame);vFemaleSpeed = vFemaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,13) = mean(vFemaleSpeed);%mean female lateral velocity
        vFemaleSpeed = fRV(nStartBinFrame:nEndBinFrame);vFemaleSpeed = vFemaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,14) = mean(vFemaleSpeed);%mean rotational velocity
        vFemaleSpeed = fAbsoluteVelocity(nStartBinFrame:nEndBinFrame);vFemaleSpeed = vFemaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,15) = mean(vFemaleSpeed);%mean female absolute velocity
        
        %Male speed
        vMaleSpeed = mFV(nStartBinFrame:nEndBinFrame);vMaleSpeed = vMaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,16) = mean(vMaleSpeed);%mean female FWD velocity
        vMaleSpeed = mLV(nStartBinFrame:nEndBinFrame);vMaleSpeed = vMaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,17) = mean(vMaleSpeed);%mean female lateral velocity
        vMaleSpeed = mRV(nStartBinFrame:nEndBinFrame);vMaleSpeed = vMaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,18) = mean(vMaleSpeed);%mean rotational velocity
        vMaleSpeed = mAbsoluteVelocity(nStartBinFrame:nEndBinFrame);vMaleSpeed = vMaleSpeed(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,19) = mean(vMaleSpeed);%mean female absolute velocity
        
        %Male-Female distance
        vMaleFemaleDis = dis(nStartBinFrame:nEndBinFrame);vMaleFemaleDis = vMaleFemaleDis(vNoSong_Indexes);
        mMaleSong_FliesKinematics_Binned(nBin_Index,20) = mean(vMaleFemaleDis);%mean female absolute velocity
        
        mMaleSong_FliesKinematics_Binned(nBin_Index,21) = ii;%Index of the trial in AllData
        mMaleSong_FliesKinematics_Binned(nBin_Index,22) = nStartBinFrame;%Frame number in the beginning of the bin
        mMaleSong_FliesKinematics_Binned(nBin_Index,23) = nEndBinFrame;%Frame number in the end of the bin
    end
    
    %Don't take bins with <5% or >95% song - same criterion is used in Clemans et al.
    mMaleSong_FliesKinematics_Binned = mMaleSong_FliesKinematics_Binned(mMaleSong_FliesKinematics_Binned(:,3)>0.05 &...
        mMaleSong_FliesKinematics_Binned(:,3)<0.95,:);
    mMaleSong_FliesKinematics_Binned = mMaleSong_FliesKinematics_Binned(mMaleSong_FliesKinematics_Binned(:,21)>0,:);%Just in case there are extra lines with zeros
    
    %Z-Score every fly individualy [not doing, just an optioned]
    %mMaleSong_FliesKinematics_Binned(:,1:20) = zscore(mMaleSong_FliesKinematics_Binned(:,1:20));
    
    MaleSong_FliesKinematics(nGroup).MaleSong_FemaleSpeed_Binned =...
        [MaleSong_FliesKinematics(nGroup).MaleSong_FemaleSpeed_Binned;mMaleSong_FliesKinematics_Binned];
    
    
    %Looking in a short window after the beginning of Pulse/sine/song bout
    Shifted = [0 allData(ii).Audio.MASK(1:end-1)];
    PulseBoutStart = find(allData(ii).Audio.MASK==1 & Shifted==0);
    TEMP = [diff(PulseBoutStart) 0];
    Window = (param.IgnoreSeconds + param.TransientSeconds+ param.SteadyStateSeconds)*param.nSamplesPerSeconds;
    PulseBoutStart = PulseBoutStart(TEMP>Window);
    vSamples_TransientAfterPulse = zeros(size(allData(ii).Audio.MASK));
    vSamples_SteadyStateAfterPulse = zeros(size(allData(ii).Audio.MASK));
    nRangeTransient = param.nSamplesPerSeconds*[param.IgnoreSeconds param.IgnoreSeconds+param.TransientSeconds];
    nRange_SteadyState = param.nSamplesPerSeconds*...
        [param.IgnoreSeconds+param.TransientSeconds param.IgnoreSeconds+param.TransientSeconds+param.SteadyStateSeconds];
    for nPulseBoutStart_Sample = PulseBoutStart
        vSamples_TransientAfterPulse(nPulseBoutStart_Sample+(nRangeTransient(1):nRangeTransient(2))) = 1;
        vSamples_SteadyStateAfterPulse(nPulseBoutStart_Sample+(nRange_SteadyState(1):nRange_SteadyState(2))) = 1;
    end
end


cd(DataFolder)


%Song features and kinematics
%11 song features
MaleSong_FliesKinematics(1).SongFeatures = {'Pulse amount','Sine amount','Bout amount','Pulse onsets',...
    'Sine onsets','Bout onsets','IPI','Pulse duration','Sine duration','Bout duration','Pause'};

%9 kinetic variables (4 female speed, 4 male speed, 1 distance)
MaleSong_FliesKinematics(1).FliesKinematics = {'Female forward velocity','Female lateral velocity','Female rotational velocity',...
    'Female absolute velocity','Male forward velocity','Male lateral velocity','Male rotational velocity',...
    'Male absolute velocity','Male-Female distance'};

%Saving parameters (including param.CurrTime, the time the function started
%to run) and output parameters
% disp(['Saving in ',pwd,filesep(),SaveTo])
save('/home/kyle/Desktop/output.mat','param','Folders','MaleSong_FliesKinematics')
cd(DataFolder)

disp(['Done at ',datestr(datenum(clock))])




