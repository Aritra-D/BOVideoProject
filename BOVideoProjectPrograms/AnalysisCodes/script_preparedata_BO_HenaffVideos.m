% first run script_processalbarcodes_singledataset to get timing
%information.
clear all;

addpath(genpath('E:\Spikesorting\bombcell'))

clear D;

mergeddatasets_order={'VideoAndBo'}; %order of files as in merge
currdate='20240523';
datapath_ksdata = ['E:\Data_local\' currdate '\VideoAndBo\kilosort4']; %folder of kilosort output
data_local = ['E:\Data_local\' currdate '\VideoAndBo']; %folder with data. Should contain bhv2 file and ap.bin file and ap.meta file, and *eventcodes*.mat file with the processed eventcodes (from script_processallbarcodes_singledataset.m)

n=0;

try
    load([data_local '\alleventcodes_np.mat'],'alleventcodes_np','totalNsamples','currdate');

catch ME
    load([data_local '\eventcodes_' mergeddatasets_order{1} '.mat'],'eventcodes_np','nFileSamp'); %this line is suitable for data sets with only a single task (i.e. no merge done)
    alleventcodes_np=eventcodes_np;
    totalNsamples=nFileSamp;
end

MLfilename=dir([data_local '\*bhv2']);
if numel(MLfilename)~=1, error('no unique bhv file'); end
MLfilename=MLfilename.name;
[b,a]=mlread([data_local '\' MLfilename]);
ppd=a.PixelsPerDegree(1);
allcodes=[];
for i=1:numel(b)
    allcodes=[allcodes;b(i).BehavioralCodes.CodeNumbers];
end
startpos=strfind(alleventcodes_np(:,1)',allcodes');
if numel(startpos)~=1, error('codes not found'); end
alleventcodes_np=alleventcodes_np(startpos:startpos+numel(allcodes)-1,:); %only retain codes of the task of interest
allcodes_np_pertrial=[];
starttrial_pos=find(alleventcodes_np(:,1)==9);
if numel(starttrial_pos)~=numel(b), error('number of trials do not seem to be the same'); end
for i=1:numel(starttrial_pos)

     if i==1
         disp('Reproducing stimuli...')
        m=mlexportstim([data_local],[data_local '\' MLfilename]);
     end

    if b(i).TrialError==0


    if i==numel(starttrial_pos)
        endpos=size(alleventcodes_np,1);
    else
        endpos=starttrial_pos(i+1)-1;
    end
    tempdata=[alleventcodes_np(starttrial_pos(i):endpos,:)];

    %video trial or BO trial
    clear trialtype;
    if numel(strfind(b(i).TaskObject.Attribute{3}{2},'squareR'))>0
        trialtype='BO';
    elseif numel(strfind(b(i).TaskObject.Attribute{4}{2},'video'))>0
        trialtype='video';
    end

    conditions=[];
    if strcmpi(trialtype,'BO')
        %get conditions of the sequential boards
        [rowsi,colsi]=find(b(i).ObjectStatusRecord.Status(:,2:end-2)==1);
        [~,sorti]=sort(rowsi);
        objectindices=colsi(sorti)+1; %order in which the different scenes appear
        tempi=find(tempdata(:,1)==101);
        tempi_ends=find(tempdata(:,1)==26);
        for nsi=1:numel(tempi)
            currtime_start=tempdata(tempi(nsi),2);
            currtime_end=tempdata(tempi_ends(nsi),2);
            currobject=b(i).TaskObject.Attribute{objectindices(nsi)}{2};
            currobject=currobject(max(strfind(currobject,'\'))+1:strfind(currobject,'.')-1);
            im=imread([data_local '\' currobject '.bmp']);
            bglum=mode(reshape(im,[],1));
            squarelum=setdiff(unique(im),bglum);
            squareor=str2num(currobject(strfind(currobject,'squareR')+7:strfind(currobject,'S')-1));

            clear tempstruct;
            tempstruct.starttime=currtime_start;
            tempstruct.endtime=currtime_end;
            tempstruct.squarelum=squarelum;
            tempstruct.bglum=bglum;
            tempstruct.squareor=squareor;

            conditions=[conditions {tempstruct}];
        end
    elseif strcmpi(trialtype,'video')
        currtime_start=tempdata(tempdata(:,1)==23,2);
        currtime_end=tempdata(tempdata(:,1)==117,2);
        temp=b(i).TaskObject.Attribute{4}{2};
        currobject=temp(max(strfind(temp,'\'))+1:end);

        clear tempstruct;
        tempstruct.starttime=currtime_start;
        tempstruct.endtime=currtime_end;
        tempstruct.filename=currobject;

        conditions=[conditions {tempstruct}];
    end

    allcodes_np_pertrial=[allcodes_np_pertrial; {i} {tempdata} {conditions} {trialtype}];

    end
end

allcodes_np_pertrial0=allcodes_np_pertrial;

%get spike times
spike_times_phy=readNPY([datapath_ksdata '\spike_times.npy']);
spike_clusters_phy=readNPY([datapath_ksdata '\spike_clusters.npy']);
spike_templates_phy=readNPY([datapath_ksdata '\spike_templates.npy']); %this will be the same as spike_clusters_phy if not manually curated in phy
channelpos=readNPY([datapath_ksdata '\channel_positions.npy']); %coordinates of electrodes on probe

try
    % find units labeled good by bombcell
    [param, qMetric] = bc_loadSavedMetrics(datapath_ksdata);
    unitType = bc_getQualityUnitType(param, qMetric, datapath_ksdata);
    goodUnits = unitType == 1; %SU
    muaUnits = unitType == 2;
    noiseUnits = unitType == 0;
    nonSomaticUnits = unitType == 3;

    % whichUnits = goodUnits|muaUnits|nonSomaticUnits;
    whichUnits = goodUnits|nonSomaticUnits

    selectedUnits_phyinds=qMetric.phy_clusterID(whichUnits); %(phy:0-indexed)
    selectedUnits_maxchannels=qMetric.maxChannels(whichUnits);
    selectedUnits_depths=channelpos(selectedUnits_maxchannels, 2); %depth position on probe of max channel for selected units

catch ME
    warning('Bombcell output not found, using kilosort output directly');
    t = readtable([datapath_ksdata '\cluster_KSLabel.tsv'], "FileType","text",'Delimiter', '\t'); %read KS cluster labels

    spike_positions=readNPY([datapath_ksdata '\spike_positions.npy']); %this does not work with Kilosort 2.5

    selectedUnits_phyinds=t{strcmpi(t{:,2},'good'),1}; %only include units labeled as 'good' by KS
    selectedUnits_depths=nan(size(selectedUnits_phyinds));
    for tempi=1:numel(selectedUnits_phyinds)
        allspikes=(find(spike_clusters_phy==selectedUnits_phyinds(tempi)));
        currdepth=double(mean(spike_positions(allspikes,2))); %compute depth as average depth (?) across all spike occurrences. Probably these are not identical numbers because of drift correction?
        selectedUnits_depths(tempi)=currdepth;
    end

end


Fs=30000;

pre_interval=0; %ms, how long to look prior to board comes on
post_interval_board=500; %ms, how long to look after board comes on
post_interval_video=3000; %ms, how long to look after board comes on
pre_interval_s=round(pre_interval.*1e-3.*Fs);
post_interval_board_s=round(post_interval_board.*1e-3.*Fs);
post_interval_video_s=round(post_interval_video.*1e-3.*Fs);

selectedtrials=allcodes_np_pertrial;

allspiketimes_perunit=[];
for ni=1:numel(selectedUnits_phyinds)
            %get spike times from the right unit
            allspiketimes=double(spike_times_phy(spike_clusters_phy==selectedUnits_phyinds(ni)));
            allspiketimes_perunit=[allspiketimes_perunit {allspiketimes}];
end

for trialindex=1:size(selectedtrials,1)
    disp(['Trial ' num2str(trialindex) ' out of ' num2str(size(selectedtrials,1))])

    currtrial=selectedtrials(trialindex,:);
    currcodes=currtrial{2};
    currconditions=currtrial{3};
    currtrialtype=currtrial{4};

    clear currpost_interval_s;
    if strcmpi(currtrialtype,'BO')
    currpost_interval_s = post_interval_board_s;
    elseif strcmpi(currtrialtype,'video')
    currpost_interval_s = post_interval_video_s;
    end

    for ci=1:numel(currconditions)
        currcondition=currconditions{ci};

        startpos=currcondition.starttime-pre_interval_s;
        endpos=currcondition.starttime+currpost_interval_s;

        binnedspikes=zeros(numel(selectedUnits_phyinds),round(((currpost_interval_s-pre_interval_s)./Fs).*1e3)+1);
        for ni=1:numel(selectedUnits_phyinds)
            %get spike times from the right unit
            allspiketimes=allspiketimes_perunit{ni};
            currdepth=selectedUnits_depths(ni);
            temp=allspiketimes(allspiketimes>=startpos & allspiketimes<endpos);
            temp_t=round(((temp-startpos)./Fs).*1e3); %round spike to nearest ms relative to startpos
            binnedspikes(ni,(temp_t+1))=1; %bin 0 ms is the first bin, thus do +1 for index
        end

        n=n+1;
        D(n).data=binnedspikes;

        currfieldnames=fieldnames(currcondition);

        for cfi=1:numel(currfieldnames)
            if ~strcmpi(currfieldnames{cfi},'starttime') && ~strcmpi(currfieldnames{cfi},'endtime')
            eval(['D(n).' currfieldnames{cfi} '=currcondition.'  currfieldnames{cfi} ';']);
            end
        end
    end
end

save([data_local '\' currdate '_processed.mat'],'D','selectedUnits_depths','selectedUnits_phyinds');
% save([data_local '\' currdate '_processed.mat'],'D','selectedUnits_depths','selectedUnits_phyinds','-v7.3');