% script will load barcodes for all data sets from Intan, and from ap.bin
% files, compare them, and uses them to generate time stamps in Neuropixel
% time for all MonkeyLogic eventcodes. Result is saved in the matrix
% eventcodes_np.
% This version intended for sort of a single data set (not merged).
% See E:\NeuropixelDataprocessing\ for version for merged datasets

% currdataset='PassiveReg'; %name that will be given to output file
% currdate = '20240329'; %recording date
% data_path = ['E:\Data_local\' currdate]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file
% 
% currdataset='PassiveArgyle'; %name that will be given to output file
% currdate = '20240417'; %recording date
% data_path = ['E:\Data_local\' currdate]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file
% 
% currdataset='PassiveArgyle'; %name that will be given to output file
% currdate = '20240423'; %recording date
% data_path = ['E:\Data_local\' currdate '\' currdataset]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file
% 
% currdataset='PassiveReg'; %name that will be given to output file
% currdate = '20240223'; %recording date
% data_path = ['E:\Data_local\' currdate '\' currdataset]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file
% 
% currdataset='VideoAndBo2'; %name that will be given to output file
% currdate = '20240516'; %recording date
% data_path = ['E:\Data_local\' currdate '\' currdataset]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file

currdataset='VideoAndBo'; %name that will be given to output file
currdate = '20240523'; %recording date
folderSourceString = 'C:\Users\Aritra\Downloads\Lab Workbench\Projects\NaturalVideo-BO\data\rawData\';
data_path = [folderSourceString currdate '\' currdataset]; %folder with data. Should contain intan files, bhv2 file and ap.bin file and ap.meta file

%% process barcodes and extract event codes in neuropixels time for each file
    try
        %check if the process was done already
        load([data_path '\eventcodes_' currdataset '.mat'],'eventcodes_intan');
    catch ME
        if strcmpi(ME.identifier,'MATLAB:load:couldNotReadFile')
            % load MonkeyLogic data and extract codes
            filenameML=dir([data_path '\*.bhv2']);
            if numel(filenameML)~=1, error('no or more than 1 file found'); end
            [allcodes_ML,rising,code_start,code_end] = getMonkeyLogicCodes(filenameML.folder,filenameML.name);
            % load Intan data, extract barcodes and eventcodes
            [barcodes_intan, eventcodes_intan] = readIntanCodes([data_path '\'],currdate,rising);
            eventcodes_intan=eventcodes_intan(1:size(allcodes_ML,1),:);
            if sum(mod(eventcodes_intan(:,1),128)==mod(allcodes_ML(:,1),128)) ~= size(allcodes_ML,1), error('event codes do not match'); end
            save([data_path '\eventcodes_' currdataset '.mat'],'eventcodes_intan','barcodes_intan');
        else
            error('unknown error')
        end
        %get barcodes from NPX file by reading in chunks
        file_np=dir([data_path '\*ap.bin']);
        if numel(file_np)~=1, error('no or more than 1 file found'); end
        metadata=SGLX_readMeta.ReadMeta(file_np.name,file_np.folder);
        Fs = str2num(metadata.imSampRate);
        nChan = str2num(metadata.nSavedChans);
        nFileSamp = str2double(metadata.fileSizeBytes) / (2 * nChan); %formula from from SGLX_readMeta.m
        chunksize = 0.2; %in sec
        nSamp0=round(chunksize*Fs);
        nchunks = ceil(nFileSamp/nSamp0);
        barcodechannel=nan(nFileSamp,1);
        fid = fopen([file_np.folder '\' file_np.name], 'r');
        samp0=0;
        for j=1:nchunks
            if mod(j,100)==0, disp([num2str(j) ' out of ' num2str(nchunks)]); end
            samp0 = max(samp0, 0);
            nSamp = min(nSamp0, nFileSamp - samp0);
            fseek(fid, samp0 * 2 * nChan, 'bof');
            dataArray = fread(fid, [nChan, nSamp], 'int16=>double');
            barcodechannel(samp0+1:samp0+1+size(dataArray,2)-1)=dataArray(end,:);
            samp0=samp0+size(dataArray,2);
        end
        fclose(fid);
        
        %correction necessary for data 20240419: two samples with high values in barcodechannel that mess up barcode detection -> remove those 
        f=find(barcodechannel>99);
        for fi=1:numel(f)
            barcodechannel(f(fi))=barcodechannel(f(fi)-1);
        end
        
        barcodes_np = getBarCodes(barcodechannel',Fs);
        timestamps_np=alignTimestampsFromBarcodes(eventcodes_intan(:,2),barcodes_np,barcodes_intan);
        eventcodes_np=[eventcodes_intan(:,1) timestamps_np]; %first column is event codes, second column is time stamps in neuropixels time values (samples)
        save([data_path '\eventcodes_' currdataset '.mat'],'barcodes_np','eventcodes_np','nFileSamp','-append');
    end