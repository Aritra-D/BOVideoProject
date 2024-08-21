function [barcodes, eventcodes] = readIntanCodes(datapath,basefilename,rising)

photodiodech=1;
barcodech=2;
strobech=8;
eventcodechs=1:7;
fileInfo = dir(fullfile(datapath,[basefilename '*.rhd']));
%sort filenames by date and time
allnumbers=[];
for j=1:numel(fileInfo)
    currname=fileInfo(j).name;
    underscores=strfind(currname,'_');
    period=strfind(currname,'.');
    currnumber=str2num([currname(underscores(1)+1:underscores(2)-1) currname(underscores(2)+1:period-1)]);
    allnumbers=[allnumbers; currnumber];
end
[~,sorti]=sort(allnumbers,'ascend');
fileInfo=fileInfo(sorti);
%read barcodes, photodiode data and event codes
photodiodedata=[];barcodedata=[];eventcodes=[];
clear laststrobeval_prev;
currtimeoffset=0;
for j=1:numel(fileInfo)
    [board_dig_in_data,board_adc_data,Fs]=read_Intan_RHD2000_file([fileInfo(j).folder '\'],fileInfo(j).name);
    [allBits,laststrobeval] = DecodeBits(board_dig_in_data,eventcodechs,strobech,rising);
    eventtemp=allBits.bitStr;
    if rising %make sure detected strobes are actually rising
        if j>1
            if laststrobeval_prev==1 && eventtemp(1,2)==1 %then not rising
                eventtemp=eventtemp(2:end,:);
            end
        else
            if eventtemp(1,2)==1 %then not rising
                eventtemp=eventtemp(2:end,:);
            end
        end
    else
        error('not implemented yet for falling')
    end
    photodiodedata=[photodiodedata board_adc_data(photodiodech,:)];
    barcodedata=[barcodedata board_adc_data(barcodech,:)];
    eventcodes=[eventcodes;eventtemp(:,1) [eventtemp(:,2)+currtimeoffset]];
    laststrobeval_prev=laststrobeval;
    currtimeoffset=currtimeoffset+size(board_dig_in_data,2);
    clear board_adc_data board_dig_in_data;
end

% read barcodes
barcodes = getBarCodes(barcodedata,Fs);

end