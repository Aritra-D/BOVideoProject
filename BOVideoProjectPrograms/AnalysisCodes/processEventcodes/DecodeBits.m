function [allBits,laststrobeval] = DecodeBits(board_dig_in_data,eventcodechs,strobech,rising)

pulseTimes = find(board_dig_in_data(strobech,:) == 1);
codeData = [];
if isempty(pulseTimes)
    trialData.codes = {[]};
    trialData.times = {[]};
    allBits.bitStr = trialData;
    return
end
count = 1;
for i = 1:length(pulseTimes)

    if ~rising %registers at last strobe value (strobe on falling edge)
        if isequal(i,length(pulseTimes))
            newPulse(count) = pulseTimes(i); %! this can result in double codes if a pulsetime is split across two rhd files
        else
            if isequal(pulseTimes(i)+1,pulseTimes(i+1)) 
            else
                newPulse(count) = pulseTimes(i);
                count = count+1;
            end
        end

    else
        if i>1
            if isequal(pulseTimes(i)-1,pulseTimes(i-1))
            else
                newPulse(count) = pulseTimes(i);
                count = count+1;
            end
        else
            newPulse(count) = pulseTimes(i);
            count = count+1; %! this can result in double codes if a pulsetime is split across two rhd files
        end
    end
end
pulseTimes = newPulse;
for i = 1:length(pulseTimes)
    for j=1:numel(eventcodechs)
    bit(i,j) = num2str(board_dig_in_data(eventcodechs(numel(eventcodechs)-j+1),pulseTimes(i)));
    end
end
for i = 1:length(pulseTimes)
    bitStr(i) = bin2dec(bit(i,:));
end

% count = 1;
% for i = 1:length(pulseTimes)
%     if~isequal(i,length(pulseTimes))
%         if ~isequal(bitStr(i),bitStr(i+1))
%             codeData(count,1) = bitStr(i);
%             codeData(count,2) = pulseTimes(i)+myData.length;
%             count = count+1;
%         end
%     end
% end
codeData(:,1) = bitStr;
codeData(:,2) = pulseTimes;
allBits.bitStr = codeData;

laststrobeval=board_dig_in_data(strobech,end);

end
