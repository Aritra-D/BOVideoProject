% 20240216: corrections in bar code reading

function allbarcodes = getBarCodes(barcodedata,Fs)

global_tolerance=0.2; %the fraction (percentage) of tolerance allowed for duration measurements.
ind_wrap_duration = 20; %In milliseconds - duration of OFF-ON-OFF wrapper at start of each barcode (i.e. 3x ind_wrap_duration)
ind_bar_duration = 29; %ms - duration of each bar
notherbits=24; %not including ON wrapper
mindur_betweencodes = 10000; %ms
barcodethr=mean([min(barcodedata) max(barcodedata)]);
barcodedata_digit=barcodedata>barcodethr;
f=find(barcodedata_digit==1);
f2=find(diff(f)>1);
bitchanges=sort([f([1 f2]) f(f2+1)]); %all position where barcode channel changes from 0 to 1 or from 1 to 0

min_wrap_duration = ind_wrap_duration - ind_wrap_duration * global_tolerance;
max_wrap_duration = ind_wrap_duration + ind_wrap_duration * global_tolerance;
min_wrap_duration_s=min_wrap_duration.*1e-3.*Fs;
max_wrap_duration_s=max_wrap_duration.*1e-3.*Fs;
wrap_duration_s = round(ind_wrap_duration.*1e-3.*Fs);
min_bar_duration = ind_bar_duration - ind_bar_duration  * global_tolerance;
max_bar_duration = ind_bar_duration  + ind_bar_duration  * global_tolerance;
min_bar_duration_s = round(min_bar_duration .*1e-3.*Fs);
max_bar_duration_s =round(max_bar_duration.*1e-3.*Fs);
bar_duration_s = round(ind_bar_duration.*1e-3.*Fs);

mindur_betweencodes_s =mindur_betweencodes.*1e-3.*Fs;

t1=diff(bitchanges)>min_wrap_duration_s & diff(bitchanges)<max_wrap_duration_s;
t2=diff(bitchanges)>mindur_betweencodes_s;
wrapONstarts=round(bitchanges([1 find(t1 & [0 t2(1:end-1)])]));

t1=diff(bitchanges)>min_bar_duration_s & diff(bitchanges)<max_bar_duration_s;
t2=diff(bitchanges)>mindur_betweencodes_s;
wrapOFFstarts=round(bitchanges([find(t1 & [t2(2:end) 0]) numel(bitchanges)]));

if numel(wrapOFFstarts)>numel(wrapONstarts)
wrapOFFstarts=wrapOFFstarts(1:numel(wrapONstarts));
end

% figure;plot(barcodedata_digit);hold on;plot(wrapONstarts,1,'go');hold on;plot(wrapOFFstarts,1,'rx');

wrap_duration_s=bitchanges(find(bitchanges==wrapONstarts(1))+1)-wrapONstarts(1);

%read barcodes
allbarcodes=[];
for j=1:numel(wrapONstarts)
    currbarcode=[0 1 0]; %barcode wrapper starts with 010
    code_end=wrapOFFstarts(j);

    temp=bitchanges(bitchanges>wrapONstarts(j) & bitchanges<wrapOFFstarts(j));
    for ni=2:numel(temp) %start at 2 because the first crossing will be change from 1 to 0 of wrapON
        if ni==2
            tempcount=round((temp(ni)-temp(ni-1)-wrap_duration_s)/bar_duration_s); %subtract 1x wrap duration because the last '0' phase of the 010 wrapON needs to be removed
        else
            tempcount=round((temp(ni)-temp(ni-1))/bar_duration_s);
        end
        if tempcount>0
            currval=barcodedata_digit(temp(ni)-1);
            tempcodes=ones(1,tempcount).*currval;
            currbarcode=[currbarcode tempcodes];
        end
    end
    tempcount=round((wrapOFFstarts(j)-temp(ni))/bar_duration_s); %code duration at OFF is same as in the middle of the barcode -> wrap duration subtraction not required
    currval=barcodedata_digit(wrapOFFstarts(j)-1);
    tempcodes=ones(1,tempcount).*currval;
    currbarcode=[currbarcode tempcodes barcodedata_digit(wrapOFFstarts(j)+1)];
%     if numel(currbarcode)~=(notherbits+3)
%         error('wrong number of digits in barcode');
%     end
    % figure;plot(barcodedata_digit);hold on;plot(wrapONstarts,1,'go');hold on;plot(wrapOFFstarts,1,'rx'); hold on;plot(temp,0.5,'mx');hold on;plot([temp(1)+wrap_duration_s:bar_duration_s:wrapOFFstarts(j)],0.5,'ko'); xlim([wrapONstarts(j) wrapOFFstarts(j)])
    allbarcodes=[allbarcodes;wrapONstarts(j) bin2dec(strrep(num2str([currbarcode]),' ',''))]; %first column is time start in samples of wrap ON sequence
end

end