function timestamps_new=alignTimestampsFromBarcodes(timestamps_orig,barcodes_new,barcodes_orig)

%limit barcodes to common codes
commonbarcodes=intersect(barcodes_new(:,2),barcodes_orig(:,2)); %second column are barcodes
barcodes_new=barcodes_new(ismember(barcodes_new(:,2),commonbarcodes),:);
barcodes_orig=barcodes_orig(ismember(barcodes_orig(:,2),commonbarcodes),:);
%Determine slope (m) between main/secondary timestamps
m=(barcodes_new(end,1)-barcodes_new(1,1))/(barcodes_orig(end,1)-barcodes_orig(1,1)); %first column in barcodes are time stamps
%Determine offset (b) between main and secondary barcode timestamps
b = barcodes_new(1,1) - barcodes_orig(1,1) * m;
%use slope and offset to align timestamps
timestamps_new=round((timestamps_orig.*m)+b);