% make 3D RGB image matrix from 2D template with luminance values
function imfin=makeimfin(im,currcolor,rigid,lums,colorlum)

maxlum=max(lums);
minlum=min(lums);

%make colormatrix
colormatrix=ones(size(im,1),size(im,2),3).*0;
if colorlum==maxlum
    colormatrix(:,:,1)=colormatrix(:,:,1)+(im==max(max(im))).*currcolor(1);
    colormatrix(:,:,2)=colormatrix(:,:,2)+(im==max(max(im))).*currcolor(2);
    colormatrix(:,:,3)=colormatrix(:,:,3)+(im==max(max(im))).*currcolor(3);
    colormatrix(:,:,1)=colormatrix(:,:,1)+(im==min(min(im))).*1;
    colormatrix(:,:,2)=colormatrix(:,:,2)+(im==min(min(im))).*1;
    colormatrix(:,:,3)=colormatrix(:,:,3)+(im==min(min(im))).*1;
elseif colorlum==minlum
    colormatrix(:,:,1)=colormatrix(:,:,1)+(im==max(max(im))).*1;
    colormatrix(:,:,2)=colormatrix(:,:,2)+(im==max(max(im))).*1;
    colormatrix(:,:,3)=colormatrix(:,:,3)+(im==max(max(im))).*1;
    colormatrix(:,:,1)=colormatrix(:,:,1)+(im==min(min(im))).*currcolor(1);
    colormatrix(:,:,2)=colormatrix(:,:,2)+(im==min(min(im))).*currcolor(2);
    colormatrix(:,:,3)=colormatrix(:,:,3)+(im==min(min(im))).*currcolor(3);
else
    error('color lum should be either min or max lum');
end
clear im2 imB1 imB2 goodrows x y;
%get all unique color+luminance combinations and do
%lin_lum for each one separate (speeds up lin_lum a
%lot for large arrays)
clear temp;
temp(:,:,1)=colormatrix(:,:,1);
temp(:,:,2)=colormatrix(:,:,2);
temp(:,:,3)=colormatrix(:,:,3);
temp(:,:,4)=im;
unlumcolors=unique(reshape(temp,[size(temp,1)*size(temp,2),size(temp,3),1]),'rows');
clear relint;
relintmatrix=ones(size(im)).*NaN;
for uli=1:size(unlumcolors,1)
    tempcolor=reshape(unlumcolors(uli,1:3),1,1,3);
    templum=unlumcolors(uli,4);
    temprelint=lin_lum(templum,tempcolor,rigid);
    relintmatrix(im==templum & colormatrix(:,:,1)==tempcolor(:,:,1) & colormatrix(:,:,2)==tempcolor(:,:,2) & colormatrix(:,:,3)==tempcolor(:,:,3))=temprelint;
end
imfin(:,:,1)=colormatrix(:,:,1).*relintmatrix;
imfin(:,:,2)=colormatrix(:,:,2).*relintmatrix;
imfin(:,:,3)=colormatrix(:,:,3).*relintmatrix;

end