
%linearize luminance using measurements with spectroradiometer
%
%inputs:
%-lumvalmatrix: matrix of desired absolute luminance values (in cd/m2)
%-colormatrix: matrix of same size as lumval but in 3D, where the three
%layers resp. indicate RGB values. At least R or G or B channel should have
%a value of 1
%-rigid: string that specifies which rig you are working with.
%
%output:
%-relintmatrix: relative intensity matrix, 3D matrix of similar size as
%datacolor. Each value represents the relative intensity of the R,G,B
%channel

function [relintmatrix,intensities,Iout]=lin_lum(lumvalmatrix,colormatrix,rigid)

if(numel(size(lumvalmatrix))>2)
    error('N dimensions should be 2 or smaller');
elseif size(colormatrix)~=3
    error('colormatrix should be 3D');
elseif size(colormatrix,3)~=3
    error('3rd dimension of colormatrix should have a size of 3');
elseif unique(max(colormatrix,[],3))~=1
    error('max of color channels should be 1 at every point');
end

if(strcmpi(rigid,'propixxRB3D_Salk'))
%     l=load('Z:\Tom\scripts\params_propixxRB3D');
% x=l.fitmatrix{1,2};

l=load(['params_propixxRB3D_Salk']);
    xR=l.fitmatrix{1,2};
    l=load(['params_propixxRB3D_Salk']);
    xG=l.fitmatrix{1,2};
    l=load(['params_propixxRB3D_Salk']);
    xB=l.fitmatrix{1,2};
    
    warning('ignoring colors');
    
    %get calibration data
intensityprecision=0.001;
intensities=0:intensityprecision:1;

% commented this out on 4/17/2020. Luminance is actually only one eye's
% value, because each eye only receives one channel. No factor different
% than 1 needed to apply to the measured data.
% factor=3/2;
% Ireal_R=([xR(1)*(intensities.^xR(2))]+xR(3)).*factor; 
% %factor 3/2: 3 because for each channel R/G/B actually the data is used for all channels. And /2 because
% Ireal_G=([xG(1)*(intensities.^xG(2))]+xG(3)).*factor;
% Ireal_B=([xB(1)*(intensities.^xB(2))]+xB(3)).*factor;


factor=1;
Ireal_R=([xR(1)*(intensities.^xR(2))]+xR(3)).*factor; 
Ireal_G=Ireal_R.*0;
Ireal_B=Ireal_R.*0;

if numel(unique(colormatrix))~=1 || colormatrix(1)~=1
   error('color need to be grey scale for RB3D') 
end

elseif strcmpi(rigid,'DellLCD')


     %get calibration data
    l=load(['params_DellLCD_100.mat']);
    xR=l.fitmatrix{1,2};
    l=load(['params_DellLCD_010.mat']);
    xG=l.fitmatrix{1,2};
    l=load(['params_DellLCD_001.mat']);
    xB=l.fitmatrix{1,2};



intensityprecision=0.001;
intensities=0:intensityprecision:1;
Ireal_R=[xR(1)*(intensities.^xR(2))]+xR(3);
Ireal_G=[xG(1)*(intensities.^xG(2))]+xG(3);
Ireal_B=[xB(1)*(intensities.^xB(2))]+xB(3);

% %no gamma correction
% Ireal_R=[xR(1)*(intensities.^1)]+xR(3);
% Ireal_G=[xG(1)*(intensities.^1)]+xG(3);
% Ireal_B=[xB(1)*(intensities.^1)]+xB(3);

% figure;
% plot(intensities,Ireal_R+Ireal_G+Ireal_B,'k-')
% hold on;
% plot(intensities,Ireal_R,'r-')
% plot(intensities,Ireal_G,'g-')
% plot(intensities,Ireal_B,'b-')
% p=plot(intensities,Ireal_B+(Ireal_G.*0.7),'-');
% set(p,'Color',[0 0.7 1])
% p=plot(intensities,Ireal_B+Ireal_G,'-');
% set(p,'Color',[0 1 1])


elseif(strcmpi(rigid,'propixxRGB120'))
    l=load(['params_propixxRGB120_100']);
    xR=l.fitmatrix{1,2};
    l=load(['params_propixxRGB120_010']);
    xG=l.fitmatrix{1,2};
    l=load(['params_propixxRGB120_001']);
    xB=l.fitmatrix{1,2};
    
    %get calibration data
intensityprecision=0.001;
intensities=0:intensityprecision:1;
Ireal_R=[xR(1)*(intensities.^xR(2))]+xR(3);
Ireal_G=[xG(1)*(intensities.^xG(2))]+xG(3);
Ireal_B=[xB(1)*(intensities.^xB(2))]+xB(3);
elseif(strcmpi(rigid,'Asus_2P')) 
    temp=load(['\\nadata.snl.salk.edu\snlrdata\Tom\scripts\dimensions_Asus_2P']);
    gamma=temp.gamma;
    intensityprecision=0.001;
    intensities=0:intensityprecision:1;
    Ireal_R=(intensities.^gamma).*(temp.maxlum/3); %do not have individual luminance for separate RGB channels, thus divide total max luminance by 3
    Ireal_G=(intensities.^gamma).*(temp.maxlum/3);
    Ireal_B=(intensities.^gamma).*(temp.maxlum/3);
else
    error('Unknown rig');
end


% T=repmat(reshape(Ireal_R,1,1,numel(Ireal_R)),size(colormatrix,1),size(colormatrix,2)).*colormatrix(:,:,1);
% T=T+repmat(reshape(Ireal_G,1,1,numel(Ireal_G)),size(colormatrix,1),size(colormatrix,2)).*colormatrix(:,:,2);
% T=T+repmat(reshape(Ireal_B,1,1,numel(Ireal_B)),size(colormatrix,1),size(colormatrix,2)).*colormatrix(:,:,3);
% T=T-lumvalmatrix;

T=repmat(reshape(Ireal_R,1,1,numel(Ireal_R)),size(colormatrix,1),size(colormatrix,2)).*repmat(colormatrix(:,:,1),1,1,numel(Ireal_R));
T=T+repmat(reshape(Ireal_G,1,1,numel(Ireal_G)),size(colormatrix,1),size(colormatrix,2)).*repmat(colormatrix(:,:,2),1,1,numel(Ireal_G));
T=T+repmat(reshape(Ireal_B,1,1,numel(Ireal_B)),size(colormatrix,1),size(colormatrix,2)).*repmat(colormatrix(:,:,3),1,1,numel(Ireal_B));
T=T-repmat(lumvalmatrix,1,1,size(T,3));

T=abs(T);
[minvals,mininds]=min(T,[],3);
if max(max(minvals))>1
    error('error of more than 1 cd/m2');
end
relintmatrix=intensities(mininds);

% %generate a xx4 matrix from lumvals and colorvals to check for unique
% %elements
% totmatrix(:,:,1)=colormatrix(:,:,1);
% totmatrix(:,:,2)=colormatrix(:,:,2);
% totmatrix(:,:,3)=colormatrix(:,:,3);
% totmatrix(:,:,4)=lumvalmatrix;
% a1=totmatrix(:,:,1);
% a2=totmatrix(:,:,2);
% a3=totmatrix(:,:,3);
% a4=totmatrix(:,:,4);
% v=cell2mat(arrayfun(@(x1,x2,x3,x4) [x1 x2 x3 x4],a1(:),a2(:),a3(:),a4(:),'un',0));
% [unrows,~,~]=unique(v,'rows','stable');
% %for those unique elements, determine the relative intensity required to
% %get the desired luminance. Then fill out the output matrix
% relintmatrix=ones(size(lumvalmatrix)).*NaN;
% for i=1:size(unrows,1)
%     currrow=unrows(i,:);
%     currcalibfunction=(Ireal_R.*currrow(1))+(Ireal_G.*currrow(2))+(Ireal_B.*currrow(3));
%     [minval,minind]=min(abs(currcalibfunction-currrow(4)));
%     if minval>1
%        error('error of more than 1 cd/m2') 
%     end
%     
%     relintmatrix(totmatrix(:,:,1)==currrow(1) & totmatrix(:,:,2)==currrow(2) & totmatrix(:,:,3)==currrow(3)  & totmatrix(:,:,4)==currrow(4))=intensities(minind);
% end

Iout = Ireal_R+Ireal_G+Ireal_B;

end