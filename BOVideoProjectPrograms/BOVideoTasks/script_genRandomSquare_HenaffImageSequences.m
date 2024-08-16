% This program creates BO stimuli and Henaff Image sequences
% Videostims_Henaffetal2021.mat (~866MB) needs to be added to the 
% ..path../BOVideoTasks/ to run script_genRandomSquare_HenaffImageSequences. 

clear;

fileNameWithPath = mfilename('fullpath');
[filePath,name,ext] = fileparts(fileNameWithPath);
cd(filePath);
addpath(genpath(pwd));

displayname='propixxRGB120';
display=load(['dimensions_' displayname '.mat']);

% increase pixels while generating stimuli
display_highres=display;
% display_highres.res=[3840 2160]; %upsample images during design
magfactor=unique(display_highres.res./display.res);
[~,ppd_highres]=angle2pix(display_highres,1);
% ppd = ppd_highres;
basedatapath = filePath;

t=tic;

%% delete previous stimuli
delete(fullfile(basedatapath,'*Henaff*.bmp'));
delete(fullfile(basedatapath,'*square*.bmp'));
% delete(fullfile(basedatapath,'**.bmp'));
% delete(fullfile(basedatapath,'*video*.mp4'));
% delete(fullfile(basedatapath,'*video*.avi'));

%% USER ADJUSTABLE PARAMETERS
RFpos_real = [-1 -0.5]; %x and y ccoordinateds of RF center in ddegrees

repetitions = 50; %how often a single scene is repeated

%video stimuli parameters
timingfilename_video='fixationtask_video_sceneframework';
sizerange_dva = [2 10]; %range of possible widths/heights of a rectangle
maxspeed = 0.3; %max speed (dva per frame)
maxaccel = 0; %max acceleration
scenesize_dva = [20 20]; %dva
bgval=round(mean([0 255]));
mincontrastdifferencewithbg = 0.1;% Michelson
nsurfaces = randi([10,15],1); %number of rectangles
videodur = 3000; %duration of video (ms)
framerate = 30; %frames per second. 30 f/s is Allen video framerate.
nsurfacerange=[5,10]; %range of potential number of surfaces
nvideos = 5; %number of different videos
createnewvideos=1;
n_allen_videos = 4;

%BO parameters
colors_BO=[1 1 1];
rotangles=linspace(0,180,5); %border orientations
rotangles=rotangles(1:end-1);
% rotangles=0;
sidelengths=[... % border length
    17 17;...
    ];
lumcontrast = 0.537; %luminance contrast accross border 0.537 in O Herron and von der Heydt
scenespertrial = 6; %number of square scenes shown in each trial (must be 1 if npos>1)
stimdur = 200;
blankdur = 100; % duration of blank interval in between stims
positions = [... % x and y coordinates in dva, relative to center of RF
    0 0;...
    ].*3;
% positions = [... % x and y coordinates in dva, relative to center of RF
%     0 0;...
%     1 0;...
%     -1 0;...
%     0 1;...
%     0 -1;...
%     ].*3;
timingfilename='BOSsquaretest_basic';


%% RETRIEVE LUMINANCE INFORMATION FROM SCREEN
[~,ppd]=angle2pix(display,1);
maxdims=display.res;

l = load(['params_' displayname '_100.mat']);
xR = l.fitmatrix{1,2};
l = load(['params_' displayname '_010.mat']);
xG = l.fitmatrix{1,2};
l = load(['params_' displayname '_001.mat']);
xB = l.fitmatrix{1,2};
maxlum = xR(1)+xR(3)+xG(1)+xG(3)+xB(1)+xB(3);
minlum = xR(3)+xG(3)+xB(3);
lumsrange = [minlum maxlum];

%% generate video stimuli with randomly moving surfaces

RF_relvideo = [0 0]; %RF relative to center of video frame

forbiddencolors = round(((1-mincontrastdifferencewithbg)/(1+mincontrastdifferencewithbg))*bgval):round(((1+mincontrastdifferencewithbg)/(1-mincontrastdifferencewithbg))*bgval);

scenesize_dva0 = scenesize_dva;
nsurfacerange0 = nsurfacerange;
scenesize0 = round(scenesize_dva.*ppd);

% if createnewvideos
%     % increase size before generating stimuli, so that there is more uniform
%     % distribution of elements moving over the scene
%     scenesize_dva=scenesize_dva.*3;
%     nsurfacerange=nsurfacerange.*9;
%
%     nframes = round(videodur*framerate/1e3);
%     scenesize=round(scenesize_dva.*ppd);
%     sizerange = round(sizerange_dva.*ppd); %pixels
%     for nvi=1:nvideos
%         disp(nvi)
%         squarewidths = []; squareheights=[];
%         xs=[];ys=[];dxs=[];dys=[];colors=[];ddxs=[];ddys=[];
%         nsurfaces = randi(nsurfacerange,1); %number of rectangles
%         for i = 1:nsurfaces
%             squarewidths = cat(2,squarewidths,randi(sizerange,1));
%             squareheights = cat(2,squareheights,randi(sizerange,1));
%             colors = cat(2,colors,randsample(setdiff(0:1:255,forbiddencolors),1)); % random color, but cannot be too close to bgcolor
%             xs = cat(2,xs,randi([1 scenesize(2)],1)); % initial square position (left top corner column)
%             ys = cat(2,ys,randi([1 scenesize(1)],1)); % initial square position (left top corner row)
%             dxs = cat(2,dxs,randi(round([-1 1].*maxspeed.*ppd),1)); % how much pixels the x coordinate of a square changes in each frame
%             dys = cat(2,dys,randi(round([-1 1].*maxspeed.*ppd),1)); % how much pixels the y coordinate of a square changes in each frame
%             ddxs = cat(2,ddxs,randi(round([-1 1].*maxaccel.*ppd),1)); % how much pixels the x coordinate of a square changes in each frame
%             ddys = cat(2,ddys,randi(round([-1 1].*maxaccel.*ppd),1)); % how much pixels the y coordinate of a square changes in each frame
%         end
%         scene0 = ones(scenesize).*bgval;
%         allframes=[];
%         for nfi=1:nframes
%             disp(nfi)
%             scene=scene0; %start with an empty frame before the first square is added
%             for i=1:nsurfaces
%                 if i>1 && mod(i,3)==0
%                     dxs(i) = dxs(i) + ddxs(i); %#ok<*SAGROW>
%                     if abs(dxs(i))>=max(round(abs(maxspeed.*ppd)))
%                         ddxs(i)=ddxs(i).*-1;
%                     end
%                     dys(i) = dys(i)+ddys(i);
%                     if abs(dys(i))>=max(round(abs(maxspeed.*ppd)))
%                         ddys(i)=ddys(i).*-1;
%                     end
%                 end
%                 scene=addSquare(scene,squarewidths(i),squareheights(i),xs(i)+(dxs(i)*(nfi-1)),ys(i)+(dys(i)*(nfi-1)),colors(i));
%             end
%             scene_trimmed=scene(scenesize0(1)+1:(scenesize0(1)*2),scenesize0(1)+1:(scenesize0(1)*2));
%             allframes = cat(2,allframes,{scene_trimmed});
%         end
%         v=VideoWriter(fullfile(basedatapath,['video' num2str(nvi)],'Uncompressed AVI'));
%         v.FrameRate = framerate;
%         open(v);
%         for i=1:numel(allframes)
%             writeVideo(v,allframes{i}./255);
%         end
%         close(v);
%     end
% end

%% MAKE FIX TARGET
fixluminance_center = maxlum; %cd/m2
fixluminance_edge = mean(lumsrange); %cd/m2
fixtarget_width_total = 0.5951; %dva
fixtarget_width_center = 0.255; %dva
fixfilename = 'fixtarget';
[~, total_pix] = angle2pix(display,fixtarget_width_total);
[~, center_pix] = angle2pix(display,fixtarget_width_center);
total_pix = round(total_pix);
center_pix = round(center_pix);
edge_pix = (total_pix-center_pix)/2;
imt = ones(total_pix,total_pix,1).*fixluminance_edge;
imc = ones(center_pix,center_pix,1).*fixluminance_center;
imt((edge_pix+1):(edge_pix+1)+center_pix-1,(edge_pix+1):(edge_pix+1)+center_pix-1)=imc;
imt_blackcenter = imt;
imt_blackcenter((edge_pix+1):(edge_pix+1)+center_pix-1,(edge_pix+1):(edge_pix+1)+center_pix-1)=5;
relintmatrix = lin_lum(imt,ones(size(imt,1),size(imt,2),3),displayname);
relintmatrix_blackcenter = lin_lum(imt_blackcenter,ones(size(imt_blackcenter,1),size(imt_blackcenter,2),3),displayname);
matrix(:,:,1) = relintmatrix;
matrix(:,:,2) = relintmatrix_blackcenter;
matrix(:,:,3) = relintmatrix_blackcenter;
imwrite(matrix,fullfile(basedatapath, [fixfilename '.bmp']));
clear matrix;

%% MAKE BG FRAME
bgluminance = mean(lumsrange); %cd/m2
bgfilename = 'bgframe';
imt = ones(scenesize0(1),scenesize0(2),1).*bgluminance;
relintmatrix = lin_lum(imt,ones(size(imt,1),size(imt,2),3),displayname);
matrix(:,:,1) = relintmatrix;
matrix(:,:,2) = relintmatrix;
matrix(:,:,3) = relintmatrix;
imwrite(matrix,fullfile(basedatapath, [bgfilename '.bmp']));
clear matrix;


%% DELETE PREVIOUS BO SCENES
delete([basedatapath '*square*.bmp']);

%% MAKE BO STIMULI
timingfilename_BO = 'BOSsquaretest_withinvideotask';

shapename='square';

allscenes=[];

for ci = 1:size(colors_BO,1)

    currcolor=colors_BO(ci,:);
    maxlum=xR(1).*currcolor(1)+xR(3)+xG(1).*currcolor(2)+xG(3)+xB(1).*currcolor(3)+xB(3);
    minlum=xR(3)+xG(3)+xB(3);
    lumsrange = [minlum maxlum];

    % keep bg lum fixed
    [minlum,maxlum,lumcontrast,bglum]=getLumsAndContrast(NaN,NaN,lumcontrast,mean(lumsrange));
    lums=[maxlum minlum]; %cd/m2
    colorlums=lums;

    if sum(colors_BO==[1 1 1])==3
        maxcli=1;
    else
        maxcli=2;
    end

    for cli=1:maxcli %different luminances

        colorlum=colorlums(cli);
        colorlumname=['CL' num2str(round(colorlum))]; % luminance of the 'color' side

        for si=1:size(sidelengths,1)

            clear im imfin;
            currsidelength = sidelengths(si,1);
            currheight = currsidelength;
            lengthname=['S' num2str(currsidelength)];

            if currsidelength ~= sidelengths(si,2)
                currheight = sidelengths(si,2);
                lengthname=['S' num2str(currsidelength) 'h' num2str(currheight)];
            end

            for ri = 1:numel(rotangles)

                currrotangle = rotangles(ri);

                rotname = ['R' num2str(round(currrotangle))];

                clear imfin;

                currlum = lums(1);
                bgval = lums(2);

                placeholderval = -1;

                if currlum==colorlum
                    figname ='C'; %figure is colorside
                else
                    figname ='G'; %background is colorside
                end

                [~,sidelength_pix] = angle2pix(display,currsidelength);
                sidelength_pix = round(sidelength_pix);
                [~,currheight_pix] = angle2pix(display,currheight);
                currheight_pix = round(currheight_pix);

                im = ones(currheight_pix,sidelength_pix); % % resp H and W in terms of pixels
                im = im.*bgval;
                [x,y] = find(im(:,:,1) == bgval);
                goodrows = find((abs(x-(mean(x))))<(currheight_pix/2) & ...
                    (abs(y-(mean(y))))<(sidelength_pix/2));
                for i = 1:numel(goodrows)
                    im(x(goodrows(i)),y(goodrows(i)),:)=currlum;
                end
                im2 = ones(size(im)).*placeholderval;
                im = cat(1,im,im2);
                tofill=size(im,1)-size(im,2);
                if (tofill>0)
                    imB1 = ones(size(im,1),round(tofill/2)).*bgval;
                    imB2 = ones(size(im,1),tofill-size(imB1,2)).*bgval;
                    im = cat(2,imB1,im,imB2);
                end
                im = imrotate(im,currrotangle);
                im(im==0)=bgval;


                %trim
                if(size(im,1)>maxdims(2))
                    totrim=size(im,1)-maxdims(2);
                    im=im(round(totrim/2)+1:end,:);
                    im=im(1:end-(totrim-round(totrim/2)),:);
                end
                if (size(im,2)>maxdims(1))
                    totrim=size(im,2)-maxdims(1);
                    im=im(:,round(totrim/2)+1:end);
                    im=im(:,1:end-(totrim-round(totrim/2)));
                end

                im0=im;clear im;

                % pad image to full screen size
                im0=padimage(im0,maxdims(2),maxdims(1));

                % make sure squares are still centered
                [x,y]=find(im0~=bgval);
                minrow=min(x); maxrow=max(x);
                mincol=min(y); maxcol=max(y);
                meanrow=mean([minrow maxrow]);
                meancol=mean([mincol maxcol]);
                if abs(meanrow-0.5*size(im0,1))>1, error('need to move'); end
                if abs(meancol-0.5*size(im0,2))>1, error('need to move'); end


                % separate the two possibilities for square location
                % base im
                im = im0; im(im0 == placeholderval)=bgval;
                imfin = makeimfin(im,currcolor,displayname,lums,colorlum);


                im_other = im0;
                im_other(im0 == currlum) = bgval;
                im_other(im0 == placeholderval) = currlum;
                imfin_other = makeimfin(im_other,currcolor,displayname,lums,colorlum);

                clear im0;

                colorstring=num2str(currcolor);
                colorstring=colorstring(~isspace(colorstring));

                fn = [shapename rotname lengthname figname colorlumname '_' colorstring];
                basefn = fn;

                temppos = strfind(fn,'S');
                fn_other = [fn(1:7) num2str(mod(currrotangle+180,360)) fn(temppos:end)];


                if numel(temppos)~=1
                    error('not found');
                else
                    imfin_swap=swapColors(imfin);

                    temppos=strfind(fn,'CL');
                    if strcmpi(fn(temppos-1),'G')
                        newfig='C';
                    elseif strcmpi(fn(temppos-1),'C')
                        newfig='G';
                    else
                        error('fig unclear');
                    end
                    swapfn=fn;
                    swapfn(temppos-1)=newfig;
                end

                %                 if sum(colors==[1 1 1])~=3
                % if grey+non-grey, generate alternative cases by swapping colors
                % 1) swap colors of base image
                % 3) swap color of other image
                temppos=strfind(fn_other,'CL');
                if strcmpi(fn_other(temppos-1),'G')
                    newfig = 'C';
                elseif strcmpi(fn_other(temppos-1),'C')
                    newfig = 'G';
                else
                    error('fig unclear');
                end
                imfin_other_swap=swapColors(imfin_other);
                swapfn_other=fn_other;
                swapfn_other(temppos-1)=newfig;

                %                 end

                currim = imfin;currfn=fn;
                imwrite(uint8(currim.*255),fullfile(basedatapath,[currfn '.bmp']));

                allscenes = cat(2,allscenes,{currfn});

                currim = imfin_other;currfn=fn_other;
                imwrite(uint8(currim.*255),fullfile(basedatapath,[currfn '.bmp']));

                allscenes =cat(2,allscenes,{currfn});

                currim = imfin_swap;currfn=swapfn;
                imwrite(uint8(currim.*255),fullfile(basedatapath,[currfn '.bmp']));

                allscenes = cat(2,allscenes,{currfn});

                currim = imfin_other_swap;currfn=swapfn_other;
                imwrite(uint8(currim.*255),fullfile(basedatapath,[currfn '.bmp']));

                allscenes = cat(2,allscenes,{currfn});

            end
        end
    end
end

%make an neutral background: color grey, with intermediate luminance. Seemed
%that this makes more sense compared to an 'intermediate color' as used in
%vdh
lumname = 'Lintermed'; %high luminance field is color
currlum = bglum;
tempcolor = [1 1 1];
imb = ones(maxdims(2),maxdims(1)).*currlum;
colormatrix=ones(size(imb,1),size(imb,2),3);
colormatrix(:,:,1)=colormatrix(:,:,1).*tempcolor(1);
colormatrix(:,:,2)=colormatrix(:,:,2).*tempcolor(2);
colormatrix(:,:,3)=colormatrix(:,:,3).*tempcolor(3);
if numel(unique(imb))==1 %we can simplify the lin_lum algorithm if it is an isoluminant field
    if size(unique(reshape(colormatrix,[size(colormatrix,1)*size(colormatrix,2),size(colormatrix,3),1]),'rows'),1)==1
        relintmatrix=lin_lum(imb(1),colormatrix(1,1,:),displayname);
        relintmatrix=repmat(relintmatrix,size(imb,1),size(imb,2));
    else
        relintmatrix=lin_lum(imb,colormatrix,displayname);
    end
else
    relintmatrix=lin_lum(imb,colormatrix,displayname);
end
imbfin(:,:,1)=colormatrix(:,:,1).*relintmatrix;
imbfin(:,:,2)=colormatrix(:,:,2).*relintmatrix;
imbfin(:,:,3)=colormatrix(:,:,3).*relintmatrix;
imwrite(imbfin,fullfile(basedatapath,['imb' lumname '_111.bmp']));

im=imread(fullfile(basedatapath, ['imb' lumname '_111.bmp']));
temp=unique([unique(im(1,1,:)) unique(im(end,end,:)) unique(im(1,end,:)) unique(im(end,1,:))]);
disp(['Color ' num2str(currcolor) ':'])
if numel(temp)>1 && ~(numel(temp)==2 && sum(temp==0)>0)
    error('non-unique background intensity');
elseif (numel(temp)==2 && sum(temp==0)>0)
    temp=temp(temp~=0);
else
    disp(['background should be ' num2str(im2double(temp))]);
    bgforBO=im2double(temp);
end

%% PROCESS ALLEN VIDEOS
% for i=1:n_allen_videos
%     if mod(i,2)==0
%         v=VideoReader('natural_movie_three.mp4'); %#ok<*TNMLP>
%     else
%         v=VideoReader('natural_movie_one.mp4');
%     end
%     startframei=randi(v.NumFrames-nframes);
%     video=read(v,[startframei startframei+nframes-1]);
%     v=VideoWriter(['video_allen_' num2str(i)],'MPEG-4');
%     v.Quality=100;
%     v.FrameRate=framerate;
%     open(v)
%     writeVideo(v,video);
%     close(v)
% end
%
% scenesize_allen=[v.Height v.Width];


%% MAKE BG FRAME FOR ALLEN
% bgfilename_allen = 'bgframe_allen';
% imt=ones(scenesize_allen(1),scenesize_allen(2),1).*bgluminance;
% relintmatrix=lin_lum(imt,ones(size(imt,1),size(imt,2),3),displayname);
% matrix(:,:,1)=relintmatrix;
% matrix(:,:,2)=relintmatrix;
% matrix(:,:,3)=relintmatrix;
% imwrite(matrix,[basedatapath bgfilename_allen '.bmp']);
% clear matrix;

%% Make Henaff et al. Stimuli
henaff_stim = load(fullfile(basedatapath,'Videostims_Henaffetal2021.mat'));
natural_movie_sequences = squeeze(henaff_stim.stim_matrix(:,1,:,:,:,:)); % Taking only the first luminance index in the 2nd dimension; image set : Zoom (2) x movie ID (10) x image pizels (x = 512) x image pixels (y = 512) x frames (11)
natural_movie_sequences_rearranged = permute(natural_movie_sequences,[1 2 5 3 4]); % Bringing the dimensions of varying conditions together; image set : Zoom (2) x movie ID (10) x frames (11) x image pizels (x = 512) x image pixels (y = 512)
dim_nms = size(natural_movie_sequences_rearranged);
num_allframes = dim_nms(1)*dim_nms(2)*dim_nms(3);

% Generating the frame IDs from the varying video frame conditions
frame_order = reshape(1:num_allframes,[dim_nms(1),dim_nms(2),dim_nms(3)]);
frame_pos = cell(1,num_allframes);
for i = 1:num_allframes
    [zoom,mov,frame] = ind2sub(size(frame_order),find(frame_order == i)); % finding the position of an element in the 3D matrix containing all video frames (zoom x movie ID x frames)
    frame_pos{1,i} = [zoom mov frame];
end

% Creating a 3D Matrix containing all frames
natural_movie_sequences_reshaped = reshape(natural_movie_sequences_rearranged,[num_allframes,dim_nms(end-1),dim_nms(end)]);

% saving each frame in basefolder/Henaff_stim % This needs to be on the
% same path as Monkeylogic timing file and condition file.
% dir_Henaff_stim = 'Henaff_stim';
% if ~exist(dir_Henaff_stim,'dir')
%     mkdir(fullfile(basedatapath,dir_Henaff_stim))
% end

bgLum_Henaff = process_HenaffImages(basedatapath,natural_movie_sequences_reshaped,'circular-crop',0,1,[],0);

for iframe = 1:size(natural_movie_sequences_reshaped,1)
    % temp_frame = squeeze(natural_movie_sequences_reshaped(iframe,:,:));
    % % temp_frame(temp_frame == 0) = bgforBO; % changing the background to match that of BO stimuli background; this approach is changing few pixels inside the image
    % imwrite(temp_frame, fullfile(basedatapath,['Henaff_stimID_' num2str(iframe) '.bmp']));
    allscenes = cat(2,allscenes,{['Henaff_stimID_' num2str(iframe) '.bmp']});
end

% background for Henaff Imageset
bgfilename_henaff='Henaff_bgLum'; % high luminance field is color
imt=ones(maxdims(2),maxdims(1),3).*bgLum_Henaff;
imwrite(imt,fullfile(basedatapath,['imb' bgfilename_henaff '.bmp']));

%% COLLECT STIMULI AND PUT IN RANDOM ORDER

% collect videos
% allvideos=[];
% for i=1:nvideos
%     d=dir([basedatapath 'video' num2str(i) '.avi']);
%     if numel(d)~=1
%         error('not found');
%     else
%         allvideos=[allvideos;{d(1).name} cell(1,4)]; %#ok<*AGROW> %empty cells are such that the cell array dimensions correspond to those for BO
%     end
% end
%
% % add Allen videos
% for i=1:n_allen_videos
%     d=dir([basedatapath 'video_allen_' num2str(i) '.mp4']);
%     if numel(d)~=1
%         error('not found');
%     else
%         allvideos=[allvideos;{d(1).name} cell(1,4)]; %empty cells are such that the cell array dimensions correspond to those for BO
%     end
% end

allscenes_forcond_all=[];
for ri=1:repetitions
    %parse BO scenes such that all stimuli are played once per trial
    allBOscenes_forcond=[];
    for pi=1:size(positions,1)
        allscenes_temp=allscenes(randperm(numel(allscenes)))';
        if mod(numel(allscenes),scenespertrial)~=0
            ntoadd=scenespertrial-mod(numel(allscenes),scenespertrial);
            allscenes_temp=[allscenes_temp; randsample(allscenes_temp,ntoadd)]; %#ok<*AGROW>
        end
        tempsize = size(allscenes_temp(reshape(1:numel(allscenes_temp),[],scenespertrial)));
        xpos=num2cell((ones(tempsize(1),1).*positions(pi,1))+RFpos_real(1));
        ypos=num2cell((ones(tempsize(1),1).*positions(pi,2))+RFpos_real(2));
        allBOscenes_forcond=[allBOscenes_forcond; allscenes_temp(reshape(1:numel(allscenes_temp),[],scenespertrial)) xpos ypos];
    end
    temp=allBOscenes_forcond;
    allscenes_forcond_all=[allscenes_forcond_all; temp(randperm(size(temp,1)),:)];
end

%% MAKE CONDITIONS FILE
%adjust position of fixation point such that RF is centered on square

savepath = basedatapath;

fileID = fopen(fullfile(savepath,'Henaff_stimAndBO.txt'),'w');

% fileID = fopen([savepath 'videosAndBO.txt'],'w');
% fileID = fopen(['C:\Users\Tom\Homework\' timingfilename '.txt'],'w');

clear TaskObject;
clear Info;

headerline = generate_condition('FID',fileID,'Header', 10);

currcondno=0;
allconditions=[];

for i=1:size(allscenes_forcond_all,1)

    currcondno=currcondno+1;

    if numel(strfind(allscenes_forcond_all{i,1},'squareR')) || numel(strfind(allscenes_forcond_all{i,1},'Henaff_stim'))>0
        %then BO stimulus trial
        currtimingfile = timingfilename_BO;

        currpos = [allscenes_forcond_all{i,end-1} allscenes_forcond_all{i,end}];

        fixbias=[0 0]-currpos;

        TOno=1;
        TaskObject(TOno).Type = 'pic'; %#ok<*SAGROW>
        TaskObject(TOno).Arg{1} = fixfilename;
        TaskObject(TOno).Arg{2} = fixbias(1);
        TaskObject(TOno).Arg{3} = fixbias(2);

        for scenei=1:scenespertrial
            TOno=TOno+1;
            TaskObject(TOno).Type = 'pic';
            TaskObject(TOno).Arg{1} = allscenes_forcond_all{i,scenei};
            TaskObject(TOno).Arg{2} = RFpos_real(1)+fixbias(1);
            TaskObject(TOno).Arg{3} = RFpos_real(1)+fixbias(1);
        end
        TOno=TOno+1;
        TaskObject(TOno).Type = 'pic';
        TaskObject(TOno).Arg{1} = ['imb' bgfilename_henaff];
        TaskObject(TOno).Arg{2} = 0;
        TaskObject(TOno).Arg{3} = 0;

        TOno=TOno+1;
        TaskObject(TOno).Type = 'sqr';
        TaskObject(TOno).Arg{1} = 100;
        TaskObject(TOno).Arg{2} = [1 1 1].*bgforBO;
        TaskObject(TOno).Arg{3} = 1;
        TaskObject(TOno).Arg{4} = 0;
        TaskObject(TOno).Arg{5} = 0;
        Info.info = 'BO';

    elseif numel(strfind(allscenes_forcond_all{i,1},'video'))>0
        %then video trial

        if numel(strfind(allscenes_forcond_all{i,1},'allen'))>0
            currbgfilename=bgfilename_allen;
        else
            currbgfilename=bgfilename;
        end

        fixbias=RF_relvideo-RFpos_real;

        currtimingfile = timingfilename_video;
        TOno=1;
        TaskObject(TOno).Type = 'pic';
        TaskObject(TOno).Arg{1} = fixfilename;
        TaskObject(TOno).Arg{2} = fixbias(1);
        TaskObject(TOno).Arg{3} = fixbias(2);
        TOno=TOno+1;
        TaskObject(TOno).Type = 'fix';
        TaskObject(TOno).Arg{1} = fixbias(1);
        TaskObject(TOno).Arg{2} = fixbias(2);
        TOno=TOno+1;
        TaskObject(TOno).Type = 'sqr'; %test stim to  check positioning
        TaskObject(TOno).Arg{1} = 0.2;
        TaskObject(TOno).Arg{2} = [1 1 1];
        TaskObject(TOno).Arg{3} = 1;
        TaskObject(TOno).Arg{4} = RFpos_real(1)+fixbias(1);
        TaskObject(TOno).Arg{5} = RFpos_real(2)+fixbias(2);

        currfn=allscenes_forcond_all{i,1};
        TOno=TOno+1;
        TaskObject(TOno).Type = 'mov';
        TaskObject(TOno).Arg{1} = currfn;
        TaskObject(TOno).Arg{2} = 0;
        TaskObject(TOno).Arg{3} = 0;
        TOno=TOno+1;
        TaskObject(TOno).Type = 'pic';
        TaskObject(TOno).Arg{1} = currbgfilename;
        TaskObject(TOno).Arg{2} = 0;
        TaskObject(TOno).Arg{3} = 0;
        TOno=TOno+1;
        TaskObject(TOno).Type = 'sqr';
        TaskObject(TOno).Arg{1} = 100;
        TaskObject(TOno).Arg{2} = [0 0 0];
        TaskObject(TOno).Arg{3} = 1;
        TaskObject(TOno).Arg{4} = 0;
        TaskObject(TOno).Arg{5} = 0;
        Info.info = currfn;
    end


    textline = generate_condition('FID',fileID,'Condition', currcondno, ...
        'Frequency', 1, 'Block',1,'TimingFile', currtimingfile, 'Info', Info,...
        'TaskObject', TaskObject);
    clear TaskObject;
    allconditions=[allconditions;...
        currcondno NaN];

end
fclose(fileID);

%% SAVE STIM PARAMS
save stimparams_BO stimdur scenespertrial blankdur;

%% GENERATE CONDITION PLAY ORDER
condplayorder=[];
for i=1:size(allscenes_forcond_all,1)
    condplayorder=[condplayorder;i];
end
save condplayorder condplayorder;