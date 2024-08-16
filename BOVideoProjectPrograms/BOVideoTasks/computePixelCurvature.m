% This function computes the discrete pixel curvatures at different time
% frames of a video and global pixel curvatures across frames in an
% image trajectory. Run this code with appropriate/default arguments
% with the Henaff Data/Data folder in the path and current folder set 
% to the path of this program.


function [ct_pixel, globalCurvature_imageSet, ct_pixel_pca, globalCurvature_pca_imageSet] = computePixelCurvature(data,imageType,imagematrix,frameFormat,numPC)

if ~exist('data','var'); data = 'HenaffStim2021'; end
if ~exist('imageType', 'var'); imageType = 'artificial'; end
if ~exist('imagematrix', 'var'); imagematrix = 'HQ'; end
if ~exist('frameFormat', 'var'); frameFormat = 'All Frames'; end
if ~exist('numPC', 'var'); numPC = 2; end

if strcmp(data,'HenaffStim2021')

    tic;
    %% load data
    fileNameWithPath = mfilename('fullpath');
    [filePath,~,~] = fileparts(fileNameWithPath);
    cd(filePath);
    addpath(genpath(pwd));

    HenaffImageData = load(fullfile(filePath,'Videostims_Henaffetal2021.mat')); % This program assumes that  the image data is located in the same path as the function!

    %% Natural Movie frames
    naturalImageLabels = HenaffImageData.natural_movie_labels;
    if strcmp(imagematrix,'HQ')
        naturalImageSet = squeeze(HenaffImageData.stim_matrix(:,1,:,:,:,:)); % setting dimension 2 index as 1 selects the natural images
    elseif strcmp(imagematrix,'blurred') % blurred version has pixalated images with slighly poorer image resolution
        naturalImageSet = squeeze(HenaffImageData.stim_matrix_blurred(:,1,:,:,:,:)); % setting dimension 2 index as 1 selects the natural images
    end

    % Arificial/Synthetic Movie frames
    artificalImageLabels = HenaffImageData.artificial_movie_labels;
    if strcmp(imagematrix,'HQ')
        artificalImageSet = squeeze(HenaffImageData.stim_matrix(:,2,:,:,:,:)); % setting dimension 2 index as 2 selects the synthetic images
    elseif strcmp(imagematrix,'blurred')
        artificalImageSet = squeeze(HenaffImageData.stim_matrix_blurred(:,2,:,:,:,:));
    end

    if strcmp(imageType,'natural')
        imageLabels = naturalImageLabels;
        switch(frameFormat)
            case 'First 6-Frame'
                frameNums = 1:6;
            case 'Alternate Frames'
                frameNums = 1:2:11;
            case 'All Frames'
                frameNums = 1:11;
        end
        imageSet = naturalImageSet(:,:,:,:,frameNums);


    elseif strcmp(imageType,'artificial')
        imageLabels = artificalImageLabels;
        switch(frameFormat)
            case 'First 6-Frame'
                frameNums = 1:6;
            case 'Alternate Frames'
                frameNums = 1:2:11;
            case 'All Frames'
                frameNums = 1:11;
        end
        imageSet = artificalImageSet(:,:,:,:,frameNums);
    end


    %% Arrange Image Sequences by Sequence ID (ID x Zoom (10 x 2) x pixel
    % (512) x pixel (512) x frames (11) to make it similar to the order in
    % Image Labels variable
    count = 1;
    arraySize = size(imageSet);
    sortedImageSet = zeros([length(imageLabels) arraySize(3:5)]);
    for iZoom = 1:size(imageSet,1) % Zoom Condition (Spatial Scale)
        for iImageSequence = 1:size(imageSet,2) % Video Frame
            sortedImageSet(count,:,:,:) = imageSet(iZoom,iImageSequence,:,:,:);
            count = count + 1;
        end
    end

    %% create and save gif file of images for quick Reference
    for iVideo = 1:size(sortedImageSet,1) % image
        if strcmp(imageType,'natural')
            folderSave =  fullfile(filePath,'Stimuli','Natural_Video-Henaff_2021');
        elseif strcmp(imageType,'artificial')
            folderSave =  fullfile(filePath,'Stimuli','Artificial_Video-Henaff_2021');
        end
        if ~exist(folderSave,'dir')
            mkdir(folderSave)
        end
        fileName = fullfile(folderSave,[imageLabels{iVideo} '.gif']);
        if ~exist('fileName','var')
        else
            for iFrame = 1:size(sortedImageSet,4) % Video Frames
                imageFrame = uint8(squeeze(sortedImageSet(iVideo,:,:,iFrame)).*255);
                if iFrame == 1
                    imwrite(imageFrame,fileName,"gif","LoopCount",Inf,"DelayTime",0.1);
                else
                    imwrite(imageFrame,fileName,"gif","WriteMode","append","DelayTime",0.1);
                end
            end
        end
    end

    %% Calculating the discrete curvature at diffferent time frames

    % reshape 2-D x-y pixel data into 1D array (N-length vector)
    numVideos = size(sortedImageSet,1);
    pixel_vectorDimension = size(sortedImageSet,2)*size(sortedImageSet,3);
    numFrames = size(sortedImageSet,4);
    sortedImageSet_vectorFormat = reshape(sortedImageSet,numVideos,pixel_vectorDimension,numFrames);

    for iVideo = 1:size(sortedImageSet,1) % image ID
        clear temp_Video delta_xt_pixel temp_vt_pixel
        temp_Video =  squeeze(sortedImageSet_vectorFormat(iVideo,:,:)); % selects one video at a time
        delta_xt_pixel = diff(temp_Video,1,2); % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
        norm_vt_pixel = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel,1),'UniformOutput',false); % computing the unit displacement vectors
        moving_dotProduct_norm_vt_pixel = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel(1:end-1),norm_vt_pixel(2:end),'UniformOutput',false)); % dot product of sequential displacement vectors
        ct_pixel(iVideo,:) = rad2deg(acos(moving_dotProduct_norm_vt_pixel)); %#ok<*AGROW>
    end

    % Calculate global curvature for each image sequence/video
    globalCurvature_imageSet = mean(ct_pixel,2);
    toc;

    % Displaying discrete curvature and compute global curvature

    disp(['Summary of Hénaff Imageset: ' upper(imageType) ' imageset, ' upper(imagematrix) ', ' frameFormat])
    tab = table(imageLabels',ct_pixel,globalCurvature_imageSet);
    tab.Properties.VariableNames = {'Image Label', 'Discrere Curvatures', 'Global Curvature'};
    disp(tab)

    % %% PCA of digital images to compress them using n Principal Components
    % % and calculating discrete and global curvature of compressed images
    % warning('off')
    % for iImage = 1:size(sortedImageSet,1)
    %     for iFrame = 1:size(sortedImageSet,4)
    %         temp_image = squeeze(sortedImageSet(iImage,:,:,iFrame));
    %         [pc,score,~,~,~,mu] = pca(temp_image,'NumComponents',numPC); % Principal Component Analysis on each image
    %         temp_pca_image_nPC = score(:,1:numPC)*pc(:,1:numPC)' + mu; % Image Compression using the first 2 Principal Components
    %         pca_imageset_nPC(iImage,:,iFrame) = temp_pca_image_nPC(:); % Collapsing the 2 dimensional image vector into 1D
    %     end
    % end
    % warning('on')
    %
    % % Calculating discrete curvature on the PC-constrcuted images
    % for iImage = 1:size(sortedImageSet,1) % image ID
    %     clear temp_image delta_xt_pixel norm_vt_pixel moving_dotProduct_norm_vt_pixel
    %     temp_image =  squeeze(pca_imageset_nPC(iImage,:,:)); % selects one image at a time
    %     delta_xt_pixel = diff(temp_image,1,2); % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
    %     norm_vt_pixel = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel,1),'UniformOutput',false);
    %     moving_dotProduct_norm_vt_pixel = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel(1:end-1),norm_vt_pixel(2:end),'UniformOutput',false));
    %     ct_pixel_pca(iImage,:) = rad2deg(acos(moving_dotProduct_norm_vt_pixel)); %#ok<*AGROW>
    % end
    %
    % globalCurvature_pca_imageSet = mean(ct_pixel_pca,2);
    %
    % % Displaying discrete curvature and compute global curvature
    %
    % disp(['Summary of 5-PC Hénaff Imageset: ' upper(imageType) ' imageset, ' upper(imagematrix) ', ' frameFormat])
    % tab = table(imageLabels',ct_pixel_pca,globalCurvature_pca_imageSet);
    % tab.Properties.VariableNames = {'Image Label', 'Discrere Curvatures', 'Global Curvature'};
    % disp(tab)

    %% PCA of Video Trajectory
    frameSize = size(sortedImageSet,2,3);
    frameNum = size(sortedImageSet,4);
    numPCList = [frameNum, numPC];
    for iVideo = 1:size(sortedImageSet,1)
        for iPCList = 1:length(numPCList)
            temp_video = reshape(squeeze(sortedImageSet(iVideo,:,:,:)),[prod(frameSize) frameNum]); % Collapsing the 2 dimensional image vector into 1D for each frame
            [pc,score,~,~,~,mu] = pca(temp_video,'NumComponents',numPCList(iPCList)); % Principal Component Analysis on each Video across all frames (Video Trajectory)
            reconstructedVideo = score(:,1:numPCList(iPCList))*pc(:,1:numPCList(iPCList))' + mu; % Video reconstruction using the first n Principal Components (dimensionality reduction)
            delta_xt_pixel = diff(reconstructedVideo,1,2); % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
            norm_vt_pixel = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel,1),'UniformOutput',false);
            moving_dotProduct_norm_vt_pixel = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel(1:end-1),norm_vt_pixel(2:end),'UniformOutput',false));
            ct_pixel_pca(iVideo,iPCList,:) = rad2deg(acos(moving_dotProduct_norm_vt_pixel)); %#ok<*AGROW>
            globalCurvature_pca_imageSet(iVideo,iPCList) = mean(squeeze(ct_pixel_pca(iVideo,iPCList,:)));
        end
    end

elseif isempty(data)||~exist('data','var')||strcmp(data,'')
    error('Please specify data specifics; data not found!')
else
    error('Please specify data specifics; data not found!')
end
end

