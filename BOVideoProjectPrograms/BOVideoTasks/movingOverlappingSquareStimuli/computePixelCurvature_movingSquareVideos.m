% This function takes all the frames of the moving sqaure stimuli and
% generates the pixel curvature for a central circular portion of the
% image (which should supposedly overlap with the RF of recorded site)

function [ct_pixel,globalCurvature_CroppedVideo] = computePixelCurvature_movingSquareVideos(allframes,processingMode,manualCroppingMode,imageID,demoMode)

if ~exist('imageID','var');   imageID = 1; end
if ~exist('processingMode','var');   processingMode = 'circular-crop'; end
if ~exist('manualCroppingMode','var');   manualCroppingMode = 0; end
% if ~exist('bgforBO','var');   bgforBO = []; end
if ~exist('demoMode','var');   demoMode = 0; end

scene_size_allFrames = cellfun(@size,allframes,'UniformOutput',false);
if all(cellfun(@(x) isequal(x, scene_size_allFrames{1}), scene_size_allFrames))
    scene_size = scene_size_allFrames{1};
else
    error('Video Frame diemnsion Mismatch!')
end

numFrames = length(allframes);

if strcmp(processingMode,'original') % No processing done on the original set

    % Full frames will be processed unaltered!

elseif strcmp(processingMode,'circular-crop') % A circular region will be cropped for all frames for further curvature computation!
    % Get selected image
    selectedOriginalImage = allframes{imageID};
    cx = round(scene_size(2)/2); % x-coordinate of the circle center
    cy = round(scene_size(1)/2); % y-coordinate of the circle center


    if manualCroppingMode
        % Maximize the window to make it easier to draw.
        g = gcf;
        g.WindowState = 'maximized';
        % Ask user to draw a circle:
        uiwait(helpdlg('Please click and drag out a circle.'));
        h.Radius = 0;
        while h.Radius == 0
            h = drawcircle('Color','r','FaceAlpha',0.4);
            if h.Radius == 0
                uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
            end
        end
    else
        % manual circle handle
        h.Radius = 315; % radius = 234.7734 from drawn circle in Image 1
        h.Center = [cy cx]; %center = [259.5296 258.3935] from drawn circle in Image 1
    end

    % Get coordinates of the circle.
    angles = linspace(0, 2*pi, 10000);
    x = cos(angles) * h.Radius + h.Center(2);
    y = sin(angles) * h.Radius + h.Center(1);

    % Image dimensions
    [rows, columns, numberOfColorChannels] = size(selectedOriginalImage);

    % process all frames in a video with the set circle radius and center
    % Get a mask of the circle - Mask Region is identical across all
    % images in the set
    mask = poly2mask(x, y, rows, columns);

    if demoMode

        % Mask the image with the circle.
        if numberOfColorChannels == 1
            maskedImage = selectedOriginalImage; % Initialize with the entire image.
            maskedImage(mask==0) = 1; % Zero image outside the circle mask.
        else
            % Mask the image. % Processes images with color channels
            maskedImage = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
        end

        % Demo to have the user click and draw a circle over an image,
        % then blacken outside the circle and crop out the circular
        % portion into a new image.
        fontSize = 14;
        figure;
        title('Image masked with the circle', 'FontSize', fontSize);
        % fprintf('Done running %s.m ...\n', mfilename);
        subplot(2, 2, 1);
        imshow(selectedOriginalImage./255);
        axis('on', 'image');
        title('Original Image', 'FontSize', fontSize);

        % Show circle over image.
        subplot(2, 2, 2);
        imshow(selectedOriginalImage./255);
        axis('on', 'image');
        hold on;
        plot(x, y, 'r-', 'LineWidth', 2);
        title('Original image with circle mask overlaid', 'FontSize', fontSize);

        subplot(2, 2, 3);
        imshow(mask./255);
        axis('on', 'image');
        title('Circle Mask', 'FontSize', fontSize);

        % Crop the image to the bounding box.
        % props = regionprops(mask, 'BoundingBox');
        % maskedImage = imcrop(maskedImage, props.BoundingBox);
        % Display it in the lower right plot.
        subplot(2, 2, 4);
        imshow(maskedImage./255, []);
    end



    for iFrame = 1:numFrames
        tempFrame = allframes{iFrame};
        [rows, columns, numberOfColorChannels] = size(tempFrame);

        if rows ~= scene_size(1) || columns ~= scene_size(2)
            error (['Image dimension of image ID: ' num2str(iFrame) 'is not ' num2str(scene_size(1)) 'x' num2str(scene_size(2)) '! Check Image dataset!'])
        else

            % Mask the image with the circle.
            if numberOfColorChannels == 1
                % tempMaskedFrame = tempFrame; % Initialize with the entire image.
                % tempMaskedFrame(mask==0) = 1; % Zero image outside the circle mask.
                croppedFrames(:,iFrame) = tempFrame(mask~=0)./255;
            else
                % Mask the image. % Processes images with color
                % channels % needs to be worked on!
                tempMaskedFrame = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
                tempMaskedFrame(mask==0) = 1; % Zero image outside the circle mask.
                croppedFrames(:,iFrame) = tempMaskedFrame(mask~=0)./255;
            end
        end
    end

    [ct_pixel,globalCurvature_CroppedVideo] = computePixelCurvature_croppedFrame(croppedFrames);

    % grandMeanPixelLuminance = mean(luminanceList_CircleImageSet);
    % bgforHenaffStimSet = grandMeanPixelLuminance;
    %
    % for iImage = 1:numFrames
    %     tempFrame = allframes{iFrame};
    %
    %     % Mask the image with the circle.
    %     if numberOfColorChannels == 1
    %         tempMaskedFrame = tempFrame; % Initialize with the entire image.
    %         if ~isempty(bgforBO)
    %             if iImage == 1
    %                 display(['bg for Henaff Stim set is set as bg for BO stim:' num2str(bgforBO)])
    %             end
    %             tempMaskedFrame(mask==0) = bgforBO; % Zero image outside the circle mask.
    %
    %         else
    %             if iImage == 1
    %                 display(['bg for Henaff Stim set is set as the grand avg of the circular cropped portion of all images in the set:' num2str(bgforHenaffStimSet)])
    %             end
    %             tempMaskedFrame(mask==0) = bgforHenaffStimSet; % Zero image outside the circle mask.
    %         end
    %     else
    %         % Mask the image. % Processes images with color channels
    %         tempMaskedFrame = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
    %     end
    %     imwrite(tempMaskedFrame, fullfile(folderSourceString,['Henaff_stimID_' num2str(iImage) '.bmp']));
    % end

end
end



%% Accessory function

function [ct_pixel,globalCurvature_CroppedVideo] = computePixelCurvature_croppedFrame(croppedFrames)
% reshape 2-D x-y pixel data into 1D array (N-length vector)
% numVideos = size(sortedImageSet,1);
% pixel_vectorDimension = size(croppedFrames,2);
% numFrames = size(croppedFrames,1);
% sortedImageSet_vectorFormat = reshape(sortedImageSet,numVideos,pixel_vectorDimension,numFrames);

% for iVideo = 1:size(sortedImageSet,1) % image ID
% clear temp_Video delta_xt_pixel temp_vt_pixel
temp_Video =  croppedFrames; % selects one video at a time
delta_xt_pixel = diff(temp_Video,1,2); % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
norm_vt_pixel = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel,1),'UniformOutput',false); % computing the unit displacement vectors
moving_dotProduct_norm_vt_pixel = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel(1:end-1),norm_vt_pixel(2:end),'UniformOutput',false)); % dot product of sequential displacement vectors
ct_pixel = rad2deg(acos(moving_dotProduct_norm_vt_pixel)); %#ok<*AGROW>
% end

% Calculate global curvature for each image sequence/video
globalCurvature_CroppedVideo = mean(ct_pixel,2);

% Displaying discrete curvature and compute global curvature

% disp(['Summary of Hénaff Imageset: ' upper(imageType) ' imageset, ' upper(imagematrix) ', ' frameFormat])
% tab = table(imageLabels',ct_pixel,globalCurvature_imageSet);
% tab.Properties.VariableNames = {'Image Label', 'Discrere Curvatures', 'Global Curvature'};
% disp(tab)
end