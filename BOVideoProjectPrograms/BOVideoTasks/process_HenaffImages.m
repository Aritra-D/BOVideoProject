% Process the Henaff Images- Crop Only the central portion of Images and
%  make the background identical to BO stimuli

% Manual Cropping Mode helps you choose an image to draw a circle and use
% the radius and center of the drawn circle to crop the image

% demoMode shows you how the cropped image will look like
% setting demoMode to zero processes the

function bgforHenaffStimSet = process_HenaffImages(folderSourceString,imageSet,processingMode,manualCroppingMode,imageID,bgforBO,demoMode)

if ~exist('folderSourceString','var')|| isempty(folderSourceString)
    folderSourceString = 'C:\Users\Aritra\OneDrive - Washington University in St. Louis\Lab Workbench\FrankenLab-WashU\Projects\NaturalVideo-BO\Programs\BOvideo_scripts';
end
if ~exist('imageID','var');   imageID = 1; end
if ~exist('imageSet','var')|| isempty(imageSet)
temp = load(fullfile(folderSourceString,'imageData.mat'));
imageSet = temp.natural_movie_sequences_reshaped;
end

numImages = size(imageSet,1);

if strcmp(processingMode,'original') % No processing done on the original set
    for iImage = 1: numImages
        tempImage = squeeze(imageSet(iImage,:,:));
        imwrite(tempImage, fullfile(folderSourceString,['Henaff_stimID_' num2str(iImage) '.bmp']));
    end

elseif strcmp(processingMode,'circular-crop')
    % Get selected image
    selectedOriginalImage = squeeze(imageSet(imageID,:,:));

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
        h.Radius = 234; % radius = 234.7734 from drawn circle in Image 1
        h.Center = size(selectedOriginalImage)/2; %center = [259.5296 258.3935] from drawn circle in Image 1
    end

    % Get coordinates of the circle.
    angles = linspace(0, 2*pi, 10000);
    x = cos(angles) * h.Radius + h.Center(1);
    y = sin(angles) * h.Radius + h.Center(2);

    % Image dimensions
    [rows, columns, numberOfColorChannels] = size(selectedOriginalImage);

    % Get a mask of the circle
    mask = poly2mask(x, y, rows, columns);

    if demoMode

    % Mask the image with the circle.
    if numberOfColorChannels == 1
        maskedImage = selectedOriginalImage; % Initialize with the entire image.
        maskedImage(mask==0) = bgforBO; % Zero image outside the circle mask.
    else
        % Mask the image. % Processes images with color channels
        maskedImage = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
    end

    % Demo to have the user click and draw a circle over an image,
    % then blacken outside the circle and crop out the circular
    % portion into a new image.

    title('Image masked with the circle', 'FontSize', fontSize);
    % fprintf('Done running %s.m ...\n', mfilename);
    subplot(2, 2, 1);
    imshow(selectedOriginalImage);
    axis('on', 'image');
    title('Original Image', 'FontSize', fontSize);

    % Show circle over image.
    subplot(2, 2, 2);
    imshow(selectedOriginalImage);
    axis('on', 'image');
    hold on;
    plot(x, y, 'r-', 'LineWidth', 2);
    title('Original image with circle mask overlaid', 'FontSize', fontSize);

    subplot(2, 2, 3);
    imshow(mask);
    axis('on', 'image');
    title('Circle Mask', 'FontSize', fontSize);

    % Crop the image to the bounding box.
    % props = regionprops(mask, 'BoundingBox');
    % maskedImage = imcrop(maskedImage, props.BoundingBox);
    % Display it in the lower right plot.
    subplot(2, 2, 4);
    imshow(maskedImage, []);

    else % process all images in the set with the set circle radius and center
        % Get a mask of the circle - Mask Region is identical across all
        % images in the set
        mask = poly2mask(x, y, rows, columns);

        for iImage = 1:numImages
            tempImage = squeeze(imageSet(iImage,:,:));
            [rows, columns, numberOfColorChannels] = size(tempImage);

            if rows ~= 512 || columns ~= 512
                error (['Image dimension of image ID: ' num2str(iImage) 'is not 512 x 512! Check Image dataset!'])
            else
 
                % Mask the image with the circle.
                if numberOfColorChannels == 1
                    tempMaskedImage = tempImage; % Initialize with the entire image.
                    % tempMaskedImage(mask==0) = bgforBO; % Zero image outside the circle mask.
                else
                    % Mask the image. % Processes images with color channels
                    tempMaskedImage = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
                end
                luminanceList_CircleImageSet(iImage) = mean(tempMaskedImage(mask~=0)); %#ok<*AGROW>
            end
        end

        grandMeanPixelLuminance = mean(luminanceList_CircleImageSet);
        bgforHenaffStimSet = grandMeanPixelLuminance;

        for iImage = 1:numImages
            tempImage = squeeze(imageSet(iImage,:,:));
            
            % Mask the image with the circle.
            if numberOfColorChannels == 1
                tempMaskedImage = tempImage; % Initialize with the entire image.
                if ~isempty(bgforBO)
                    if iImage == 1
                    display(['bg for Henaff Stim set is set as bg for BO stim:' num2str(bgforBO)])
                    end
                    tempMaskedImage(mask==0) = bgforBO; % Zero image outside the circle mask.
        
                else
                    if iImage == 1
                    display(['bg for Henaff Stim set is set as the grand avg of the circular cropped portion of all images in the set:' num2str(bgforHenaffStimSet)])
                    end
                    tempMaskedImage(mask==0) = bgforHenaffStimSet; % Zero image outside the circle mask.
                end
            else
                % Mask the image. % Processes images with color channels
                tempMaskedImage = bsxfun(@times, selectedOriginalImage, cast(mask, class(selectedOriginalImage)));
            end
            imwrite(tempMaskedImage, fullfile(folderSourceString,['Henaff_stimID_' num2str(iImage) '.bmp']));
        end

    end
end
end