% PCA tutorial 
% https://www.mathworks.com/matlabcentral/answers/881862-problem-understanding-pca-and-eigenvectors-of-covariance-matrix

% x = rand(10,3);
% 
% n = size(x,1);
% mu = mean(x);
% xcm = x - mu;               % cm values for each column
% covx = (xcm'*xcm)/(n-1);    % covariance matrix
% [v, lambda] = eig(covx);
% % for consistency with pca, sort the eigenvalues and permute
% % the columns of the eigenvector matrix to match.
% [lambda, ind] = sort(diag(lambda),'descend');
% v = v(:,ind);
% 
% coeff = v;           % each column of coeff describes a principal vector
% score = xcm*v;      % scores are projections of each observation onto the principal vectors
% latent = lambda;     % magnitude of each principal vector
% 
% % compare with pca
% [coeff1,score1,latent1] = pca(x);
% % the principal vectors can differ by a factor of -1 between methods, so
% % the coeff ratio below may have either +1 or -1 down columns. 
% % However, the score ratio below will have matching -1 down its columns, so the desription
% % of observations in terms of principal vectors is unchanged.
% % latent values, being eigenvalues, always match.
% coeffComparison = coeff./coeff1;
% scoreComparison = score./score1;
% latentComparison = latent./latent1;


% PCA to compress image to n Principal Components
% filePath = 'C:\Users\Aritra\OneDrive - Washington University in St. Louis\Lab Workbench\FrankenLab-WashU\Projects\BOVideoProject\BOVideoProjectPrograms\BOVideoTasks';
% HenaffImageData = load(fullfile(filePath,'Videostims_Henaffetal2021.mat'));
% 
% % step by step PCA procedure 
% exampleImage = squeeze(HenaffImageData.stim_matrix(1,1,2,:,:,1));
% x = exampleImage;
% n = size(x,1);
% mu = mean(x,1);
% xcm = x - mu;
% covx = (xcm'*xcm)/(n-1);
% [v, lambda] = eig(covx);
% [lambda, ind] = sort(diag(lambda),'descend');
% v = v(:,ind);
% coeff = v;           % each column of coeff describes a principal vector
% score = xcm*v;      % scores are projections of each observation onto the principal vectors
% latent = lambda;
% reconstructedImage = (score*coeff')+mu;
% 
% figure(1);imshow(exampleImage)
% figure(2);imshow(reconstructedImage)
% 
% % PCA using the MATLAB inbuilt function
% [coeff1,score1,latent1,tsquared1,explained1,mu1] = pca(x);
% reconstructedImage2 = score1*coeff1'+ mu1;
% figure(3);imshow(reconstructedImage2)


% PCA on the video tracjectory
clear; close all;
filePath = 'C:\Users\Aritra\OneDrive - Washington University in St. Louis\Lab Workbench\FrankenLab-WashU\Projects\BOVideoProject\BOVideoProjectPrograms\BOVideoTasks';
HenaffImageData = load(fullfile(filePath,'Videostims_Henaffetal2021.mat'));


frameSize = size(HenaffImageData.stim_matrix,4,5);
frameNum = 6; %size(HenaffImageData.stim_matrix,6);


numPCList = 2;
example_video = reshape(squeeze(HenaffImageData.stim_matrix(1,2,1,:,:,1:2:11)),[prod(frameSize) frameNum]);
ct_pixel_pca = zeros(length(numPCList),frameNum-2);
globalCurvature_pca_imageSet =  zeros(1,length(numPCList));

for iPC = 1:length(numPCList)
    [pca_vt,score_vt,~,~,exVar_vt,mu_vt] = pca(example_video,'NumComponents',numPCList(iPC));
    temp_reconstructedVideo = score_vt*pca_vt'+mu_vt;
    
    reconstructedVideo = reshape(temp_reconstructedVideo,frameSize(1),frameSize(2),[]);
    implay(reconstructedVideo)


    

    % Calculating discrete curvature on the PC-constrcuted images
    delta_xt_pixel = diff(temp_reconstructedVideo,1,2); % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
    norm_vt_pixel = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel,1),'UniformOutput',false);
    moving_dotProduct_norm_vt_pixel = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel(1:end-1),norm_vt_pixel(2:end),'UniformOutput',false));
    ct_pixel_pca(iPC,:) = rad2deg(acos(moving_dotProduct_norm_vt_pixel)); %#ok<*AGROW>
    globalCurvature_pca_imageSet(iPC) = mean(ct_pixel_pca(iPC,:),2);
end

scatter(pca_vt(:,1),pca_vt(:,2),'filled')
plot(pca_vt(:,1),pca_vt(:,2),'-o'); 

delta_xt_pixel2 = diff(pca_vt,1)'; % difference in vectors of pixel intensities sequence of vectors representing sequential time frames; note reduction of frame dimension from n to n-1
norm_vt_pixel2 = cellfun(@(x) x./norm(x), num2cell(delta_xt_pixel2,1),'UniformOutput',false);
moving_dotProduct_norm_vt_pixel2 = cell2mat(cellfun(@(x,y) dot(x,y), norm_vt_pixel2(1:end-1),norm_vt_pixel2(2:end),'UniformOutput',false));
ct_pixel_pca2 = rad2deg(acos(moving_dotProduct_norm_vt_pixel2)); %#ok<*AGROW>
globalCurvature_pca_imageSet2 = mean(ct_pixel_pca2);
