function [goodStack,goodFlags,training] = tom_os3_classifyImageStack(imageStack,dimension,numberOfEigenImages,numberOfCentroidImages,mask,normflag,sigma,trainingReference,alignReference,desperateFlag,numberOfClusters)
%tom_os3_classifyImageStack
%
%
%
%   tom_os3_classifyImageStack(imageStack,dimension,numberOfEigenImages,numberOfCentroidImages,numberOfClusters,mask,normflag,sigma,trainingReference,alignReference)
%
%PARAMETERS
%
%  INPUT
%       imageStack - each layer must be normed by 'phase' with tom_norm
%       dimension - the dimension of the particles (2 or 3 d)
%       numberOfEigenImages - the number of eigenimages(vectors) used for
%       reduction
%
%       numberOfCentroidImages - number of images contributing to centroid
%       calculation
%
%       mask - a binary mask used to cut out the gegion of interest
%       normflag - an additional normalizetion of each stack layer
%       sigma - the deviation from the centroid
%       trainingReference - structure containing PCA Eigenvectors and data
%       alignReference - aligment reference for each stack layer
%       desperateFlag - use additional aligment information for
%       classification
%       numberOfClusters - numberOfClusters
%lowest frequ (in pixels)
%  OUTPUT
%       goodStack
%       goodFlags
%       training
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by TH2 07/07/07
%   updated by ..
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if(ischar(imageStack))

    imageStack = tom_emread(imageStack);
    imageStack = imageStack.Value;

end;
% if(nargin <2)
%     error('Please set the dimensio of the stack.');
% end;

if(nargin < 3)
    numberOfEigenImages = 5;
end;

if(nargin < 4)
    numberOfCentroidImages =8;
    %hier kann man noch mit adaptiver bestimmung den wert festlegen.
    %ab wann fangen die ersten werte an zu varieren, da setze den
    %schwellwert
    mask = tom_os3_sphereMask(1);
    normflag = 1;
    trainingReference = [];
    alignReference = [];
    sigma = 1;
end;


if(nargin < 5)
    mask = tom_os3_sphereMask(1);
    normflag = 0;
    trainingReference = [];
    alignReference = [];
    sigma = 1;
end;

if(nargin < 6)
    normflag = 0;
    trainingReference = [];
    alignReference = [];
    sigma = 1;
end;

if(nargin < 7)
    trainingReference = [];
    alignReference = [];
    sigma = 1;
end;
if(nargin < 8)
    trainingReference = [];
    alignReference = [];

end;
if(nargin < 9)
    alignReference = [];
end;

if(nargin <10)
    desperateFlag = 0;
end;

if(nargin < 11)
    if(~isempty(trainingReference) && isfield(trainingReference,'numberOfClusters'))
        numberOfClusters = trainingReference.numberOfClusters;
    else
        numberOfClusters = 1;
    end;
end;
originalStack = imageStack;
stackSize = size(imageStack);
goodFlags = ones([stackSize(3),1],'single');
%%  align the particle stack if align reference is given
if(~isempty(alignReference))
    stck = [];
    stckCounter = 1;
    reference = alignReference;
    filter_param.Apply = 0;

    reference = tom_norm((reference+100)*100,'phase');

    if(exist('trainingReference') && ~isempty(trainingReference))
        sphere = trainingReference.alignTranslationMask;
        transMask = trainingReference.alignRotationMask;
    else
        sphere = tom_os3_sphereMask(reference);
        transMask = tom_os3_pasteCenter(zeros(size(reference),'single'),tom_os3_sphereMask(zeros(20)));
    end;

    iterations = 1;
    %     figure;

    if(desperateFlag == 1)
        sphere = ones(size(reference));
        transMask = ones(size(reference));
    end;

    for i= 1:size(imageStack,3)

        [a b c d]=tom_av2_align(reference,imageStack(:,:,i),ones(size(reference)),[],transMask,filter_param,iterations,0);


        if(desperateFlag == 1)
            if(((abs(a(1)) > 10 && abs(a(1)) < 170) || (abs(a(1)) > 190 && abs(a(1)) < 350)) || (abs(b(1)) > 1 || abs(b(2)) > 1))
                goodFlags(i) = 0;
            end;

            %             if(goodFlags(i) == 0)
            %                 tom_imagesc(imageStack(:,:,i));
            %             end;
        end;
      %  imageStack(:,:,i)= tom_norm(imageStack(:,:,i),'phase');
        imageStack(:,:,i)= tom_rotate(imageStack(:,:,i),a,'linear');% .* mask;

    end;

end;

%allocate memory for reshaped stack
reshapedStack = zeros([stackSize(end),prod(stackSize(1:end-1))],'single');



%% reshape image stack
if(dimension == 2)

    for i=1:stackSize(3)


        im = imageStack(:,:,i);
        im = reshape(im,1,[]);

        if (normflag == 1)
            im = tom_norm(im+100,'phase');
        end
        reshapedStack(i,:)=im;
    end;

else
    %3d


end;

%%  execute pca analysis if no eigenvectors are given
if(isempty(trainingReference))
    [scores,coefs,eigenValues]=tom_calc_pca(double(reshapedStack),numberOfEigenImages);
    % determine the center in reduced space of the first n images

%     [classification clusterPosition sumDistance clusterDistances] = kmeans(scores',numberOfClusters);
% 
%     avgDistances = zeros(numberOfClusters,2);
%     for i =1:numberOfClusters
% 
%         dist = [];
%         for j =1:numel(classification)
% 
%             if(classification(j) == i)
%                 dist(numel(dist)+1) = clusterDistances(j,i);
%             end;
%         end;
% 
%         avgDistances(i,1) = mean(dist);
%         avgDistances(i,2) = std(dist);
%     end;



        % calculate the centroid of each transformed image
        center = zeros(numberOfEigenImages,1,'single');
        for i=1:numberOfEigenImages
            center(i) = mean(scores(i,1:numberOfCentroidImages));
        end;
    
        %calculate the average distance of each eigenimage to the centroid
        avgDist = zeros(numberOfCentroidImages,1);
        for i=1:numberOfCentroidImages
            d = sqrt(sum((scores(1:numberOfEigenImages,i)-center).^2));
            avgDist(i) = d;
        end;
    
        avgDistMean = mean(avgDist);
        avgDistStd = std(avgDist);

else
    %use calculated training coefficients for classification
    coefs= trainingReference.coefs;
    eigenValues = trainingReference.eigenValues;
    scores=coefs'*reshapedStack';
%     clusterPosition = trainingReference.clusterPosition;
%     avgDistances = trainingReference.avgDistances;
    
    center = trainingReference.center;
    avgDistMean = trainingReference.avgDistMean;
    avgDistStd = trainingReference.avgDistStd;    
    
end;


%% classify images according to their distance from the centroid

goodStack = [];
goodCounter = 1;

for i=1:size(imageStack,3)


    
%     cluster = 0;
%     %     disp([ mat2str(scores(:,i)')]);
%     for j = 1:numberOfClusters
%         %         disp([mat2str(clusterPosition(j,:))]);
%         d = sqrt(sum((scores(:,i)'-clusterPosition(j,:)).^2));
%         
%         
%         if(d <= avgDistances(j,1) + sigma*avgDistances(j,2))
%             cluster = j;
%             %tom_imagesc(imageStack(:,:,i));drawnow;
%         end;
%     end;
%     
%     if(cluster > 0)
%         goodStack(:,:,goodCounter)= tom_norm(originalStack(:,:,i),'mean0+1std');% /d;
%         goodCounter = goodCounter+1;
%         goodFlags(i) = 1;
%     else
%         goodFlags(i) = 0;
%     end;
    if(length(numberOfEigenImages) >1)
        d = sqrt(sum((scores(numberOfEigenImages(1):numberOfEigenImages(2),i)-center(numberOfEigenImages(1):numberOfEigenImages(2))).^2));
    else
        d = sqrt(sum((scores(1:numberOfEigenImages,i)-center(1:numberOfEigenImages)).^2));
    end;
    if(d <= avgDistMean+avgDistStd && goodFlags(i) == 1)
        %plot(i,d,'gx');

        goodStack(:,:,goodCounter)= originalStack(:,:,i);% /d;
        goodCounter = goodCounter+1;
        goodFlags(i) = 1;
    else
        goodFlags(i) = 0;
    end;
end;

% training.clusterPosition = clusterPosition;
% training.avgDistances = avgDistances;
%training.numberOfClusters = numberOfClusters;

training.center = center;
training.avgDistMean = avgDistMean;
training.avgDistStd = avgDistStd;
training.coefs = coefs;
training.eigenValues = eigenValues;
training.originalSize = size(imageStack);
training.numberOfEigenImages = numberOfEigenImages;
training.alignReference = reference;
training.alignTranslationMask = transMask;
training.alignRotationMask = sphere;
