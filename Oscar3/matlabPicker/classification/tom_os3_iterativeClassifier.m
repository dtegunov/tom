function [ goodFlags goodStack] = tom_os3_iterativeClassifier(stack,trainingStack,options,plotON)
%(stack,numberTrainingImages,numberOfClusters,sigmaScale)
%tom_os3_iterativeClassifier
%   
%   tom_os3_iterativeClassifier
%
%PARAMETERS
%
%  INPUT
%   stack       - the image stack          
%   training    - the training stack of preselected particles  
%   NOTE        - both stacks must be masked and normalised BEFORE
%   classification. reshaping is applied here if necessary
%   options     - the options structure of the particle picker
%
%  OUTPUT
%   
%   goodFlags   - array of ones and zeros marking good particles by one
%   goodStack   - the stack of "good" particles
%   
%
%%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by TH2 3/12/07
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


%% if set, online plotting of classification
if(~exist('plotON'))
    plotON = false;
end;

if(plotON)
    figHandle = figure;
end;

%%

numberOfClusters = options.classification.numberOfClusters;
eigenStart = options.classification.eigenStart;
eigenEnd = options.classification.eigenEnd;
sigmaScale = options.classification.sigmaScale;

%% if k means clustering could not be processed on the dataset, return the
%  original stack

if(numberOfClusters<1)
    goodStack = stack;
    goodFlags = ones(ones(size(stack,3),1));
    return;
end;

% eigenStart = 1;
% eigenEnd = 20;
 eigenCounter = 0;
%oldStackShow = stack;

%% reshape image stack
img =[];
if(numel(size(stack))>2)
    for i = 1:size(stack,3)
        img(i,:) = reshape(tom_bin(stack(:,:,i)),1,[]);
    end;
else
    img = stack;
end;
%% wurde veraendert. vorher war der trainingstack teil des stacks, ersten
% 100 partikel, jetzt separat
% training = img(1:100,:);
% trainingStack = stack(:,:,1:100);
if(numel(size(trainingStack))>2)
    for i = 1:size(trainingStack,3)
        training(i,:) = reshape(tom_bin(trainingStack(:,:,i)),1,[]);
    end;
else
    training = stack;
end;


goodStack = [];
goodFlags = ones(size(stack,3),1);


%% determine the number of steps
if(length(sigmaScale) <3)
    step = -(sigmaScale(1)-sigmaScale(2))/4;
else
    %make sure step is negative. iteration loop would fail otherwise
    if(sigmaScale(1)>sigmaScale(2) && sigmaScale(3) > 0)
        step =-sigmaScale(3);
    else
        step =sigmaScale(3);
    end;
end;

%% perform iterative classification here
for distance = sigmaScale(1):step:sigmaScale(2)
    
%   concatenate training and image stack (both are reshaped now) for pca  
    contStack = cat(1,training,img);
    
%   calculate pca for the image stack and training stack
    [scores,coefs, eigenvalues]=tom_calc_pca(double(contStack),eigenEnd,'pca','','','','','','',0);

    if(plotON)
        tmp = reshape(coefs,sqrt(size(coefs,1)),sqrt(size(coefs,1)),[]);
        subplot(5,1,1,'replace');tom_dspcub(tmp);drawnow;
    end;
%   repeat k means clustering for 10 times. if not successfull, try a classification with less clusters    
    kmeansProcessed = false;
    i=1;
    while(~kmeansProcessed && i < 10)
        try
            %k means goes here
            [classifications,centroids,distances,interClusterDistances] =kmeans(scores(eigenStart:eigenEnd,1:size(training,1))',numberOfClusters);
            distances = sqrt(distances);
            interClusterDistances = sqrt(interClusterDistances);
            kmeansProcessed = true;
        catch
            disp('Error during clustering. Will retry. Choose lower number of clusters?');
            i = i+1;
        end;
    end;
    
    if(~kmeansProcessed)    
        disp(['NumberOfClasses ' mat2str(numberOfClusters) 'not possible']);
        %recursion starts here
        options.classification.numberOfClusters = numberOfClusters - 1;
        [goodStack goodFlags] = iterativePCA(stack,trainingStack,options);    
        return;
    end;
    
%%  calculate the centroids and deviations for each cluster


    %for each cluster
    for clusterIterator=1:options.classification.numberOfClusters
        
        %set distance vector
        dist =[];
        
        %find members of i.th cluster
        memberIndex = find(classifications == clusterIterator);
        
        dist = interClusterDistances(memberIndex,clusterIterator);
        
        avgDistances(clusterIterator,1) = mean(dist);
        avgDistances(clusterIterator,2) = std(dist,1);
        
    end;

%     newReshapedStack = training;
%     newStack = trainingStack;
    positionInNewStack = size(training,1) + 1;
    newStack = [];
    
    
    if(plotON)
        
        subplot(5,1,2,'replace');hold on;
        for(i=1:positionInNewStack-1)      
            plot(scores(1,i),scores(2,i),'b.','MarkerSize',1);
        end;
        
        for i=1:numberOfClusters
            plot(centroids(i,1),centroids(i,2),'k.','MarkerSize',15);
        end;

        subplot(5,1,5,'replace');hold on;
        for i=1:size(trainingStack,3)
            plot(i,interClusterDistances(i,classifications(i)),'.g');
        end;
    end;
%%  classify the image stack here
    for i=1:length(goodFlags)
        if(goodFlags(i) == 1)
            cluster = 0;
            %     disp([ mat2str(scores(:,i)')]);
            for j = 1:numberOfClusters
                %         disp([mat2str(clusterPosition(j,:))]);
                
                if(cluster == 0)
                    d = sqrt(sum((centroids(j,:)'-scores(eigenStart:eigenEnd,positionInNewStack)).^2));

                    if((d <= avgDistances(j,1) + distance*avgDistances(j,2)))
                        cluster = j;
                        %tom_imagesc(imageStack(:,:,i));drawnow;
                    end;
                end;
            end;

            if(cluster >0)
                 
%                 newReshapedStack(size(newReshapedStack,1)+1,:) = img(positionInNewStack,:);
%                 newStack(:,:,size(newStack,3)+1) = stack(:,:,positionInNewStack);
                  newStack(size(newStack,1)+1,:) = contStack(positionInNewStack,:);
                  if(plotON)
                    subplot(5,1,5);plot(i+size(trainingStack,3),d,'go');
                    subplot(5,1,2);plot(scores(1,positionInNewStack),scores(2,positionInNewStack),'g.','MarkerSize',5);
                  end;
            else
                 
%                 oldStackShow(:,:,i) = oldStackShow(:,:,i) *0;
                  goodFlags(i) = 0; 
                  if(plotON)
                    subplot(5,1,5);plot(i+size(trainingStack,3),d,'rx');
                    subplot(5,1,2);plot(scores(1,positionInNewStack),scores(2,positionInNewStack),'r.','MarkerSize',5);
                  end;
            end;
            positionInNewStack  = positionInNewStack +1;
        end;
    end;
%    hold off;
    img = newStack;    
    if(plotON)

        tmp = reshape(img',sqrt(size(coefs,1)),sqrt(size(coefs,1)),[]);
        subplot(5,1,3,'replace');tom_dspcub(tmp);drawnow; 
        subplot(5,1,4,'replace');tom_imagesc(sum(tmp,3));drawnow;
        
        ginput(1);
    end;
    clear newStack;
% figure;tom_dspcub(stack(:,:,101:end))    

%     img = newReshapedStack;
%     clear newReshapedStack;
%     figure;tom_dspcub(stack(:,:,101:end));

    eigenCounter = eigenCounter +1;
end;


goodFlags = goodFlags(1:end);
goodStack = [];
