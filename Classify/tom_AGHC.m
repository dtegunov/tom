function [patterns, targets] = tom_AGHC(train_patterns, train_targets, params, plot_on)
%TOM_AGHC creates ...
%
%   [patterns, targets] = tom_AGHC(train_patterns, train_targets, params, plot_on)
%
%PARAMETERS
%
%  INPUT
%   train_patterns      ...
%   train_targets       ...
%   params              ...
%   plot_on             ...
%  
%  OUTPUT
%   patterns    		...
%   patterns    		...
%
%EXAMPLE
%   ... = tom_AGHC(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

error(nargchk(0, 4, nargin, 'struct'))

if (nargin < 4),
    plot_on = 0;
end

[c, method] = process_params(params);
[D,c_hat]	= size(train_patterns);
label       = 1:c_hat;
n           = ones(1,c_hat);

%Compute distances
N           = size(train_patterns,2);
temp        = repmat(train_patterns,[1 1 N]);
dist        = sqrt(squeeze(sum((temp - permute(temp, [1 3 2])).^2)));

while (c_hat > c),
    Uc       = unique(label);
    Nc       = length(Uc);
    new_dist = zeros(Nc);

    switch method
    case 'min'
        %Find minimum distance between vectors from different clusters
        
        %For each two clusters, find the shortest distance between vectors
        for i = 1:Nc,
            i_in = find(label == Uc(i));
            for j = 1:Nc,
                j_in = find(label == Uc(j));
                new_dist(i,j) = min(min(dist(i_in,j_in)));
            end
        end
        new_dist    = new_dist + eye(Nc)*1e33;
        [i,j]   = find(new_dist == min(min(new_dist)));
        i = Uc(i(1)); j = Uc(j(1));
    case 'max'
        %Find maximum distance between vectors from different clusters
        
        %For each two clusters, find the longest distance between vectors
        for i = 1:Nc,
            i_in = find(label == Uc(i));
            for j = 1:Nc,
                j_in = find(label == Uc(j));
                new_dist(i,j) = max(max(dist(i_in,j_in)));
            end
        end
        new_dist = new_dist .* (ones(Nc)-eye(Nc));
        [i,j]   = find(new_dist == max(max(new_dist)));
        i = Uc(i(1)); j = Uc(j(1));
        
    case 'avg'
        %Find average distance between vectors from different clusters
        
        %For each two clusters, find the average distance between vectors in one cluster to each vector in the other cluster
        for i = 1:Nc,
            i_in = find(label == Uc(i));
            for j = 1:Nc,
                j_in = find(label == Uc(j));
                new_dist(i,j) = mean(mean(dist(i_in,j_in)))/(length(j_in)*length(i_in));
            end
        end
        new_dist = new_dist .* (ones(Nc)-eye(Nc));
        [i,j]   = find(new_dist == max(max(new_dist)));
        i = Uc(i(1)); j = Uc(j(1));
        
    case 'mean'
        %Find mean distance between cluster centers 
        
        %For each two clusters, find the average distance between vectors in one cluster to each vector in the other cluster
        for i = 1:Nc,
            i_in = find(label == Uc(i));
            for j = 1:Nc,
                j_in = find(label == Uc(j));
                new_dist(i,j) = sum((mean(train_patterns(:,i_in)')'-mean(train_patterns(:,j_in)')').^2);
            end
        end
        new_dist    = new_dist + eye(Nc)*1e33;
        [i,j]   = find(new_dist == min(min(new_dist)));
        i = Uc(i(1)); j = Uc(j(1));
    otherwise
        error('Distance method unknown')
    end
      
    %Merge cluster i with cluster j
    label(find(label == j)) = i;
    
    c_hat = c_hat - 1;
    
    %Computer cluster centers
    Uc       = unique(label);
    Nc       = length(Uc);
    patterns = zeros(D,Nc);
    for i = 1:Nc,
        in            = find(label == Uc(i));
        if (length(in) == 1)
            patterns(:,i) = train_patterns(:,in);
        else
            patterns(:,i) = mean(train_patterns(:,in)')';
        end
    end
    
    %Plot the centers during the process 
    plot_process(patterns, plot_on)
    
end
 
%Label the data
targets = zeros(1,c);
Uc      = unique(label);
Ut      = unique(train_targets);
targets = zeros(1,c);
for i = 1:c,
    indices    = find(label == Uc(i));
    N          = hist(train_targets(:,indices), Ut);
    [m, max_l] = max(N);
    targets(i) = Ut(max_l);
    if (length(indices) == 1)
        patterns(:,i) = train_patterns(:,indices);
    else
        patterns(:,i) = mean(train_patterns(:,indices)')';
    end
end

