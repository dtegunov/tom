editfunction [patterns, targets, label, W] = tom_leader_follower(train_patterns, train_targets, params, plot_on)
%TOM_LEADER_FOLLOWER creates ...
%
%   [patterns, targets, label, W] = tom_leader_follower(filelabel,threshold, label, color, transformmatrix, iconposition, host)
%
%PARAMETERS
%
%Reduce the number of data points using the basic leader-follower clustering algorithm
%
%  INPUT
%	train_patterns		Input patterns
%	train_targets		Input targets
%	params 				Algorithm parameters: [Min distance to connect, Rate of convergence]
%   plot_on             Plot stages of the algorithm
%
%  OUTPUT
%	patterns			New patterns
%	targets				New targets
%	label				The labels given for each of the original patterns
%   W               	Weights matrice
%
%EXAMPLE
%   ... = tom_amira_createisosurface(...);
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
[theta, eta] = tom_process_params(params);

[D,L]	= size(train_patterns);

%Preprocessing
patterns = [train_patterns ; ones(1,L)];
patterns = patterns ./ (ones(D+1,1) * sqrt(sum(patterns.^2)));

%w1 <- x
w       = patterns(:,1);

for i = 2:L,
    %Accept new pattern x
    x   = patterns(:,i);
    
    %j <- argmin||x-wj|| (Find nearest cluster)
    dist  = sqrt(sum((w - x*ones(1,size(w,2))).^2));
    j     = find(min(dist) == dist);
    
    %if ||x-wj|| < theta
    if dist(j) < theta,
        %wj <- wj + eta*x
        w(:,j) = w(:,j) + eta*x;
    else
        %Add new w <- x
        w(:,end+1) = x;
    end
    
    w   = w ./ (ones(D+1,1) * sqrt(sum(w.^2)));
    
    if (plot_on > 0),
        %Assign each of the patterns to a center
        dist        = w'*patterns;
        [m, label]  = max(dist);
        centers     = zeros(D,size(w,2));
        for i = 1:size(w,2),
            in = find(label == i);
            if ~isempty(in)
                centers(:,i) = mean(train_patterns(1:2,find(label==i))')';
            else
                centers(:,i) = nan;
            end
        end
        %Plot centers during training
        plot_process(centers, plot_on)
    else
        disp(['There are ' num2str(size(w,2)) ' clusters so far'])
    end
    
end 
       
%Assign each of the patterns to a center
N           = size(w,2);
dist        = w'*patterns;
[m, label]  = max(dist);
patterns    = zeros(D,N);
targets     = zeros(1,N);
Uc          = unique(train_targets);
for i = 1:N,
    in = find(label == i);
    if ~isempty(in),
        h            = hist(train_targets(in), Uc);
        [m, best]    = max(h);
        targets(i)	 = Uc(best);
        if length(in) == 1,
            patterns(:,i)	= train_patterns(:,in);
        else
            patterns(:,i)  = mean(train_patterns(:,in)')';
        end
    else
        patterns(:,i) = nan;
    end   
end

