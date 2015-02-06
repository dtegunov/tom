function [patterns, targets, label] = tom_kohonen_SOFM(train_patterns, train_targets, params, plot_on)
%TOM_KOHONEN_SOFM creates ...
%
%   [patterns, targets, label] = tom_kohonen_SOFM(train_patterns, train_targets, params, plot_on)
%
%PARAMETERS
%
%  INPUT
%	train_patterns		Input patterns
%	train_targets		Input targets
%	params				[Number of output data points, Window width]
%   plot_on         	Plot stages of the algorithm
%  
%  OUTPUT
%	patterns			New patterns
%	targets				New targets
%	label				The labels given for each of the original patterns
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

[Nmu, win_width] = tom_process_params(params);
if (nargin < 4),
    plot_on = 0;
end

[D,L]	= size(train_patterns);
dist	= zeros(Nmu,L);
label   = zeros(1,L);

%Initialize W
W			= sqrtm(cov(train_patterns',1))*randn(D,Nmu);
W			= W ./ (ones(D,1)*sqrt(sum(W.^2)));
dW			= 1;

%Learning rate
eta		= 0.5;
deta		= 0.995;
iter		= 0;

while (dW > 1e-15),
   %Choose a sample randomally
   i		= randperm(L);
   phi	= train_patterns(:,i(1));
   
   net_k = W'*phi;
   y_star= find(net_k == max(net_k));
   y_star= y_star(1); %Just in case two have the same weights!
   
   oldW	= W;
   W		= W + eta*phi*gamma(win_width*abs(net_k - y_star))';
   W		= W ./ (ones(D,1)*sqrt(sum(W.^2)));
   
   eta	= eta * deta;
   
   dW		= sum(sum(abs(oldW-W)));
   iter	= iter + 1;   
   
   if (plot_on > 0),
      %Assign each of the patterns to a center
      dist        = W'*train_patterns;
      [m, label]  = max(dist);
      centers     = zeros(D,Nmu);
      for i = 1:Nmu,
         in = find(label == i);
         if ~isempty(in)
            centers(:,i) = mean(train_patterns(:,find(label==i))')';
         else
            centers(:,i) = nan;
         end
      end

      %Plot centers during training
     % plot_process(centers, plot_on)
   end
   
   if (iter/100 == floor(iter/100)),
      disp(['Iteration number ' num2str(iter)])
   end
   
end

%Assign a weight to each pattern
label = zeros(1,L);
for i = 1:L,
   net_k 	= W'*train_patterns(:,i);
   label(i) = find(net_k == max(net_k));
end

%Find the target for each weight and the new patterns
targets 	= zeros(1,Nmu);
patterns	= zeros(D, Nmu);
Uc          = unique(train_targets);

for i = 1:Nmu,
    in				= find(label == i);
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


function G = gamma(dist)
%The activation function for the SOFM
G = exp(-dist);