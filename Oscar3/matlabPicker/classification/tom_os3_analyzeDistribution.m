function newPickList = tom_os3_analyzeDistribution(pickList,sigmaWidth)
%tom_os3_analyzeDistribution
%   
% 
% 
%   tom_os3_analyzeDistribution(pickList,sigmaWidth)
%
%PARAMETERS
%
%  INPUT
%       pickList    - the old picklist
%       sigmaWidth  - sigmaWidth 
%
%  OUTPUT
%       newPickList
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


curve = ones(length(pickList),1);
for i=1:length(curve)
    pick = pickList{i};
    curve(i) = pick.value;
end;

d1=diff(curve);

meanGradient = mean(d1);
stdGradient = std(d1);

newPickList = {};
lastGoodIndex = numel(pickList);

for i =1:length(d1)
    
    if(abs(d1(i)) > abs(meanGradient) + sigmaWidth*stdGradient)
        lastGoodIndex = i;
    end;
    
end;

for i=1:lastGoodIndex
    newPickList{length(newPickList)+1} = pickList{i};
end;
