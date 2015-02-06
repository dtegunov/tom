function numberPicks = tom_os3_costFunction(results,currentPosition,align2d)
%tom_os3_costFunction
%
% 
%   tom_os3_costFunction(results,currentPosition,align2d,path)
%
%   Cost function used for simulated annealing. Returns the number of picks
%   until each particle in align2d is found
%PARAMETERS
%
%  INPUT
%   results
%   currentPosition
%   align2d
%   path
%  OUTPUT
%   results - 
%
%EXAMPLE
%   
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_os3_createJobList, tom_os3_corr , tom_os3_
%
%   created by TH2 07/11/07
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

peaks = currentPosition(1) * results.ccc + currentPosition(2) * results.psr + currentPosition(3) * results.autoc;
pickList = tom_os3_returnPicklist(peaks,results.angles,results.job,results.job.options.filter.numberParticles);
picksAl = tom_os3_pickList2Align2d(pickList);


numberPicks = length(pickList);

found = zeros(length(align2d));

radius = align2d(1).radius;

for i=1:length(pickList)
    
    p = picksAl(i);
    
    x1 = p.position.x/2^(results.job.options.modifications.binning+1);
    y1 = p.position.y/2^(results.job.options.modifications.binning+1);

    for j=1:length(align2d)
        if(~found(j) == 1)
            a = align2d(j);

            x2 = a.position.x;
            y2 = a.position.y;

            d = sqrt((x1-x2)^2+(y1-y2)^2);

            if(d <= radius && numberPicks == length(pickList))
                numberPicks = i;
            end;
        end;
    end;
    
    if(sum(found(:) == length(numberPicks)))
        numberPicks = i;
    end;    
end;







