function dist=tom_calc_euler_distance(euler1,euler2)
% TOM_CALC_EULER_DISTANCE 
%     tom_calc_euler_distance(euler1,euler2)
%  
%  PARAMETERS
%  
%    INPUT
%     euler1               euler triple 1
%     euler2               euler triple 2
%
%    
%    OUTPUT
%     dist                distance in degree
%  
%  EXAMPLE
%
%   
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB/ mid july
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom


%transfer 2 quat
 
if (euler1(2) >= 180)
    euler1(2)=179.9;
end;

if (euler1(2) == 0)
    euler1(2)=1;
end;

if (euler2(2) >= 180)
    euler2(2)=179.9;
end;

if (euler2(2) == 0)
    euler2(2)=1;
end;

try
    q1=SpinCalc('EA323toQ',euler1,0.01,0);
catch
    disp(['Error: ' num2str(euler1) ]);
    dist=360;
    return;
end;

try
    q2=SpinCalc('EA323toQ',euler2,0.01,0);
catch
    disp(['Error: ' num2str(euler2) ]);
    dist=360;
    return;
end;


sc=q1(1).*q2(1);
sc=sc+q1(2).*q2(2);
sc=sc+q1(3).*q2(3);
sc=sc+q1(4).*q2(4);

if sc >= 1
    sc = 1.0;
end;

if sc <= -1
    sc = -1.0;
end;
  
dist=acos(abs(sc)).*(180./pi).*2;
