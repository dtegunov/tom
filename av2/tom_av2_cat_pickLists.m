function test=tom_av2_cat_pickLists(wildcard,outputname)
%TOM_AV2_CAT_PICKLISTS cats tom picklists
%
%
%   align2d=tom_av2_cat_pickLists(wildcard,outputname)
%
%PARAMETERS
%
%  INPUT
%   wildcard           unix wildcard  
%   outputname         name of new picklist                   
%  
%  OUTPUT
%    test           picklist in mem  
%
%EXAMPLE
%  
%  tom_av2_cat_pickLists('align-*.mat','../../auto-pl.mat');
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_bin_sel,tom_xmippsellread
%
%   created by fb 01/08/07
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


dd=dir(wildcard);
bp=fileparts(wildcard);
%disp(['reading ' num2str(length(dd)) ' files']);
if (isempty(bp))
    load([dd(1).name]);
else
    load([bp '/' dd(1).name]);
end;

test=align2d;
for i=2:length(dd)
    if (isempty(bp))
        load([dd(i).name]);
    else
        load([bp '/' dd(i).name]);
    end;
    if (size(align2d,2)>0)
        test=cat(2,test,align2d);
    else
        disp(['warnig ' dd(i).name ' is empty']);
    end;
    
    if (mod(i,100)==0)
 %       disp([num2str(i) ' lists loaded']);
    end;
end;


if (exist('outputname','var') )
    align2d=test;
    save(outputname,'align2d','-v7.3');
end;



