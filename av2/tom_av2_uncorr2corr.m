function tom_av2_uncorr2corr(f_align,new_f_align)
%tom_av2_uncorr2corr changes an filenames in alignment-list form uncorr 2 corr
%
%   tom_av2_uncorr2corr(f_align,new_f_align)
%    append a _corr 2 every filename in the list
%    corrected images have 2 be in a folder XXX_corr
%
%
%PARAMETERS
%
%  INPUT
%   f_align                 filename for input align2d 
%   new_f_align             output filename
%
%EXAMPLE
%
% tom_av2_uncorr2corr('13_uncorr_high_128_clean_reorg.mat','test.mat');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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




tmp=load(f_align);
align2d=tmp.align2d;
clear('tmp');

h = waitbar(0,'Please wait...');

f_name_old=align2d(1,1).filename;

[a b c]=fileparts(f_name_old);
tmp=['/' a(max(strfind(a,'/'))+1:end) '/'];


zz=1;
img_count=1;
for i=1:size(align2d,2)
    if (strcmp(f_name_old,align2d(1,i).filename)==0 || i==1  )
        clear('st_out');
        try
            load([align2d(1,i).filename '.mat']);
            f_name_old=align2d(1,i).filename;
            if (st_out.sel.selected==0)
                disp([align2d(1,i).filename ' not accepted!']);
                img_count=img_count+1;
            end;
        catch ME
            st_out.sel.selected=1;
            disp([align2d(1,i).filename '.mat not found!']);
            disp([ME.identifier ME.message]);
        end;
    end;
    if (st_out.sel.selected==1)
        align2d(1,zz)=align2d(1,i);
        align2d(1,zz).filename=strrep(align2d(1,i).filename,tmp,[tmp(1:end-1) '_corr/']);
        zz=zz+1;
    end;
    
    waitbar(i./size(align2d,2),h);
end;

disp(' ');
disp(' ');
disp([num2str(size(align2d,2)-zz+1) ' particles sorted out!']);
disp([num2str(img_count-1) ' imges sorted out!']);


align2d=align2d(1,1:(zz-1));

save(new_f_align,'align2d');


disp('done..!');

try
close(h);
catch
end;