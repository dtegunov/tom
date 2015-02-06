function tom_av2_prepare_refitctf(flag,filenames,astiborder)

%TOM_AV2_REFITCTF checks the refit-folder for already existing files and
%moves them to the respective original folder; then it checks the .mat-files from a
%previous ctf-fit and moves images to the refit-folder for refitting, if
%one of the following is true: fitted defocus is equal to the lower or
%upper boundary of the ctf-fit-search criteria; the fitted astigmatism is
%greater then 0.15e-6 or the given value (asti_border)
%
%   tom_av2_prepare_refitctf(flag,filenames,astiborder)
%
%PARAMETERS
%
%  INPUT
%   flag                    'low' or 'high'
%   filenames               wildcard for number in filename
%   astiborder              value of maximal astigmatism, standard is 1.8e-7
%
%EXAMPLE
%
%   tom_av2_prepare_refitctf('low','low/sc26s_rpn1gfp_*.em.mat',2.2e-7);
%
%REFERENCES
%
%SEE ALSO
%
%   created by SB 10/02/24
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

if nargin<3;
    astiborder=1.8e-7;
end;

dd=dir(filenames);

if strcmp(flag,'low')==1
    
    c=dir('lowrefit');
    if length(c)==0;
        unix('mkdir lowrefit');
    end;
    if size(c,1)>2;
        unix('mv lowrefit/* low');
        dd=dir(filenames);
    end;
    for k=1:length(dd);
        try
            [a b c]=fileparts(dd(k).name);
            [bb bc]=fileparts(b);
            load(['low/' dd(k).name]);
            if st_out.Fit.Dz_delta_det>astiborder && st_out.sel.accepted==1;
                disp(['mv low/' bc '.* lowrefit']);
                unix(['mv low/' bc '.* lowrefit']);
            end;
            if st_out.Search.Dz_search(1)==st_out.Fit.Dz_det && st_out.sel.accepted==1;
                disp(['mv low/' bc '.* lowrefit']);
                unix(['mv low/' bc '.* lowrefit']);
            end;
            if st_out.Search.Dz_search(size(st_out.Search.Dz_search,2))==st_out.Fit.Dz_det && st_out.sel.accepted==1;
                disp(['mv low/' bc '.* lowrefit']);
                unix(['mv low/' bc '.* lowrefit']);
            end;
        catch
            disp(['Image' num2str(k) ' not existent.']);
        end;
    end;
else

    c=dir('highrefit');
    if length(c)==0;
        unix('mkdir highrefit');
    end;
    if size(c,1)>2;
        unix('mv highrefit/* low');
        dd=dir(filenames);
    end;
    for k=1:length(dd);
        try
            [a b c]=fileparts(dd(k).name);
            [bb bc]=fileparts(b);
            load(['high/' dd(k).name]);
            if st_out.Fit.Dz_delta_det>astiborder && st_out.sel.accepted==1;
                disp(['mv high/' bc '.* highrefit']);
                unix(['mv high/' bc '.* highrefit']);
            end;
            if st_out.Search.Dz_search(1)==st_out.Fit.Dz_det && st_out.sel.accepted==1;
                disp(['mv high/' bc '.* highrefit']);
                unix(['mv high/' bc '.* highrefit']);
            end;
            if st_out.Search.Dz_search(size(st_out.Search.Dz_search,2))==st_out.Fit.Dz_det && st_out.sel.accepted==1;
                disp(['mv high/' bc '.* highrefit']);
                unix(['mv high/' bc '.* highrefit']);
            end;
        catch
            disp(['Image' num2str(k) ' not existent.']);
        end;
    end;
end;