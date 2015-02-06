function tom_av2_xmipp_create_sel_file_from_dir(parts_path,selfile,wildcard,subset)

%tom_av2_xmipp_create_sel_file_from_dir creates a xmipp sel-file from a
%directory file list/
%   
%
%
%PARAMETERS
%
%  INPUT
%   parts_path        directory path with particles
%   selfile           name of selfile
%   wildcard          take wildcard files only, default('*.spi')
%   subset            number of files, random
%
%  OUTPUT
%  
%EXAMPLE
%     
%  tom_av2_xmipp_create_sel_file_from_dir('./parts','./parts.sel','*.spi')
%
%  % random subset of 10000 files
%  tom_av2_xmipp_create_sel_file_from_dir('./parts','./parts.sel','*.spi',10000)
%
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by SN 04/07/10
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
%

if nargin<3
    wildcard='*.spi';
end;

if nargin<2
    error('Syntax: tom_av2_xmipp_create_sel_file_from_dir(parts_path,selfile,[wildcard])');
end;

if parts_path(1)~='/'
    warning('this is a relative path, not an absolute path!');
end;

d=dir([parts_path '/' wildcard]);

try
    fp=fopen(selfile,'wt');
catch
    error(['tom_av2_xmipp_create_sel_file_from_dir: cannot write selfile: ' selfile]);
end;

if nargin==4
    r=randperm(length(d));
    r=r(1:subset);
    for i=1:length(r)
        fprintf(fp,[parts_path '/' d(r(i)).name ' 1 ' '\n']);
    end;
    disp([selfile ' : ' num2str(length(r)) ' files in total.']);
else
    for i=1:length(d)
        fprintf(fp,[parts_path '/' d(i).name ' 1 ' '\n']);
    end;
    disp([selfile ' : ' num2str(length(d)) ' files in total.']);
end;
fclose(fp);

