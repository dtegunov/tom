function stack=tom_xmippsellread(sel_file_name,part_idx,flag,binning)
% tom_xmippsellread reads xmipp .sel file
%  
%  
%  tom_xmippsellread(sel_file_name)
%  
%  PARAMETERS
%  
%    INPUT
%     filename     filename of selfile
%     part_idx     nr of particle in selfile
%     flag         'part_nr'
%     binning      read stack binned
%
%    OUTPUT
%     stack        stack
%  
%  EXAMPLE
%      tom_xmippsellread('./data/26S.sel');
%  
%  REFERENCES
%  
%  SEE ALSO
%     tom_mrc2emstack, tom_emread, tom_emwrite
%  
%     created by fb 01/08/07
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

if (nargin<2)
    part_idx='';
end;

if (nargin<3)
    flag='entry';
end;

if (nargin<4)
    binning=0;
end;

i=0;
name='o';
try
    fp=fopen(sel_file_name,'rt');
catch
    error('Cannot open sel-file');
end;

if (isempty(part_idx)==0 && strcmp(flag,'part_nr') )
    part_nr=zeros(1000000,1);
end;
zz=1;
while ~isequal(name,'')
    try
        [name c]= fscanf(fp,'%s\n',1);
        if (isempty(part_idx)==0 && strcmp(flag,'part_nr') )
            [a b c]=fileparts(name);
            [tok rest]=strtok(b,'_');
            rest=strrep(rest,'_','');
            part_nr(zz)=str2double(rest);
            zz=zz+1;
        end;
        
        [nr c]= fscanf(fp,'%s\n',1);
        i=i+1;
    catch
        fclose(fp);
        break;
    end;
end;

if (isempty(part_idx)==0 && strcmp(flag,'part_nr') )
    tmpp=part_nr;
    part_nr=tmpp(1:zz-1);
end;

fp=fopen(sel_file_name,'rt');
[name c]= fscanf(fp,'%s\n',1);

try
    in=tom_spiderread(['./' name]);
    path_pre='./';
catch Me
    in=tom_spiderread(name);
    path_pre='';
end;


if (max(size(binning)) > 1 )
    in.Value=imresize(in.Value,binning);
else
    in.Value=tom_bin(in.Value,binning);
end;

fclose(fp);
szz=size(in.Value);

%tom_emwritec(output_em_file,[size(in.Value,1) size(in.Value,2) i-1],'new','single');

if (isempty(part_idx)==0 && strcmp(flag,'part_nr') )
    stack=zeros(szz(1),szz(2),length(part_idx));
end;

if (isempty(part_idx))
    stack=zeros(szz(1),szz(2),i-1);
end;

fp=fopen(sel_file_name,'rt');
for idx=1:i-1
    if idx./1000==round(idx./1000) disp([num2str(idx) ' particles done']);end;
    [name c]= fscanf(fp,'%s\n',1);
    [nr c]= fscanf(fp,'%s\n',1);
    tmp=tom_spiderread([path_pre name]);
    
    if (max(size(binning)) > 1 )
        stack(:,:,idx)=imresize(tmp.Value,binning);
    else
        stack(:,:,idx)=tom_bin(tmp.Value,binning);
    end;
    
    
    %disp(['index ' num2str(idx) '  name: ' name  ]);
    %tom_emwritec(output_em_file,in.Value,'subregion',[1 1 idx],[size(in.Value,1) size(in.Value,2) 1]);
    
end;
stack=tom_emheader(stack);

fclose(fp);