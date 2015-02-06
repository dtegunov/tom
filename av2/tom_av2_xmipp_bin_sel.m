function tom_av2_xmipp_bin_sel(f_sel,repl,repl_to,bin,new_sel,cut_size)
%TOM_AV2_XMIPP_BIN_SEL bins particles according to input sel
%
%   tom_av2_xmipp_bin_sel(f_sel,repl,repl_to,bin,new_sel)
%
%  TOM_AV2_XMIPP_BIN_SEL bins particles according to input sel and
%  writes out a new folder with binned particel and new sel-file
%  
%
%PARAMETERS
%
%  INPUT
%   f_sel       filename of the input sel
%   repl        string in the sel which should be replaced
%   repl_to     new string (according 2 find and replace in a text editor)
%   bin         binning of the images or new size
%   new_sel     name of the new sel file 
%   cut_size    (opt.) cut out the image before and rescale or bin it after                     
%
%EXAMPLE
%  
%  %bin   
%  tom_av2_xmipp_bin_sel('13_corrf_high_128_clean.sel','/parts_high_128','/parts_high_64',1,'13_corrf_high_64_clean.sel'); 
%   
%  %resale   
%  tom_av2_xmipp_bin_sel('13_corrf_high_128_clean.sel','/parts_high_128','/parts_high_64',[128 128],'13_corrf_high_64_clean.sel'); 
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_aling2d
%
%   created by fb ...ole !!
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

if (nargin < 6)
    cut_size=-1;
end;


d=importdata(f_sel);
tmp=tom_spiderread(d.textdata{1});  

if (bin(1) > size(tmp.Value,1) )
    mask=tom_spheremask(ones(bin(1),bin(2)),round(size(tmp.Value,1)./2)-4,3);
end;

fid=fopen(new_sel,'wt');
tic;
for i=1:length(d.textdata)
    if(exist(d.textdata{i},'file'))
        f_name=d.textdata{i};
        if (isempty(strfind(f_name,repl)) )
            disp(['Token ' repl ' not found in ' f_name  '...skipping']);
            continue;
        end;
        
        [foldername filename f_ext]=fileparts(f_name);
        new_fold_name=[foldername(1:max(strfind(foldername,'/'))-1 ) strrep([foldername(max(strfind(foldername,'/')):end) '/'],repl,repl_to)];
        if(exist(new_fold_name,'dir')==0)
            mkdir(new_fold_name);
        end;
        im=tom_spiderread(f_name);    
        
        if (cut_size(1)~=-1)
            im.Value=tom_cut_out(im.Value,'center',cut_size);
        end;
        
        if (length(bin) > 1)
           % if (exist('mask','var'))
                %im=tom_rescale(im.Value,[bin(1) bin(2)],mask);
                im=imresize(im.Value,[bin(1) bin(2)]);
                %else
             %   im=tom_rescale(im.Value,[bin(1) bin(2)],mask); 
            %end;
        else
            im=tom_bin(im.Value,bin);
        end;
        %new_name=strrep(f_name,repl,repl_to);
        new_name=[new_fold_name filename f_ext];
        if (strcmp(new_name,f_name)==1)
            error(['new name: ' new_name ' == old name: ' f_name]);
        end;
        tom_spiderwrite(new_name,im);
        fprintf(fid,[new_name ' 1 \n']);
        
    else
        disp([d.textdata(i) 'not found ...skipping']);
    end;

    if (mod(i,1000)==0)
        toc;
        disp(num2str(i));
        tic;
    end;
    
end;

fclose(fid);
