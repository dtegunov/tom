function [avg_out number_mat]=tom_av2_xmipp_somvs2stack(vs_file,parts_file,output_fold,class)
%  TOM_AV2_XMIPP_SOMVS2STACK transforms xmipp som vs 2 em stack,sel and doc files
%                            for visualisation and further processing 
%  
%     avg_out=tom_av2_xmipp_somvs2stack(vs_file,parts_file,output_fold,class)
%  
%  PARAMETERS
%  
%    INPUT
%     vs_file        .vs form xmipp-kerdensom
%     parts_file     .sel or .doc file with the used particles
%     output_fold    ('out_som') folder for file output 
%     class          ('all') 2d vect for the position in the grid 
%                    or use 'all' to process all centroids   
%    OUTPUT
%     avg_out        class average 
%     number_mat     number of particles per class as matrix
%
%  EXAMPLE
%  
%  matlabpool open local 8;
%
%  [avg numb]=tom_av2_xmipp_somvs2stack('som_9.vs','ml2d_ref000004.doc','/fs/scratch/out_som','all'); 
%  out=tom_dspcub2(avg,[9 10]); figure; tom_imagesc(out);
%  
%  %to select certain classes use tom_av2_stackbrowser and gen filtered
%  %align2d struct of out_som/avg_st.em and out_som/avg_st.mat 
%  %use tom_av2_xmipp_filter_som 2 generate a sel and doc for further
%  %processing
%
%
%  %for one class use
%  avg=tom_av2_xmipp_somvs2stack('som_9.vs','ml2d_ref000004.doc','/fs/scratch/out_som',[0 0]); 
%  figure; tom_imagesc(avg);
%  
%  
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_av2_xmipp_filter_som,tom_av2_stackbrowser
%  
%     created by FB 03/02/10
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

[a b c]=fileparts(parts_file);

if (strcmp(c,'.doc'))
    doc_flag=1;
else
    doc_flag=0;
end;    

warning off; mkdir(output_fold); warning on;

[aa bb cc]=fileparts(vs_file);
if (isempty(aa))
    aa='.';
end;

if (doc_flag)
    unix(['head -3 ' parts_file ' > tmp_muratXXX.doc']);
    doc=tom_xmippdocread('tmp_muratXXX.doc');
    unix('rm tmp_muratXXX.doc');
end;

im_tmp=tom_spiderread(doc(1).name);
sz=size(im_tmp.Value);

if (strcmp(class,'all'))
     
     call=['cat ' aa '/' bb  '.inf | grep "Horizontal dimension (Xdim) = " | awk ''{print $5}''' ];
     [tmp num_h]=unix(call);
     cl1=0:(str2double(num_h)-1);
     call=['cat ' aa '/' bb  '.inf | grep "Vertical dimension (Ydim) = " | awk ''{print $5}''' ];
     [tmp num_l]=unix(call);
     cl2=0:(str2double(num_l)-1); 
     unix(['head -10 ' aa '/' bb  '.inf']);   
 else
    cl1=class(1);
    cl2=class(2);
end;
 
zz=1;
avg_out=zeros(sz(1),sz(2),length(cl2)*length(cl1));
tom_emwritec([output_fold '/avg_st.em'],avg_out);
align2d=tom_av2_create_alignfromstack([output_fold '/avg_st.em']);
number_mat=zeros(length(cl2),length(cl1));

for iy=1:length(cl2)
    for ix=1:length(cl1)
        if (ix==1 && iy==2)
            disp(' ');
        end;
        disp(['Processing: class_' num2str(cl1(ix)) '_' num2str(cl2(iy))]);
        output_base=[output_fold '/cl_' num2str(cl1(ix)) '_' num2str(cl2(iy)) ];
        call=['cat ' vs_file ' | grep "^' num2str(cl1(ix)) ' ' num2str(cl2(iy)) ' " | awk ''{split($4,a,"/"); print a[2]}''' ' > ' output_base '.sel'];
        [a b]=unix(call);
        if (a==1)
            error(b);
        end;
        [a b]=unix(['head -1 ' parts_file ' > ' output_base '.doc']);
        call=['fgrep -f ' output_base '.sel '  parts_file ' -A1 | grep -v ''\--''' ' >> ' output_base '.doc'];
        unix(call);
        if (a==1)
            error(b);
        end;
        avg=tom_av2_xmipp_doc2em([output_base '.doc']);
        disp([num2str(size(avg,3)) ' particles found ']);
        avg_out(:,:,zz)=tom_norm((sum(avg,3)./size(avg,3)),'mean0+1std');
        align2d(1,zz).filename=[output_fold '/cl_' num2str(cl1(ix)) '_' num2str(cl2(iy)) ];
        align2d(1,zz).dataset=zz;
        zz=zz+1;
        number_mat(iy,ix)=size(avg,3);
        clear('avg');
        if (doc_flag==1)
            tom_xmippdoc2sel([output_base '.doc'],[output_base '.sel']);
        end;
       
    end;
end;
disp('Writing stacks');
tom_emwritec([output_fold '/avg_st.em'],avg_out);
save([output_fold '/avg_st.mat'],'align2d');

unix(['cat ' aa '/' bb  '.inf']);   


