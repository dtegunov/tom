function tom_av2_picker(train_st,images,outputdir)
%tom_av2_picker_train picks particles with a trained struct
%
%   tom_av2_picker(train_st,image_list,outputdir)
%
%PARAMETERS
%
%  INPUT
%   train_st            structure containing the training from tom_av2_picker_train
%   images              can be:  cell of images, a wildcard or a filefilter  
%   outputdir           picklist are written 2 this folder
%
%
%  OUTPUT
%
%
%EXAMPLE
%    
% matlabpool open local 8;
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_picker_train, tom_av2_compare_picklists 
%
%   created by FB 11/11/10
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


if (isstruct(train_st)==0)
    load(train_st);
end;

if (iscell(images))
    image_list=images;
end;

if (ischar(images))
    [a b c]=fileparts(images);
    dd=dir(images);
    for i=1:length(dd)
        image_list{i}=[a '/' dd(i).name];
    end;
end;

parfor i=1:length(image_list)
%for i=1:length(image_list)
  match(image_list{i},train_st,outputdir)
end;



function match(img_name,options,outputdir)

%read image
disp(['Reading ' img_name]);
im=tom_emreadc(img_name);
im=im.Value;
if (options.micrograph.invert==1)
    im=-im;
end;
%reacale image 2 library
scaling_f=(options.micrograph.pixelsize./options.template.pixelsize);
new_size=round(size(im).*scaling_f);
im=single(imresize(im,new_size));

%do the matching
disp(['Matching ' img_name]);
try
    peaks=tom_os3_match_templ(im,options.template.dir,options.peaks);
catch Me
    disp(['Error Matching Image: ' img_name ' ...skipping!!']);
    disp(Me.message);
    return;
end;

clear('im');
[coord,align2d]=tom_os3_extract_peaks2(peaks,options.peaks.num_of_peaks,options.template.size,img_name,(1./scaling_f),options.template.zero_rad);

%filter by cc
org_length_align=size(align2d,2);
align2d=align2d(1,find([align2d(1,:).ccc] > options.peaks.cc_thresh));
disp(['Filtering by matching ' img_name ' org ' num2str(org_length_align) ' filtered ' num2str(size(align2d,2)) ]);


%refine matchign by classification
if (options.mult_cl.do_multiref==1)
    %generate particle stack from picklist
    scaling_f=(options.micrograph.pixelsize./options.mult_cl.pixelsize);
    sz_mult_lib=[size(options.mult_cl.cents,1) size(options.mult_cl.cents,2)];
    rad=round((sz_mult_lib(1) .* (1./scaling_f))./2);
    norm_str=['local_' num2str(round(rad./1.9)) '&mean0+1std'];
    stack=tom_av2_xmipp_picklist2stack(align2d,'','',rad,norm_str,options.micrograph.invert,sz_mult_lib);
    
    %filter the stack
    if (strcmp(options.mult_cl.filter_type,'none')==0)
        if (strcmp(options.mult_cl.filter_type,'kernel') )
            for i=1:size(stack,3);
                stack(:,:,i)=tom_filter(stack(:,:,i),options.mult_cl.filter_values(1));
            end;
        end;
    end;
    %do the classification
    disp(['Classifying ' img_name  ]);
    [classes ccc_out]=tom_av2_cl_multiref_classify(options.mult_cl,stack);
    %transfer values
    for i=1:length(align2d)
        align2d(1,i).ccc_match=align2d(1,i).ccc;
        align2d(1,i).ccc=ccc_out(i);
    end;
    %filter the picklist
    org_length_align=size(align2d,2);
    align2d=align2d(1,find(classes==1));
    disp(['Filtering by multiref classification ' img_name ' org ' num2str(org_length_align) ' filtered ' num2str(size(align2d,2)) ]);
else
   %transfer values
    for i=1:length(align2d)
        align2d(1,i).ccc_match=align2d(1,i).ccc;
    end; 
end;

%write output
[a b c]=fileparts(img_name);
save([outputdir '/pick_' b c '.mat'],'align2d');






