function  tom_av2_createstack2(inalignname,outalign,outfilestruct,binning,radius,rad_fact,norm_str,align_flag,taper_flag,mask,num_of_packages,file_buffer_size,gappy_write,verbose)
%TOM_AV2_CREATESTACK2 creates a particle stack from an alignment file
%
%    tom_av2_createstack2(inalignname,outalign,outfilestruct,binning,radius,rad_fact,norm_str,align_flag,taper_flag,mask,num_of_packages,file_buffer_size,gappy_write)
%
%PARAMETERS
%
%  INPUT
%   inalignname         full path to input alignment file (.mat)
%   outalign            output align struct
%   outfilestruct       structure for fileoutput
%                       outfilestruct.method             file structure '1stack' 'Nstacks' 'singleFiles'  
%                       outfilestruct.num_of_entries     number of entries
%                       outfilestruct.path               path for output             
%                       outfilestruct.filename           filename 
%                       outfilestruct.ext                extension     
%                       outfilestruct.fileformat         used fileformat (em,spider,imagic)
%                                               
%
%
%   binning             binning
%   radius              radius the particles should be cut out
%   rad_fact            factor for cutting the particles ...for writing
%                       aligned particles without cropping
%   norm_str            string for norming the particle ('mean0+1std','phase'...) 
%   align_flag          1: save particles aligned, 0: save unaligned (optional)
%   taper_flag          flag for tapering particles
%   mask                mask in memory use 'no_mask' not ones() for no mask ...speed!
%   num_of_packages     nuber of packages for paraell computing 
%   file_buffer_size    size of the filebuffer       
%   gappy_write         0: no jumps in the numbering sceme 
%                       1: position in the inalignment file is equal position                             
%                       2: position in the inalignment file is equal
%                          position gap are filled with zeros 
%                          
%  OUTPUT
%   -
%
%EXAMPLE
%   
%
% tom_av2_createstack2('pick_small.mat','',outfilestruct,2,160,1,'mean0+1std',0,0,'no_mask',4);
%
%
%
%
%
%
%
%
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%   cre0ted by FB 06/16/08
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




if nargin < 14 || isempty(verbose)
    verbose=0;
end;

h=waitbar_wrapper(1,1,'reading alignment file ...',verbose,'');

if ~isstruct(inalignname)
    try
        s = load(inalignname);
    catch
        error('Could not load alignment file.');
    end
    alignstruct = s.align2d;
else
    alignstruct = inalignname;
end;


if nargin < 2 || isempty(outalign)
    outalign_filename='part.mat';
end;


if nargin < 3 || isempty(outfilestruct)
    outfilestruct.method='singleFiles'; 
    outfilestruct.num_of_entries=10000;
    outfilestruct.path='part_fold';
    outfilestruct.filename='part_';
    outfilestruct.ext='.em';
    outfilestruct.fileformat='em';
end;

if nargin < 4 || isempty(binning)
    binning=0;
end;

if nargin < 5 || isempty(radius)
    radius = alignstruct(1,1).radius;
end;

if nargin < 6 || isempty(rad_fact)
    rad_fact=1;
end;

if nargin < 7 || isempty(binning)
    norm_str='mean0+1std';
end;

if nargin < 8 || isempty(align_flag)
    align_flag=0;
end;

if nargin < 9 || isempty(taper_flag)
    taper_flag=1;
end;

if nargin < 10 || isempty(mask)
    mask='no_mask';
end;

if nargin < 11 || isempty(num_of_packages)
    num_of_packages=1;
end;

if nargin < 12 || isempty(file_buffer_size)
    file_buffer_size=750;
end;

if nargin < 13 || isempty(num_of_packages)
    gappy_write=1;
end;

%used to cope with parfor parsing outputvariable of the parfor
%is written to hd and read angai ugaahhhh!!
%used to update alignment file
tmp_h_name=[outfilestruct.path '_tmppp.em'];


radius_bin=round(radius./2^binning);
packages=tom_calc_packages(num_of_packages,size(alignstruct,2));


tom_emwritec([tmp_h_name],[3 size(alignstruct,2) 1],'new');


for ii=1:num_of_packages
%parfor ii=1:num_of_packages

    partcounter = 1;
    stackcounter = 1;
    al_shift=[0 0];
    al_rot=0;
    used_part=packages(ii,1);
    tmpstack = zeros(radius_bin*2,radius_bin*2,file_buffer_size,'single');
    h_array=zeros(3,packages(ii,3));
    
    for i=packages(ii,1):packages(ii,2)
     
    
       
    
        %h=waitbar_wrapper(packages(ii,2),i,'cutting particles ',verbose,h);
        
        h_array=update_h_array(h_array,i-packages(ii,1)+1,0,align_flag,norm_str);
        
        if (alignstruct(1,i).use==0)
            continue;
        end;
        
        
        
       %check if image is readable
        try
            header = tom_reademheader(alignstruct(1,i).filename);
        catch
            warning(['error reading file ' alignstruct(1,i).filename]);
            continue;
        end
        imagesize = header.Header.Size;

        x = alignstruct(1,i).position.x;
        y = alignstruct(1,i).position.y;

       
        
        %check coordinates
        [inside_flag,start_cut,size_cut]=check_cutout(x,y,radius,rad_fact,imagesize);
        if (inside_flag==0 && taper_flag==0)
            continue;
        end;

        %cut out particle
        try
            particle = tom_emreadc(alignstruct(1,i).filename,'subregion',start_cut,size_cut,'binning',binning);
        catch
           warning(['error cutting particle: ' alignstruct(1,i).filename ' position: ' num2str(x) ' ' num2str(y)]);
           continue;
        end;
        particle = single(particle.Value);


        %process particle
        if (align_flag==1)
            al_shift=[alignstruct(1,i).shift.x alignstruct(1,i).shift.y];
            al_rot=alignstruct(1,i).angle;
        end;

        try
            particle=process_particle(particle,al_shift,al_rot,radius,rad_fact,binning,mask,norm_str);
        catch
            warning(['error processing particle: ' alignstruct(1,i).filename ' position: ' num2str(x) ' ' num2str(y)]);
            continue;
        end;

        tmpstack(:,:,partcounter)=particle;
        
        if (gappy_write==1)
            used_part_idx(partcounter)=i;
        else
            used_part_idx(partcounter)=packages(ii,1)+used_part-1;
        end;
        
        partcounter = partcounter + 1;
        used_part=used_part+1;
        
        h_array=update_h_array(h_array,i-packages(ii,1)+1,1);
        
        %partcounter
        if (partcounter > file_buffer_size || (i==packages(ii,2) && partcounter > 2))
            if (i==packages(ii,2) && partcounter > 2)
                tmpstack = tmpstack(:,:,1:partcounter-1);
            end;
            %h=waitbar_wrapper(packages(ii,2),i,'writing buffered files...',verbose,h);
            tom_av2_stack_write(tmpstack,outfilestruct,used_part_idx);
            partcounter = 1;
            used_part_idx=[];
            stackcounter = stackcounter + file_buffer_size;
        end
    
        if (i==packages(ii,2))
            tom_emwritec(tmp_h_name,h_array,'subregion',[1 packages(ii,1) 1],[3 packages(ii,3) 1]);
       end;

    %disp(num2str(i));
    end; % end of for i

    
    
 end; % end of parfor

 
 %save align
 align2d=alignstruct; clear('alignstruct');
 h_array=tom_emread(tmp_h_name); h_array=h_array.Value;
 align2d=update_align2d(align2d,h_array);
 save(outalign_filename,'align2d');
 
 outfilestruct.size=[radius_bin*2 radius_bin*2 size(align2d,2)];
 save(outfilestruct.path,outfilestruct)
 
waitbar_wrapper(packages(1,2),packages(1,2)+2,'writing buffered files...',verbose,h);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Process Particle                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function particle=process_particle(particle,shift,angle,radius,rad_fact,binning,mask,norm_str)


%align particle
if (sum(shift==[0 0])~=2) || (angle~=0) 
     particle = tom_rotate(particle,angle);
     particle = tom_shift(particle,shift);
end

%cut out particle to final size
if (rad_fact > 1)
    sz_part=size(particle);
    try
        cut_start=round((sz_part-(2*radius/2^binning))./2);
        cut_stop=cut_start+(radius/2^binning)*2;
        particle = particle(cut_start(1):cut_stop(1)-1,cut_start(2):cut_stop(2)-1);
    catch
        particle = tom_taper(particle,[((radius/2^binning)*rad_fact*2) ((radius/2^binning)*rad_fact*2)]);
        particle = particle(radius/2^binning:3*radius/2^binning-1,radius/2^binning:3*radius/2^binning-1);
    end;
end;

%mask is not set ones if no mask ...SPEED!
if strcmp(mask,'no_mask')==0
    particle = particle.*mask;
end

if (size(particle,1) < 2*(radius./2^binning) || size(particle,2) < 2*(radius./2^binning) )
    particle=tom_taper(particle,[2*radius  2*radius]);
end;


%norm particle
if strcmp(norm_str,'no_norm')==0
    if strcmp(mask,'no_mask')==0
        particle = tom_norm((particle+5).*2,norm_str,mask);
    else
        particle = tom_norm((particle+5).*2,norm_str);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  check cutout                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inside_flag,start_cut,size_cut]=check_cutout(x,y,radius,rad_fact,imsize)

inside_flag=1;


start_cut(1)=x-(rad_fact*radius);
start_cut(2)=y-(rad_fact*radius);
start_cut(3)=1;

size_cut(1)=2.*(rad_fact*radius)-1;
size_cut(2)=2.*(rad_fact*radius)-1;
size_cut(3)=0;


if (start_cut(1) < 1)
    inside_flag=0;
    start_cut(1)=1;
end;

if (start_cut(2) < 1)
    inside_flag=0;
    start_cut(2)=1;
end;

if (x+rad_fact*radius > imsize(1))
    inside_flag=0;
    size_cut(1)=rad_fact*radius+(imsize(1)-x)-1;
end;

if ((y+rad_fact*radius> imsize(2)) )
    inside_flag=0;
    size_cut(2)=rad_fact*radius+(imsize(2)-y)-1;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Waitbar wrapper%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=waitbar_wrapper(max,progress,message,verbose,h)

if (verbose==0)
    return;
end;

if (progress > max)
    close(h);
    return;
end;

if (isempty(h))
    h = waitbar(progress./max,message);
end;


h = waitbar(progress./max,h,message);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update_h_array%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_array=update_h_array(h_array,count,iswritten,isaligned,normstr)
    

h_array(1,count)=iswritten;


if nargin < 4
    return; %hack for speed !
end;

h_array(2,count)=isaligned;

switch normstr
        case 'phase'
            normed_num=1;
        case '2std'
            normed_num=2;
        case '3std'
            normed_num=3;
        case 'mean0+1std'
            normed_num=4;
        otherwise
            error('norm method not implemeted');
end;


h_array(3,count)=normed_num;



function align2d=update_align2d(align2d,h_array)

for i=1:size(align2d,2)
    
    align2d(1,i).iswritten=h_array(1,i);
    align2d(1,i).isaligned=h_array(2,i);
    switch h_array(3,i)
        case 1
            hstr='phase';
        case 2
            hstr='2std';
        case 3
            hstr='3std';
        case 4
            hstr='mean0+1std';
        otherwise
            error('norm method not implemeted');
    end;
    align2d(1,i).normed=hstr;
    
    
end;





