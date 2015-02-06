function stack = tom_av2_createstack(inalignname, classes, outstackname, outalignname, alignedflag, normflag, radius)
%TOM_AV2_CREATESTACK creates a particle stack from an alignment file
%
%   stack = tom_av2_createstack(inalignname, classes, outstackname, outalignname, alignedflag, normflag, radius)
%
%PARAMETERS
%
%  INPUT
%   inalignname         full path to input alignment file (.mat)
%   classes             cell of classnames to export, set to '' if all classes should be exported
%   outstackname        full path to output stack file (.em)
%   outalignname        full path to output alignment file (.mat)
%   alignedflag         1: save particles aligned, 0: save unaligned (optional)
%   normflag            0: no norming, 1: phase norming, 2: 2std norming (optional)
%   radius              radius of the particles in the stack (optional)
%  
%  OUTPUT
%   stack               ...
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%   cre0ted by AK 02/14/06
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


if nargin < 6
    normflag = 0;
end
if nargin < 5
    alignedflag = 0;
end

if ~isstruct(inalignname)
    try
        s = load(inalignname);
    catch
        error('Could not load alignment file.');
    end

    alignstruct = s.align2d;
else
    alignstruct = inalignname;
end

numpoints = size(alignstruct,2);


if strcmp(classes,'')
    [classes] = get_classes(alignstruct);
end

%determine radius
if nargin < 7
    [xxxxxx,radii] = get_classes(alignstruct);
    radius = max(cell2mat(radii));
end



%determine stack size
stacksize = 0;
for i=1:numpoints
    if strmatch(alignstruct(1,i).class,classes,'exact')
            stacksize = stacksize + 1;
    end
end
if stacksize == 0
    error('Particle stack would contain no particles!');
end

h = waitbar(0,'Creating particle stack, please be patient...');
tom_emwritec(outstackname,[radius*2,radius*2,stacksize],'new','single');
%stack = zeros(radius*2,radius*2,stacksize,'single');
selectedparticles = zeros(stacksize,1);

lauf = 1;

partcounter = 1;
stackcounter = 1;
tmpstack = zeros(alignstruct(1,1).radius*2,alignstruct(1,1).radius*2,5000,'single');

try
    
    for i=1:numpoints
        
        
        
        %add particle to stack if classname matches
        if strmatch(alignstruct(1,i).class,classes,'exact');
            %classname = storage_av2_particlepicker.align(1,j).class;
            x = alignstruct(1,i).position.x;
            y = alignstruct(1,i).position.y;
            
            %aligned particles
            if alignedflag == 1
                
                %particle is already aligned
                if alignstruct(1,i).isaligned == 0
                    %FIXME: align remaining unaligned particles
                end
                
            end
            
            %cut particles
            %taper mode
            
            %alignstruct(1,i).filename =  strrep(alignstruct(1,i).filename,'high','high_corr');
            
            try
                header = tom_reademheader(alignstruct(1,i).filename);
            catch
                warning(['error reading file ' alignstruct(1,i).filename]);
                continue;
            end
            imagesize = header.Header.Size;
            if x <= radius*2 || y <= radius*2 || x > imagesize(1)-radius*2-1 || y > imagesize(2)-radius*2-1
                
                lowx = x-2*radius;
                lowy = y-2*radius;
                highx = x+2*radius-1;
                highy = y+2*radius-1;
                
                %set x or y to 1 if particle is in the left or upper edge
                if x <= radius*2
                    lowx = 1;
                end
                if y <= radius*2
                    lowy = 1;
                end
                
                %set x or y to size of image if particle is in the right or lower edge
                if x > imagesize(1)-radius*2-1
                    highx = imagesize(1);
                end
                if y > imagesize(2)-radius*2-1
                    highy = imagesize(2);
                end
                
                %cut out particle, this will give a non quadratic matrix
                try
                    part_box = tom_emreadc2([alignstruct(1,i).filename],'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
                catch
                    continue;
                end
                part_box = single(part_box.Value);
                
                %taper in x direction
                if size(part_box,1) < radius*4
                    if lowx == 1
                        stripe = part_box(1,:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,stripe,part_box);
                        end
                    else
                        stripe = part_box(size(part_box,1),:);
                        while size(part_box,1) < radius*4
                            part_box = cat(1,part_box,stripe);
                        end
                    end
                end
                %taper in y direction
                if size(part_box,2) < radius*4
                    if lowy == 1
                        stripe = part_box(:,1);
                        while size(part_box,2) < radius*4
                            part_box = cat(2,stripe,part_box);
                        end
                    else
                        stripe = part_box(:,size(part_box,2));
                        while size(part_box,2) < radius*4
                            part_box = cat(2,part_box,stripe);
                        end
                    end
                end
                
                particle = part_box;
            else
                particle = tom_emreadc(alignstruct(1,i).filename,'subregion',[x-2*radius y-2*radius 1],[4*radius-1 4*radius-1 0]);
                particle = single(particle.Value);
            end
            
            %rotate if particle is aligned
            if alignedflag == 1 && alignstruct(1,i).angle ~= 0
                particle = tom_rotate(particle,alignstruct(1,i).angle);
            end
            %shift particle
            if alignedflag == 1
                try
                    particle = tom_shift(particle,[alignstruct(1,i).shift(1) alignstruct(1,i).shift(2)]);
                catch
                    particle = tom_shift(particle,[alignstruct(1,i).shift.x alignstruct(1,i).shift.y]);
                end
            end
            
            %reduce to final size
            particle = particle(radius:3*radius-1,radius:3*radius-1);
            %norm particle
            if normflag == 2
                particle = tom_norm(particle+100000,'phase');
                alignstruct(1,i).normed='phase';
            elseif normflag == 3
                particle = tom_norm(particle,'3std');
                alignstruct(1,i).normed='3std';
            end
            %apply to stack
            
            tmpstack(:,:,partcounter) = particle;
            partcounter = partcounter + 1;
            if partcounter > 5000
                partcounter = 1;
                tom_emwritec(outstackname,tmpstack,'subregion',[1 1 stackcounter],[size(tmpstack,1) size(tmpstack,2) size(tmpstack,3)]);
                stackcounter = stackcounter + 5000;
            end
            
            if i==numpoints
                tmpstack = tmpstack(:,:,1:partcounter-1);
                tom_emwritec(outstackname,tmpstack,'subregion',[1 1 stackcounter],[size(tmpstack,1) size(tmpstack,2) size(tmpstack,3)]);
            end
            
            %        tom_emwritec(outstackname,particle,'subregion',[1 1 lauf],[size(particle,1) size(particle,2) 1]);
            %stack(:,:,i) = particle;
            selectedparticles(lauf) = i;
            lauf = lauf + 1;
        end
        
        waitbar(i./numpoints,h);
        %    i = i + 1;
    end
    
catch
    disp([alignstruct(1,i).filename ' ' num2str(i)]);
end;

%if ~isempty(outstackname)
    %tom_emwrite(outstackname,stack);
%end
if ~isempty(outalignname)
    align2d = alignstruct(1,selectedparticles>0);
    save(outalignname,'align2d');
end

close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get classes in alignment file                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newclasses,newradii] = get_classes(align)
newclasses = {};
newradii = {};
for i=1:size(align,2)
    if strmatch(align(1,i).class,newclasses,'exact')
    else
        newclasses{size(newclasses,2)+1} = align(1,i).class;
        newradii{size(newradii,2)+1} = align(1,i).radius;
    end
end