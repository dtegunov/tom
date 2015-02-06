function tom_av2_createstack(inalignname, classes, outstackname, outalignname, alignedflag, normflag, radius)
%TOM_AV2_CREATESTACK creates a particle stack from an alignment file
%
%   tom_av2_createstack(inalignname, classes, outstackname, outalignname, alignedflag, normflag, radius)
%   
%
%PARAMETERS 
%  IN
%  inalignname      full path to input alignment file (.mat)
%  classes          cell of classnames to export, set to '' if all classes should be exported
%  outstackname     full path to output stack file (.em)
%  outalignname     full path to output alignment file (.mat)
%  alignedflag      1: save particles aligned, 0: save unaligned (optional)
%  normflag         0: no norming, 1: phase norming, 2: 2std norming (optional)
%  radius           radius of the particles in the stack (optional)
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%    Copyright (c) 2006
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom
%
%   14/02/06 AK

if nargin < 6
    normflag = 0;
end
if nargin < 5
    alignedflag = 0;
end

try
    s = load(inalignname);
catch
    error('Could not load alignment file.');
end

alignstruct = s.align2d;
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
stack = zeros(radius*2,radius*2,stacksize);
selectedparticles = [];

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
        header = tom_reademheader(alignstruct(1,i).filename);
        imagesize = header.Header.Size;
        if x <= radius*2 | y <= radius*2 | x > imagesize(1)-radius*2-1 | y > imagesize(2)-radius*2-1

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
            part_box = tom_emreadc([alignstruct(1,i).filename],'subregion',[lowx lowy 1],[highx-lowx highy-lowy 0]);
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
        if alignstruct(1,i).angle ~= 0 & alignedflag == 1
            particle = tom_rotate(particle,alignstruct(1,i).angle);
        end
        %shift particle
        if alignedflag == 1
            particle = tom_shift(particle,[alignstruct(1,i).shift.x alignstruct(1,i).shift.y]);
        end

        %reduce to final size
        particle = particle(radius:3*radius-1,radius:3*radius-1);
        %norm particle
        if normflag == 1
            particle = tom_norm(particle,'phase');
            alignstruct(1,i).normed='phase';
        elseif normflag == 2
            particle = tom_norm(particle,'3std');
            alignstruct(1,i).normed='3std';
        end
        %apply to stack
        stack(:,:,i) = particle;
        selectedparticles = [selectedparticles, i];
    end

    waitbar(i./numpoints,h);
    i = i + 1;
end

tom_emwrite(outstackname,stack);
align2d = alignstruct(1,selectedparticles);
save(outalignname,'align2d');
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