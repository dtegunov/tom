function [stack_out ang_st]=tom_reshape_stack(filename_in,in_struct,binning,filename_out,normflag,maskstruct)
%TOM_RESHAPE_STACK reshape 3d volumes
%
%   stack_out=tom_reshape_stack(filename_in, in_struct, binning, filename_out, normflag, maskstruct)
%
% tom_reshape_stack('parts/','stack_1',1)
%
%PARAMETERS
%
%  INPUT
%   filename_in         ...
%   in_struct           ...
%   binning             ...
%   filename_out        ...
%   normflag            ...
%   maskstruct          ...
%  
%  OUTPUT
%   stack_out   		...
%
%EXAMPLE
%   ... = tom_reshape_stack(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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


ang_st=[0 0 0 0]';

error(nargchk(0, 6, nargin, 'struct'))

dimension = [];

if ~exist('normflag','var')
    normflag = 0;
end



if (iscell(filename_in))
   spi_flag=0;
    if (tom_isemfile(filename_in{1}))
         h=tom_reademheader(filename_in{1});
    else
        spi_flag=1;
        tmpp_sp=tom_spiderread(filename_in{1});
        h=tom_emheader(tmpp_sp.Value);
    end;
   flag_3d=1;
else
    h=tom_reademheader(filename_in);
    flag_3d=0;
end;


if (nargin==1)
    if (flag_3d==1)
        dimension(2)=h.Header.Size(1).*h.Header.Size(2).*h.Header.Size(3);
        dimension(1)=size(filename_in,2);    
    else
        dimension(2)=h.Header.Size(1).*h.Header.Size(2);
        dimension(1)=h.Header.Size(3);      
    end;
    binning=0;
    processing='memory';
end;

if (nargin==2)
    binning=0;
    processing='memory';
end;


if (nargin==3)
    processing='memory';
end;

if (nargin==4)
    stack='';
    processing='memory';
end;

if (nargin==5)
    processing='memory';
end

if (nargin==6)
    processing='memory';
end


if (isempty(dimension))
    if (flag_3d==1)
        dimension(2)=h.Header.Size(1).*h.Header.Size(2).*h.Header.Size(3);
        dimension(1)=size(filename_in,2);    
    else
        dimension(2)=h.Header.Size(1).*h.Header.Size(2);
        dimension(1)=h.Header.Size(3);      
    end;
end;


if (isempty(binning))
    binning=0;
end;

if normflag == 1
    if flag_3d==1
        maskstruct.mask.Value = h.Header.Size./2^binning;
        %mask=tom_create_mask(maskstruct.mask);    
        mask = tom_sphere([20 20 20],8,2);
    else
        maskstruct.mask.Value = [h.Header.Size(1)./2^binning, h.Header.Size(2)./2^binning];
        mask=tom_create_mask(maskstruct.mask);    
    end
end

if flag_3d==1
    maskstruct.rotmask.Value = h.Header.Size;
    maskstruct.rotmask.Apply=2;
    maskstruct.rotmask.Value.size_x=h.Header.Size(1);
    maskstruct.rotmask.Value.size_y=h.Header.Size(2);
    maskstruct.rotmask.Value.size_z=h.Header.Size(3);
    rotmask=tom_create_mask(maskstruct.rotmask);
%    rotmask = tom_sphere([128 80 80],35,2);
end


if strcmp(processing,'memory')

    if (flag_3d==1)
        %3d particle stack
        dimension(2)=dimension(2)/(8.^binning);        
        stack_out=zeros(round(dimension));
        pause(1);
        for i=1:dimension(1)
            if (tom_isemfile(filename_in{i})==1 || tom_isspiderfile(filename_in{i})==1 )
                
                if (spi_flag==0)
                    im=tom_emreadc(filename_in{i});
                else
                    im=tom_spiderread(filename_in{i});
                end;
                
                if isempty(find(isnan(im.Value)))~=1
                    disp(['NaN found at:' num2str(i)]);                    
                end;

                shiftx = in_struct(i).Shift.X;
                shifty = in_struct(i).Shift.Y;
                shiftz = in_struct(i).Shift.Z;
                rotmatrix = in_struct(i).Angle.Rotmatrix;
                wege_size_tmp=90-(((in_struct(i).Tomogram.AngleMin)+abs(in_struct(i).Tomogram.AngleMax))./2);
                ang_st(:,i)=[in_struct(i).Angle.Phi in_struct(i).Angle.Psi in_struct(i).Angle.Theta wege_size_tmp]';
                im = im.Value;
                if (isempty(rotmatrix))
                    rotmatrix=[0 0 0];
                end;
                im = tom_rotate(im,rotmatrix,'linear');
                im = tom_shift(double(im),[shiftx shifty shiftz]);
                im=im.* rotmask;
                im=tom_bin(im,binning);
                if normflag == 1
                    im = tom_norm(im,'phase',mask);
                end
                %hack
               % im(20:32,:,:)=0;
                stack_out(i,:)=reshape(im,1,[]);
            end;
        end;
    else
        %2d particle stack
        sz=h.Header.Size;
        dimension(2)=dimension(2)/(4.^binning);
        stack_out=zeros(round(dimension));
        for i=1:dimension(1)
            im=tom_emreadc(filename_in,'subregion',[1 1 i],[sz(1)-1 sz(2)-1 0]);
            im=tom_bin(im.Value,binning);
            if normflag == 1
               im = tom_norm(im,'phase',mask);
            end
            stack_out(i,:)=reshape(im,1,[]);
        end;
    end;


else

    % Harddisk ... to be done


    if (isdir(filename_in)==1)
        %3d particle stack
        laufx=0;
        d=dir(filename_in);
        for lauf=3:(size(d,1))
            if (tom_isemfile([filename_in d(lauf).name])==1)
                h=tom_reademheader([filename_in d(lauf).name]);
                laufx=laufx+1;
            end;
        end;
        h=tom_reademheader([filename_in d(3).name]);
        sz=h.Header.Size./(2.^binning);
        st=zeros(laufx,sz(1).*sz(2).*sz(3));
        laufx=0;
        for lauf=3:(size(d,1))-1 %dirty hack!!!!
            if (tom_isemfile([filename_in d(lauf).name])==1)
                laufx=laufx+1
                im=tom_emread([filename_in d(lauf).name]);
                im.Value=tom_bin(im.Value,binning);
                sz=im.Header.Size./(2.^binning);
                st(laufx,:)=reshape(im.Value,1,sz(1).*sz(2).*sz(3));
            end;
        end;
    else
        %2d particle stack
        h=tom_reademheader(filename_in);
        sz=h.Header.Size;
        st=zeros(sz(3),sz(1).*sz(2));
        %tom_emwritec(filename_out,[sz(3) (sz(1).*sz(2))],'new');

        for i=1:sz(3)
            im=tom_emreadc(filename_in,'subregion',[1 1 i],[sz(1)-1 sz(2)-1 0]);
            st(i,:)=reshape(im.Value,1,sz(1).*sz(2));
            %tom_emwritec(filename_out,im,'subregion',[i 1 1],[1 (sz(1).*sz(2)-1) 0]);
        end;
    end;

    tom_emwrite(filename_out,st);

end;