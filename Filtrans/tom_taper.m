function  out=tom_taper(in,new_size)
%TOM_TAPER performs a 2d or 3d tapering
%
%   out=tom_taper(in,new_size)
%
%PARAMETERS
%
%  INPUT
%   in                  image or volume
%   new_size            new size of image or volume
%
%  
%  OUTPUT
%   out                 tapered image or volume
%
%EXAMPLE
%   im=tom_emread('pyrodictium_1.em');
%   out=tom_taper(tom_bin(im.Value,1),[2048 2048]);
%
%REFERENCES
%
%SEE ALSO
%   tom_smooth
%
%   created by FB 10/31/05
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


if (nargin==2)
    integr=0;
end;

if (size(new_size,2)==2)
    new_size(3)=1;
end;

%allocate some memory
out=zeros(new_size);
out=tom_paste(out,in,[round( ((size(out)-size(in))./2 )+1 ) ]);

    
for z=1:new_size(3)


    diff_z=round((new_size(3)-size(in,3))./2);

    if (z>diff_z & z<=(size(in,3)+diff_z)  )
        im=in(:,:,z-diff_z);
        out_sl=out(:,:,z);

        a=im(:,1);
        b=im(:,size(im,2));
        c=im(1,:);
        d=im(size(im,1),:);

        diff_hor=round((new_size(1)-size(im,1))./2);
        diff_vert=round((new_size(2)-size(im,2))./2);
        stop_up=diff_vert;
        start_down=size(out,2)-diff_vert;
        for i=1:new_size(1)
           if (i<=diff_hor | i > (size(im,1)+diff_hor)  )
                if (i<diff_hor )
                    val_up=a(1);
                    val_low=b(1);
                else
                    val_up=a(size(a,2));
                    val_low=b(size(b,2));
                end;
            else
                val_up=a(i-diff_hor);
                val_low=b(i-diff_hor);
            end;

            out_sl(i,1:stop_up)=val_up;
            out_sl(i,start_down:size(out,2))=val_low;
        end;
        stop_left=diff_hor;
        start_right=size(out,1)-diff_hor;
        for i=1:new_size(2)
            if (i>diff_vert & i <= (size(im,2)+diff_vert)  )
                out_sl(1:stop_left,i)=c(i-diff_vert);
                out_sl(start_right:size(out,1),i)=d(i-diff_vert);
            end;
        end;
        out(:,:,z)=out_sl;
    end;
end;


if (new_size(3)==1)
    return
end;



for z=1:new_size(3)


    diff_z=round((new_size(3)-size(in,3))./2);

    if (z<=diff_z  |  z>=(size(in,3)+diff_z)  )
        if (z<=diff_z)
            out(:,:,z)=out(:,:,diff_z+1);
        else
            out(:,:,z)=out(:,:,diff_z+size(in,3));
        end;
    end;
    
end;


    






















