function vol=tom_average_ps(path,ext,start,stop,offset,window_size,split,pre_binning,post_binning)
%TOM_AVERAGE_PS creates ...
%
%   vol=tom_average_ps(path,ext,start,stop,offset,window_size,split,pre_binning,post_binning)
%
%PARAMETERS
%
%  INPUT
%   path                ...
%   ext                 ...
%   start               ...
%   stop                ...
%   offset              ...
%   window_size         ...
%   split               ...
%   pre_binning         ...
%   post_binning        ...
%  
%  OUTPUT
%   vol         		...
%
%EXAMPLE
%   ... = tom_average_ps(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

if (split==0)
    split=1;
end;

if (split~=1 & split~=4 & split~=16 & split~=64)
    error('split must be 4 16 64 ');
end;

%allocate Memory
im=tom_emread([path num2str(start) ext]);
size_im=(size(im.Value)./(2^pre_binning))./split;
z=floor(( ((stop-start)-(window_size-1))./offset));
if (z<1)
    z=1;
end;

vol=zeros((size_im(1)./(2^post_binning)),(size_im(1)./(2^post_binning)),z);
h=waitbar(0,'averaging');

 tmp=zeros(size_im);
 if (window_size>=(stop-start))
     last_file=start;
 else
     last_file=(stop-window_size-1);
 end;

 if (window_size==1)
    last_file=stop;
 end;
 z=1;

for i=start:offset:last_file

    disp(['Processing File: ' path num2str(i) ext]);
    waitbar((i./stop));

    for ii=1:(window_size)
        disp(['           Processing File Buffer: ' [path num2str(i+ii-1) ext] ]);
        try
            im=tom_emreadc([path num2str(i+ii-1) ext]);
       catch
            disp(['         ERROR: File Buffer: ' path num2str(i+ii-1) ext] );
            continue;
        end;
        im=tom_bin(im.Value,pre_binning);
        for iii=1:split
            disp(['                         Processing Split ' num2str(iii) ]);
            if (split>1)
                if (iii==1) 
                    im_old=im;
                end;
                im=split_image(im_old,split,iii);
             end;
            im=tom_smooth(im,round(size_im(1).*0.1));
            im=tom_ps(double(im));
            tmp=tmp+im;
        end;
    end;

    vol(:,:,z)=tom_bin(tmp,post_binning);
    tmp=zeros(size_im);
    z=z+1;
end;

try
    close(h);
catch
    disp('do not close the wait bar');
end;



function split_image=split_image(im,number_of_splits,split_nr)

im_sz=size(im,1);

inkre=round(im_sz./number_of_splits);
start=((split_nr-1).*inkre)+1;
stop=((split_nr).*inkre);
split_image=im(start:stop,start:stop);
    


    


