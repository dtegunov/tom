function tom_display(in,flag)
%shows 1d,2d,3d and 4d singal 
%rider's info 
% 4d pseudo time domain can also be encoded as cell-index
% 3d and 4d are chopped down 2 2d slices

if  (nargin < 2)
    flag='standard';
end;


if (iscell(in))
    tmp=in{1};
    stack=zeros(size(tmp,1),size(tmp,2),size(tmp,3).*length(in));
    zz=1;
    for i=1:length(in)
        tmp=in{i};
        stack(:,:,zz:(zz-1)+size(tmp,3))=tmp;
        zz=zz+size(tmp,3);
    end
    im_tmp=tom_gallery(stack,[length(in) size(tmp,3)]);
    clear('stack');
    figure; tom_imagesc(im_tmp);
    return;
end;

if (min(size(in))==1)
    figure; plot(in);
    return;
end;

if (length(size(in))==3)
    figure; tom_dspcub(in);
    return;
end;

if (length(size(in))==2)
    if (strcmp(flag,'standard'))
        figure; tom_imagesc(in);
    end;
    if (strcmp(flag,'single_plot'))
        
        if (size(in,2) > 20)
            l_im=20;
        else
            l_im=size(in,2);
        end;
        figure;
        d1=ceil(sqrt(l_im));
        d2=ceil(sqrt(l_im));
        for i=1:l_im
             subplot(l_im,i,1);
             plot(in(i,:),1);
        end;
    
    end;
    if (strcmp(flag,'surface'))
        figure; tom_imagesc(in);
    end;
    return;
end;

if (length(size(in))==4)
    stack=zeros(size(in,1),size(in,2),size(in,3).*size(in,4));
    zz=1;
    for i=1:size(in,4)
        stack(:,:,zz:(zz-1)+size(in,3))=in(:,:,:,i);
        zz=zz+size(in,3);
    end
    im_tmp=tom_gallery(stack,[size(in,3) size(in,4)]);
    clear('stack');
    figure; tom_imagesc(im_tmp);
    return;
end;
