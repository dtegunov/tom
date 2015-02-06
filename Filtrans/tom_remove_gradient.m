function im=tom_remove_gradient(im)

im_sz=size(im);

st_sz=round(im_sz(1)./15);

im=(im+1).*2;

background = imopen(im,strel('disk',st_sz));

im=im-background;
