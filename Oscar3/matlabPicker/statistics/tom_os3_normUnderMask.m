function [normed_vol mea st n]= tom_os3_normUnderMask(volume,mask)

mask=mask~=0;

n=sum(sum(sum(mask~=0)));
mea=sum(sum(sum((volume.*mask).*(volume~=0))))./n;
st=sqrt(sum(sum(sum((((mask==0).*mea)+(volume.*mask) -mea).^2)))./n);
normed_vol=((volume-mea)./st).*mask;