%%
function [res area]= tom_os3_getLocalPeak(img,templ,pos,angle)

%cut out area around the peak and recorrelate with the rotated template
%returns peak volume for specified template rotation

area = tom_os3_subVolume(pos,img,templ,'center');

%%
templ = tom_rotate(templ,angle,'linear');

res = tom_os3_corr(area,templ);

