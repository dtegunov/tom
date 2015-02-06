function res = tom_os3_line(source,template,angle)
dim = tom_os3_fftSelect(source);
centerX = size(template,1)/2;
centerY = size(template,2)/2;

if(length(angle) == 1)
    angle = abs(360 - angle);
else
    angle = abs([360 360 360] - angle);
end;

source = tom_rotate(source,angle,'cspline','taper');


res.xLine = source(1:size(template,1),centerY);
res.yLine = source(centerX,1:size(template,2));
res.zLine = 0;



