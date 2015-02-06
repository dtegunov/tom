function [val x y z] = tom_os3_max3d(img)

res = zeros(size(img,3),1);
a = zeros(size(img,3),1);
b = zeros(size(img,3),1);

for i=1:size(img,3)
    [res(i) a(i) b(i) ] = tom_os3_max2d(img(:,:,i));
end;


[val i] = max(res);
x = a(i);
y = b(i);
z = i;


