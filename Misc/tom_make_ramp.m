function img=tom_make_ramp(dim)


img=ones(dim,dim);
power=log(1024)./log(2);

area=2.^[0:power-1];

for i=1:2:size(area,2)-1
    img(area(i):area(i+1)-1,1:dim./2)=-1;
    
    
end;
img(1:dim,dim./2+1:end)=(img(1:dim,1:dim./2).*1);