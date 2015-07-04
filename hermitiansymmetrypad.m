function padded = hermitiansymmetrypad(z,dims)

n = size(z);
if(mod(dims(1),2)==0)
    redundant = z(2:(n(1)-1),:,:);
else
    redundant = z(2:end,:,:);
end;
redundantleft=redundant(:,1,:);
redundantright=fliplr(redundant(:,2:end,:));
redundant=horzcat(redundantleft,redundantright);
redundant=conj(redundant);
padded = vertcat(z, flipud(redundant));