function [surface_out count]=tom_remove_triangles(surface_in,thresh)

s=reshape(surface_in.faces,[3.*size(surface_in.faces,1) 1]);

[a]=unique(s);
%a=int32(a);
%b=int32(b);
%c=int32(c);

size(a,1)
count=zeros(size(a,1),1);
size(count,1)
parfor i=1:size(a,1)
    count(i)=sum(s-a(i)==0);
    if ~mod(i,1000)
        disp('.');drawnow;
    end;
end;
parfor i=1:
   
surface_out=surface_in;
surface_out.faces=new_faces;
