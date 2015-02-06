function angles=tom_create_angles_st(tiltaxis,tiltangle)

z=1;

if (tiltaxis(1)==0 & tiltaxis(2)==0 & tiltaxis(3)==0)
   tiltaxis(1)=0;
   tiltaxis(2)=1;
   tiltaxis(3)=0;
end;

for j=tiltaxis(1):tiltaxis(2):tiltaxis(3)

    for i=tiltangle(1):tiltangle(2):tiltangle(3)
        
        angles(1,z)=j;
        angles(2,z)=i;
        z=z+1;
    end

end;
