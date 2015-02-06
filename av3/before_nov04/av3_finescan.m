npart = size(motl,2);
for i=1:npart,
    particle = tom_red(volume, [x-radius, y-radius, z-radius], [2*radius+1 2*radius+1 2*radius+1]);
    x=motl(8,i);
    y=motl(9,i);
    z=motl(10,i);
    particle = tom_spheremask(particle, radius,smooth);
end;
