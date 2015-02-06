function tom_Rec3dDirectCut(vol,picklist,boxsize,particlename)

NumberOfParticles = size(picklist,1);

for p=1:NumberOfParticles

particle(1:boxsize,1:boxsize,1:boxsize) = ...
    vol(picklist(p,1)-boxsize/2+1:picklist(p,1)+boxsize/2, ...
          picklist(p,2)-boxsize/2+1:picklist(p,2)+boxsize/2, ...
          picklist(p,3)-boxsize/2+1:picklist(p,3)+boxsize/2);

particle = tom_emheader(particle);

tom_emwrite([particlename '_' num2str(p) '.vol'],particle);

clear particle;

end
