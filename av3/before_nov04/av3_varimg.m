function varvol = av3_varimg(ref,motl,particlename,iclass)

% AV3_VARIMG creates variance volume of set of particles
%
%   varvol = av3_varimg(ref,motl,particlename,iclass)
%
%   Variance map at each voxel is calculated with respect 
%   to a reference. 
%   
%   REF             3D volume - make sure dims of ref are the 
%                     same as individual particles' dims
%   MOTL            motivelist - 20 x no.part array
%   PARTICLENAME    Filename are expected to be <particlename>_#ipart.em
%   ICLASS          Class of particles - if not specified, all particles 
%                   are taken
%
%   VARVOL          Volume containg the variances at each voxel
%
% FF 08/29/03
% last change 11/05/04 FF

error(nargchk(3,4,nargin))
if (nargin < 4)
    iclass = -1;
end;
varvol = ref*0;
for indx=1:size(motl,2)
    ipart = motl(4,indx);
    jclass = motl(20,indx);
    if ((jclass == iclass) | (iclass == -1))
        name = [particlename '_' num2str(ipart) '.em'];
        part = tom_emread(name);part = part.Value;
        x=motl(14,indx);y=motl(15,indx);z=motl(16,indx);
        part = tom_shift(part,[-x -y -z]);
        phi = motl(17,indx);
        psi = motl(18,indx);
        the = motl(19,indx);
        part = double(tom_rotate(part,[-psi,-phi,-the]));
        tom_dspcub(part);
        tvar = (part-ref).^2;
        varvol = varvol + tvar;
        varvol = sqrt(varvol/size(motl,2));
        disp(['processed particle ' num2str(ipart) ' ...']);
    end;
end;
