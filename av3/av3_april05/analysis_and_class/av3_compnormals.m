function [mnangle mndev] = av3_compnormals(motlname,ind1,ind2,threshold,iclass)
% AV3_COMPNORMALS computes particles' normalvectors to determined
%       orientation 
%
%   [mnangle mndev] = av3_compnormals(motlname,ind1,ind2,threshold,iclass)
%
%   AV3_COMPNORMALS compares normalsvectors of particles according to their
%   MOTLS. The vectors of MOTLs number IND1 and IND2 are compared.
%
% PARAMETERS
%  INPUT
%   motlname        filename of corresponding motl - 'motlfilename'_#no.em
%   ind1            index of first motl to be compared
%   ind2            index of 2nd motl to be compared
%   threshold       threshold*mean(ccc) is cutoff for averaging 
%   iclass          class of particles
%
%  OUTPUT
%   mnangle         mean angle between normal IND1 and normal IND2
%   mndev           corresponding standard error
%
%   12/01/04 FF
%
% last change 03/31/05 FF - update docu

name = [motlname '_' num2str(ind1) '.em'];
motl = tom_emread(name);motl2 = motl.Value;
name = [motlname '_' num2str(ind2) '.em'];
motl = tom_emread(name);motl1 = motl.Value;
indx = find (motl2(1,:) > 0); meanv = mean(motl2(1,indx));
mnangle = 0;mndev=0;
icount =0;
for indpart = 1:size(motl2,2) 
    if (motl2(1,indpart)>threshold*meanv & motl2(20,indpart) == iclass)
        the1 = motl1(19,indpart);psi1=motl1(18,indpart);phi1=motl1(17,indpart);
        the2 = motl2(19,indpart);psi2=motl2(18,indpart);phi2=motl2(17,indpart);
        vec1 = tom_pointrotate([0 0 1],phi1,psi1,the1);
        vec2 = tom_pointrotate([0 0 1],phi2,psi2,the2);
        scpro = vec1*vec2';
        angl = 180*acos(scpro)/pi;
        mnangle = mnangle + angl;
        mndev = mndev+angl.^2;
        icount = icount+1;
    end;
end;
mnangle = mnangle/icount;
mndev = mndev / icount - mnangle.^2;
mndev = sqrt(mndev);