function [vol_weighted]=tom_exact_weight(vol,angles,thickness);
%TOM_EXACT_WEIGHT creates ...
%
%   [vol_weighted]=tom_exact_weight(vol,angles,thickness)
%
%PARAMETERS
%
%  INPUT
%   vol                 ...
%   angles              ...
%   thickness           ...
%   color               ...
%  
%  OUTPUT
%   vol_weighted		...
%
%       SUBROUTINE RANWEI (TILT,BR,NA,NB,IR,FR,MODE,LIST)
% C---------------------------------------------------------------------------
% C  MODE=3: A 3-dim. weighting function is calculated and stored in BR.
% C  MODE=4 (-4): A projection represented by its Fourier transform in BR
% C  is weighted directly. In both cases, the weighting corrects for the
% C  error arising from simple back projection. Based on SINC.
% C  The array TILT containes the tilt angles of all projections involved in 
% C  the 3-dim. reconstruction.
% C  Arguments:   IR(1) = number of column for PHI
% C               FR(2) = diameter of reconstruction volume
% C               FR(3) = THETA (In the case of THETA=0 it is assumed that 
% C                       THETA is given in the column IR(1)+1)
% C  7-Apr-1988,  R. Hegerl
% C  Update: 7-Jun-1990,  R. Hegerl
% C  Last update: 22-May-1995,  R. Hegerl (weighting of single projections)
% C---------------------------------------------------------------------------
%       INCLUDE '../common/emsys.f'
%       INCLUDE '../common/emfile.f'
% C
%       REAL      TILT(1), BR(1), FR(1)
%       INTEGER   IR(1), MODE, NA(1), NB(1)
%       LOGICAL   LIST
% C
% C  Local variables
%       REAL      A1, A2, ADD, ARG, CPHI, CTHE, PHI, SPHI, STHE, SUM, 
%      $          THETA, THICK, UST, VST, WIDTH, XST, YST, ZST, ZW
%       INTEGER   I, IND, IT, ITLT, IX, IY, IZ, NDIM
% C
% C  Check type and dimension of arrays
%       IF (NA(1).NE.4) CALL ERREND (63,*800,
%      $          'BAPRO(WEIGHT) * WRONG TYPE OF INPUT ARRAY')
%       IF (NA(2).LT.IR(1)+4) CALL ERREND (62,*800,
%      $          'BAPRO(WEIGHT) * X-DIMENSION OF INPUT ARRAY TOO SMALL')
%       IF (NB(1).NE.8) CALL ERREND (63,*800,
%      $          'BAPRO(WEIGHT) * WRONG TYPE OF OUTPUT ARRAY')
% C

% C  Define parameters
%       NDIM = (NB(2)-1)*2
%       IT = IR(1)
%       THICK = FR(2)
%       IF (FR(2).LE.0.0) THICK = NDIM
%       WIDTH = PI
%
%EXAMPLE
%   ... = tom_exact_weight(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

dimx=size(vol,1);
dimy=size(vol,2);
dimz=size(vol,3);

NDIM=dimx;
IT=1;
THICK=thickness;
if (thickness<=0) THICK=dimx; end;
WIDTH=pi; % cut off for sinc function, pi = first zero of sinc, extend to ...

BR=zeros(dimx./2.*dimy.*dimz,1);
% C
% C  Write header
%       IF (LIST .AND. MODE.NE.4) THEN
%         CALL TITL20 ('BAPRO(WEIGHT)$')
%         CALL TITLE1
%         CALL LINE_COUNT(2)
%         WRITE (NTOUT,901)
%  901    FORMAT (' TILT ANGLES USED FOR WEIGHTING FUNCTION:'/
%      $          '  NR.',8X,'PHI',7X,'THETA')
%       END IF
% C
%       IF (MODE.EQ.4) GO TO 100
% C
% C  For each projection, the projection vector is calculated from the
% C  angles PHI and THETA which are given in columns IT and IT+1 of the
% C  input array. The 3 components of this vector are stored in columns
% C  IT+2, IT+3, IT+4, thereby overwriting the previous content.

ITLT = 0;

%       DO   I = 1,NA(3)
%         IF (FR(3).NE.0.0) TILT(ITLT+IT+1) = FR(3)
%         IF (LIST) THEN
%           CALL LINE_COUNT(1)
%           WRITE (NTOUT,903) I,TILT(ITLT+IT),TILT(ITLT+IT+1)
%  903      FORMAT (I4,2(5X,F7.2))
%         END IF
%         PHI = TILT(ITLT+IT)*PI/180.
%         THETA = TILT(ITLT+IT+1)*PI/180.
%         ZW = SIN(THETA)*PI*THICK
%         TILT(ITLT+IT+2) = ZW*COS(PHI)
%         TILT(ITLT+IT+3) = ZW*SIN(PHI)
%         TILT(ITLT+IT+4) = COS(THETA)*PI*THICK
%         ITLT = ITLT + NA(2)
%       END DO
%       IF (MODE.EQ.-4) GO TO 100

dimangles=size(angles,1);
for I=1:dimangles
    PHI=angles(ITLT+IT,1).*pi./180;
    THETA=angles(ITLT+IT,2).*pi./180;
    ZW=sin(THETA).*pi.*THICK;
    TILT(I,1)=ZW.*cos(PHI);
    TILT(I,2)=ZW.*sin(PHI);
    TILT(I,3)=cos(THETA).*pi.*THICK;
    ITLT=ITLT+1;
end;

% C
% C  Main loop for calculating the weighting function for any pixel in 
% C  the 3-dim. output array.

      
%      DO 60   IZ = 1,NB(4)
%         YST = 0.0
%         DO 55   IY = 1,NB(3)
%           XST = 0.0
%           DO 50   IX = 1,NB(2)
%            SUM = 0.0
%            ITLT = IT + 2
%             DO 45   I = 1,NA(3)
%               ARG = XST*TILT(ITLT) + YST*TILT(ITLT+1) 
%      $              + ZST*TILT(ITLT+2)
%               IF (ABS(ARG).LT.1E-6) THEN
%                 ADD = 1.
%               ELSE IF (ABS(ARG).GE.1E-6 .AND. ABS(ARG).LT.WIDTH) THEN
%                 ADD = SIN(ARG)/ARG
%               ELSE
%                 ADD = 0.0
%               END IF
%               SUM = SUM + ADD
%   45          ITLT = ITLT + NA(2)
%             if (sum.gt.0.0) then
%               BR(IND) = 1./MAX(SUM,1.0)
%             else 
%               br(ind) = 0.0
%             end if
%             BR(IND+1) = 0.0
%             IND = IND + 2
%   50        XST = XST + 1./FLOAT(NDIM)
%           YST = YST + 1./FLOAT(NDIM)
%           IF (IY.EQ.NB(3)/2) YST = -YST
%   55    CONTINUE
% C        ZST = ZST + 1./FLOAT(NDIM)     ! 25-NOV-98  R. Hegerl
%         ZST = ZST + 1./FLOAT(NB(4))
%         IF (IZ.EQ.NB(4)/2) ZST = -ZST
%   60  CONTINUE
%       GO TO 800

IND = 1;
ZST = 0;
IT=0; 
for IZ=1:dimz
    YST=0;
          for IY=1:dimy
              XST=0;
              for IX=1:dimx./2
                  SUM=0;
                  ITLT=IT+1; 
                  for I=1:dimangles
                      ARG=XST.*TILT(ITLT,1)+YST.*TILT(ITLT,2)+ZST.*TILT(ITLT,3);
                      if abs(ARG)<1e-6
                          ADD=1;
                      elseif abs(ARG)>1e-6 & abs(ARG)<WIDTH
                              ADD=sin(ARG)./ARG;                          
                      else
                          ADD=0;
                      end;
                      SUM=SUM+ADD;
                      ITLT=ITLT+1;
                  end;
                  if (SUM>0) 
                      BR(IND)=1./SUM;
                      if SUM<1
                          BR(IND)=1;
                      end;
                  else
                      BR(IND)=0;
                  end;
                  IND=IND+1; % only real pointers !!! img ist set to 0 anyway
                  XST = XST + 1./(NDIM);
              end;
              YST = YST + 1./(NDIM);
              if IY==dimy./2
                  YST = -YST;
              end;
          end;
          ZST=ZST+1./dimz;
          if IZ==dimz./2
              ZST=-ZST;
          end;
      end;

      vol_weighted=reshape(BR,[dimx./2 dimy dimz]);
      vol_weighted=swap_weight(vol_weighted);


%       C
% C  Weighting of a single projection
%  100  IF (LIST) THEN
%         CALL LINE_COUNT(1)
%         WRITE (NTOUT,905) FR(5),FR(4)
%  905    FORMAT (' BAPRO(WEIGHT):    PROJECTION WITH    PHI=',F8.2,
%      $          ' DEG',6X,'THETA=',F8.2,' DEG')
%       END IF
%       SPHI = SIN(FR(5)*PI/180.)
%       CPHI = COS(FR(5)*PI/180.)
%       STHE = SIN(FR(4)*PI/180.)
%       CTHE = COS(FR(4)*PI/180.)
%       A1 = CTHE*CPHI
%       A2 = CTHE*SPHI
% C
% C  Main loop for weighting a single projection
%       IND = 1
%       VST = 0.0
%       DO 155   IY = 1,NB(3)
%         UST = 0.0
%         DO 150   IX = 1,NB(2)
%           SUM = 0.0
%           ITLT = IT + 2
%           XST = UST*A1 - VST*SPHI
%           YST = UST*A2 + VST*CPHI
%           ZST = -UST*STHE
%           DO 145   I = 1,NA(3)
%             ARG = XST*TILT(ITLT) + YST*TILT(ITLT+1) 
%      $            + ZST*TILT(ITLT+2)
%             IF (ABS(ARG).LT.1E-6) THEN
%               ADD = 1.
%             ELSE IF (ABS(ARG).GE.1E-6 .AND. ABS(ARG).LT.PI) THEN
%               ADD = SIN(ARG)/ARG
%             ELSE 
%               ADD = 0.0
%             END IF
%             SUM = SUM + ADD
%  145        ITLT = ITLT + NA(2)
%           BR(IND) = BR(IND)/SUM
%           BR(IND+1) = BR(IND+1)/SUM
%           IND = IND + 2
%  150      UST = UST + 1./FLOAT(NDIM)
%         VST = VST + 1./FLOAT(NDIM)
%         IF (IY.EQ.NB(3)/2) VST = -VST
%  155  CONTINUE
% C
%  800  RETURN
%       END


function volnewnewx=swap_weight(vol_weighted)

dimx=size(vol_weighted,1);
dimy=size(vol_weighted,2);
dimz=size(vol_weighted,3);


for i=1:dimz./2; volnew(:,:,i)=vol_weighted(:,:,dimz./2+i);end;
for i=1:dimz./2; volnew(:,:,i+dimz./2)=vol_weighted(:,:,i);end;
for i=1:dimy./2; volnewnew(:,i,:)=volnew(:,dimy./2+i,:);end;
for i=1:dimy./2; volnewnew(:,i+dimy./2,:)=volnew(:,i,:);end;
volnewnewx=zeros(dimx.*2,dimy,dimz);
ii=dimx;for i=2:dimx; volnewnewx(ii,:,:)=volnewnew(i,:,:);ii=ii-1;end;
for i=dimx+1:2*dimx; volnewnewx(i,:,:)=volnewnew(i-dimx,:,:);end;
