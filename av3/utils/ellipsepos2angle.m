function motl = ellipsepos2angle(motl,r_0,R,C,iseries)
% ELLIPSEPOS2ANGLE determines Euler angles from ellipse parameters
%
%   motl = ellipsepos2angle(motl,r_0,R,C,iseries);
%
%   Assume particles are on located on an ellipse and project normal to
%   this surface (e.g. membrane). The resulting Euler angles psi and 
%   theta are calculated and stored in motl.
%
%   Ellipse is assumed to have the form:
%
%  ( (x - x_0)/R )^2 + ( (y - y_0)/R )^2 + ( (z - z_0)/C )^2 = 1
%
% PARAMETERS
%  INPUT
%   MOTL            motive list
%   R_0             center of ellipse [x_0, y_0, z_0]
%   R               xy-radius of ellipse
%   C               z-radius of ellipse
%   ISERIES         index of tilt series (column 5 in MOTL)
%
%  OUTPUT
%   MOTL            motive list with new PSIs and THETAs
%
% SEE ALSO
%   ELLIPSEDET, NORMALVEC

for ipart =1:size(motl,2);
    if motl(5,ipart) == iseries
        x = motl(8,ipart) - r_0(1);
        y = motl(9,ipart) - r_0(2);
        z = motl(10,ipart) - r_0(3);
        xnor = x/(R^2);
        ynor = y/(R^2);
        znor = z/(C^2);
        psi = 90+180/pi*atan2(ynor,xnor);
        theta = atan2(sqrt(xnor.^2+ynor.^2),znor)*180/pi;
        motl(18,ipart) = psi;
        motl(19,ipart) = theta;
    end;
end;