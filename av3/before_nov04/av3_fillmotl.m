function motivelist = av3_fillmotl(peaklist, angfile, phi_end, psi_end, theta_end, angular_incr, phi_start, psi_start, theta_start)
% AV3_FILLMOTL creates MOTL from PEAKLIST
%
%motivelist = av3_fillmotl(peaklist, angfile, phi_end, psi_end, theta_end, angular_incr, phi_start, psi_start, theta_start)
%
%   peaklist:       x   |   y   |   z  |   CCC   
%

error(nargchk(6,9,nargin));
if (nargin < 7)
    phi_start = 0;
end;
if (nargin < 8)
    psi_start = 0;
end;
if (nargin < 9)
    theta_start = 0;
end;
n = size(peaklist,2);
motivelist = zeros(20,n);
for ind = 1:n,
    motivelist(1,ind) = peaklist(4,ind);
    x = peaklist(1,ind);y = peaklist(2,ind);z = peaklist(3,ind);
    motivelist(8,ind) = x;motivelist(9,ind) = y;motivelist(10,ind) = z;
    motivelist(2,ind) = x;motivelist(3,ind) = y;
    angleindex = tom_emread(angfile, 'subregion', [x y z],[0 0 0]);
    [phi, psi, theta] = av3_index2angle(angleindex.Value, phi_end, psi_end, theta_end, angular_incr, ...
            phi_start, psi_start, theta_start);
    disp([num2str(ind) '. Peak angles   = ' num2str(phi) ' ' num2str(psi) ' ' num2str(theta) ])
    motivelist(17,ind) = phi;
    motivelist(18,ind) = psi;
    motivelist(19,ind) = theta;
end;