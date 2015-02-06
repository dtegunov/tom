function av3_motlanalyze(motlfilename,first,n)
% AV3_MOTLANALYZE plots convergence of MOTL
%
%   av3_motlanalyze(motlfilename,first,n)
%
%   analyze convergence of motl graphically. The rms of determined shifts
%   and angles respective to the previous iteration are plotted.
%
%   motlfilename    name of motls - expected as 'motlfilename'_#.em
%   first           number # of first motl to be analyzed
%   n               number # of last motl to be analyzed
%
%   FF 08/11/03

for ind=first:n
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    if ind > first
        dist(ind,:,:) = motl.Value - prevmotl;
    end;
    prevmotl = motl.Value;
end;
figshifts = figure;
sqx = sum(sqrt(dist(first+1:n,11,:).^2),3)/size(prevmotl,2);
sqy = sum(sqrt(dist(first+1:n,12,:).^2),3)/size(prevmotl,2);
sqz = sum(sqrt(dist(first+1:n,13,:).^2),3)/size(prevmotl,2);
plot(first+1:n,squeeze(sqx),'r');
hold on;
plot(first+1:n,squeeze(sqy),'b');
plot(first+1:n,squeeze(sqz),'m');
title('RMS of difference of Shifts');
legend('x-shift','y-shift','z-shift');
figangles = figure;
sqphi = sum(sqrt(dist(first+1:n,17,:).^2),3)/size(prevmotl,2);
sqpsi = sum(sqrt(dist(first+1:n,18,:).^2),3)/size(prevmotl,2);
sqthe = sum(sqrt(dist(first+1:n,19,:).^2),3)/size(prevmotl,2);
plot(first+1:n,squeeze(sqphi),'r');
hold on;
plot(first+1:n,squeeze(sqpsi),'b');
plot(first+1:n,squeeze(sqthe),'m');
title('RMS of difference of angles');
legend('\phi','\psi','\theta');