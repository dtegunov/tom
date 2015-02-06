function tom_absinth

% tom_absinth simulates the effect of 80% absinth
% on JP and FB
%
% SN, 21/12/05

figure;

l=tom_emread('lenna.em');

set(gcf,'Position',[164   231   648   594]);
set(gcf,'DoubleBuffer','on')

for i=1:4:65; 
    imagesc(tom_filter(l.Value,i)');colormap hot
title('cheers ... ')
    drawnow;
end;

for i=65:-4:1; 
    imagesc(tom_filter(l.Value,i)');colormap cool
title('four hours later ... ')
    drawnow;
end;
