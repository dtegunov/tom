function E = Gau(a)
% GAU error of data and Gaussian
%
%   E = Gau(a)
%
%   fit to
%       Y = a(1).*exp(-((x-a(2)).^2)/(2*(a(3).^2)));
%   start values for two symmetrical Gauss functions
%   
%   last change: 
%   12/20/04 FF - fixed sigma bug

global xData yData dy
x=xData(:); y=yData(:); 
Y = a(1).*exp(-((x-a(2)).^2)/(2*(a(3).^2)));
E = sum(((y-Y).^2)./dy.^2);

