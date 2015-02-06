function E = GauGau(a)
% GAUGAU error of data and Double Gaussian
%
%   E = GauGau(a)
%   fit to
%       Y = a(1).*exp(-((x-a(2))/a(3)).^2) + a(4).*exp(-((x)/a(5)).^2);
%   start values for two symmetrical Gauss functions

global xData yData dy
x=xData(:); y=yData(:); 
Y = a(1).*exp(-((x-a(2))/a(3)).^2) + a(4).*exp(-((x)/a(5)).^2);
E = sum(((y-Y).^2)./dy.^2);

