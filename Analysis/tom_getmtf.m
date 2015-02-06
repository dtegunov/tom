function [estimated_lambda,y] = tom_getmtf(filename, filterval, split, binning, lowcutoff);
warning off;
if nargin < 5
    lowcutoff = 100;
end

if nargin < 4
    binning =1;
end
if nargin < 3
    split = 1;
end
if nargin < 2
    filterval = 1;
end

if ischar(filename)
    try
        im = tom_emreadc(filename,'binning',binning);
    catch
        error('Could not open file.');
    end

    lowcutoff = lowcutoff ./ 2^binning ./ (split./2);
    
    if filterval > 1
        im.Value = tom_filter(im.Value,filterval,'quadr','real');
    end

    %calculate and integrate power spectrum
    ps = tom_norm(tom_image2ctf(double(im.Value),1,0,split,0),1);
end

f = figure;
t = (1:length(ps))';
y = ps';
plot(t,y,'r'); hold on; h = plot(t,y,'b'); hold off;
title('Input data'); ylim([0 1]);
start = [1;0;0;0;0];
% We use an anonymous function to pass additional parameters t, y, h to the
% output function.
outputFcn = @(x,optimvalues,state) fitoutputfun(x,optimvalues,state,t,y,h);
options = optimset('OutputFcn',outputFcn,'TolX',0.01);
estimated_lambda = fminsearch(@(x)fitfun(x,t,y),start,options);

close(f);
warning on;









function err = fitfun(lambda,t,y)
%FITFUN Used by FITDEMO.
%   FITFUN(lambda,t,y) returns the error between the data and the values
%   computed by the current function of lambda.
%
%   FITFUN assumes a function of the form
%
%     y =  c(1)*exp(-lambda(1)*t) + ... + c(n)*exp(-lambda(n)*t)
%
%   with n linear parameters and n nonlinear parameters.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.8.4.1 $  $Date: 2004/11/29 23:30:50 $

A = zeros(length(t),length(lambda));
for j = 1:length(lambda)
   A(:,j) = exp(-lambda(j)*t);
end
c = A\y;
z = A*c;
err = norm(z-y);