function [envavg,noiseavg] = tom_mtf_createenvelope(e,n,x,plotflag)

if nargin < 4
    plotflag = 1;
end

if nargin < 3
    x = 512;
end

if plotflag == 1
    figure;
    hold on;
end
x=[1:x];

noisefunction = zeros(size(n,1),x(end));
envfunction = zeros(size(e,1),x(end));

for i=1:size(e,1)

   noisefunction(i,:) = noiseplotfun(x',n(i,:));
   envfunction(i,:) = envplotfun(x',e(i,:));
   
   if plotflag == 1
       plot(noisefunction(i,:));
       plot(envfunction(i,:)+noisefunction(i,:),'r');
   end
   
end

noiseavg = median(noisefunction,1);
envavg = median(envfunction,1) + noiseavg;

if plotflag == 1
    figure;hold on;
    plot(envavg,'r');
    plot(noiseavg);
    figure;hold on;

    for i=1:size(e,1)
        envfunction(i,:) = envfunction(i,:)+noiseavg;
    end
    
    envstd = std(envfunction,0,1);
    lower = envavg - envstd;
    upper = envavg + envstd;
    tom_ciplot(lower,upper,x,'r');

    noisestd = std(noisefunction,0,1);
    lower = noiseavg - noisestd;
    upper = noiseavg + noisestd;
    tom_ciplot(lower,upper,x,'b');
end

noiseavg = squeeze(noiseavg);
envavg = squeeze(envavg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise function                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = noiseplotfun(x,n)
for i=1:size(x,1)
    y(i) = (n(1)+n(2).*(double(i).^(1/2))+n(3).*double(i)+n(4).*double(i).^2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% envelope function                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = envplotfun(x,n)
for i=1:size(x,1)
    y(i) = (n(1)+n(2).*(double(i).^(1/2))+n(3).*double(i)+n(4).*double(i).^2+n(5).*double(i).^3);
end