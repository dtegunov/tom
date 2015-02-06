function tom_mark_find_gradient(marker)

msize = size(marker);

s = 0;
n = 0;
ss = [];
while (n < 10000)
    
    
    idx = ceil(rand(1,15)*msize(2));
    idx = 1:msize(2);
    idx2 = ceil(rand(1,1) * msize(3));
    
    x = marker(1,idx,idx2);
    y = marker(2,idx,idx2);

    [m, centroid] = get_perpendicular_regression_line(x, y);
    
    drawnow;
    pause(0.0);
    
    if (~isnan(m))
        ss = [ss m];
        s = s+m; 
        n = n+1;
        disp(['new angle: ' num2str(m*180/pi) '; mean (' num2str(n) '): ' num2str(s/n*180/pi)]);
    end;
    
end;


1;


% y = mx + b
function [m, b] = get_regression_line(x, y)

n = length(x);

% slope = m
sy = sum(y);
sx = sum(x);

m = (n*sum(x.*y) - sx*sy) / (n*sum(x.^2) - sx*sx);

% intercept
b = (sy - m*sx) / n;


% http://www.mathpages.com/home/kmath110.htm
% line throgh point centroid with angle q (in rad)
function [q, centroid] = get_perpendicular_regression_line(x, y)

idx = ~(isnan(x) | isnan(y));
x = x(idx);
y = y(idx);

n = length(x);


%x = 1:15;
%x = x + (rand(1,length(x))-0.5)*30
%y = -2*x+20;
%y = y + (rand(1,length(x))-0.5)*30



centroid = [sum(x), sum(y)] / n;

x = x - centroid(1);
y = y - centroid(2);

xy = x.*y;
x2_y2 = x.^2 - y.^2;

A = x2_y2 / xy;

% tan(q)^2  +  A tan(q) - 1  =  0
%( -b +- sqrt(b^2 - 4ac) ) / 2a
q1 = sqrt(A^2 + 4);
q0 = atan( ( -A + q1 ) / 2 );
q1 = atan( ( -A - q1 ) / 2 );


% y'  =  -x sin(q) + y cos(q);
y0 = -x*sin(q0) + y*cos(q0);
y1 = -x*sin(q1) + y*cos(q1);

if ((sum(y0.^2) < sum(y1.^2)))
    q = q0;
else
    q = q1;
end;

%return;

if (q1 == q)
    q1 = q0;
end;

cla;
hold on;

axis equal;

plot(x+centroid(1), y+centroid(2), 'gx');

dimension = [0 0 2048 2048; 0 2048 0 2048];
plot(dimension(1,:),dimension(2,:),'r.');
r = max(dimension);
r = max(sqrt(sum(x.^2 + y.^2, 1)));

x__ =  [r*cos(q), r*cos(q+pi())] + centroid(1);
y__ =  [r*sin(q), r*sin(q+pi())] + centroid(2);
plot(x__, y__, 'b-');

x__ =  [r*cos(q1), r*cos(q1+pi())] + centroid(1);
y__ =  [r*sin(q1), r*sin(q1+pi())] + centroid(2);
plot(x__, y__, 'y-');




