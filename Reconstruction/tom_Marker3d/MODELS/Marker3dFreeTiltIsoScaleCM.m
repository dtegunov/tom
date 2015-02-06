function [Residuum,x3d,y3d,z3d,tx,ty,alpha,...
                                          isoscale,xscale,yscale,zscale] = Marker3dFreeTiltIsoScaleCM(handles)

%--------------------------------------------------------------------------
% MARKER3D Free Tilt IsoScale CM
%--------------------------------------------------------------------------
disp(' ');
disp('------------------------------------------------------------------');
disp('Welcome to Marker3d Free Tilt IsoScale CM');
disp('------------------------------------------------------------------');
disp(' ');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% DATA
%--------------------------------------------------------------------------
global v s xm ym theta

v = handles.v;
s = handles.s;
xm = handles.xm;
ym = handles.ym;
theta = handles.theta;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% CONTROL
%--------------------------------------------------------------------------
global a

% Which algorithm
if strcmp(handles.Algorithm,'Conjugate gradients') == 1
    a = [1];
end
if strcmp(handles.Algorithm,'Simplex search') == 1
    a = [0];
end

% Which precision and optimization options
options = optimset;
options.Display = 'iter';
options.LargeScale = 'on';
options.MaxFunEvals = 50000;
options.MaxIter = 50000;
options.TolFun = str2double(handles.Precision);
options.TolX = str2double(handles.Precision);
lb = -inf;
ub = +inf;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% START VALUES
%--------------------------------------------------------------------------

% Tilt Axis Azimuth
alpha0 = handles.alphaalg3d(handles.rp);

% Marker
x3d0 = handles.x3dalg3d;
y3d0 = handles.y3dalg3d;
z3d0 = handles.z3dalg3d;

marker = [cos(alpha0) sin(alpha0) 0; -sin(alpha0) cos(alpha0) 0; 0 0 1]*...
                 [x3d0 - mean(x3d0); y3d0 - mean(y3d0); z3d0 - mean(z3d0)];
     
x3d0 = marker(1,:);
y3d0 = marker(2,:);
z3d0 = marker(3,:);

% Translation
for i=1:v
   
    zaehler = 0;
    
    tsumx=0;
    tsumy=0;
    
    for j=1:s 
      
        if isnan(xm(i,j)) == 1  % Missing data
            zaehler=zaehler+1;
            continue
        end
       
    R = [cos(alpha0) -sin(alpha0); sin(alpha0) cos(alpha0)];
    P = [1 0 0; 0 cos(theta(i)) -sin(theta(i))];
     
    translation0 = [xm(i,j);ym(i,j)]-R*P*[x3d0(j);y3d0(j);z3d0(j)];
                   
    tsumx = tsumx + translation0(1);
    tsumy = tsumy + translation0(2);

    end

    tx0(i)=(1/(s-zaehler))*tsumx;
    ty0(i)=(1/(s-zaehler))*tsumy; 

end

alpha0(1:v) = alpha0;
isoscale0(1:v) = 1.0;
StartValues = [x3d0 y3d0 z3d0 tx0 ty0 alpha0 isoscale0];
%--------------------------------------------------------------------------
clear handles;

%--------------------------------------------------------------------------
% OPTIMIZATION
%--------------------------------------------------------------------------
if (a==0)
FreeTiltIsoScaleCMresults = fminsearch(@FreeTiltIsoScaleCMfun,StartValues,options);
end
if (a==1)
FreeTiltIsoScaleCMresults = lsqnonlin(@FreeTiltIsoScaleCMfun,StartValues,lb,ub,options);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% RESULTS
%--------------------------------------------------------------------------

% Tilt axis azimuth
for i=1:v
alpha(i) = FreeTiltIsoScaleCMresults(3.*s+2.*v+i);
end
alphaMean = mean(alpha);

% 3d Marker coordinates
for j=1:s
    mx(j) = FreeTiltIsoScaleCMresults(j+0.*s);
    my(j) = FreeTiltIsoScaleCMresults(j+1.*s);
    mz(j) = FreeTiltIsoScaleCMresults(j+2.*s);
end

marker = [cos(alphaMean) -sin(alphaMean) 0; sin(alphaMean) cos(alphaMean) 0; 0 0 1]*...
                [mx;my;mz];
    
x3d = marker(1,:);
y3d = marker(2,:);
z3d = marker(3,:);

% Translations
for i=1:v
    tx(i) = FreeTiltIsoScaleCMresults(i+3.*s+0.*v);
    ty(i) = FreeTiltIsoScaleCMresults(i+3.*s+1.*v);
end

% Scale
for i=1:v
isoscale(i) = FreeTiltIsoScaleCMresults(3.*s+3.*v+i);
end
xscale(1:v) = 1;
yscale(1:v) = 1;
zscale(1:v) = 1;

% Residuum
Residuum = 0;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% MODEL FUNCTION
%--------------------------------------------------------------------------
function FreeTiltIsoScaleCMvalues = FreeTiltIsoScaleCMfun(x)

% Constant real data
global v s xm ym theta
% Control variables
global a

% Marker3d FreeTilt IsoScale CM

xcm = 0;  % Constraint
ycm = 0; 
zcm = 0;

    for j=2:s
        xcm = xcm + x(j+0.*s); 
        ycm = ycm + x(j+1.*s);
        zcm = zcm + x(j+2.*s);
    end

zaehler = 1;

for j=1:s % Cost function
    for i=1:v
       
    if isnan(xm(i,j)) == 1  % Missing data
        continue
    end

    % Rotation operator R and projection operator P
    R = [cos(x(3.*s+2.*v+i)) -sin(x(3.*s+2.*v+i)); ...
         sin(x(3.*s+2.*v+i)) cos(x(3.*s+2.*v+i))];
    P = [1 0 0;0 cos(theta(i)) -sin(theta(i))];
    
     
         if (j==1) % Constraint
         MinimizeThis = [xm(i,j);ym(i,j)]-(R*P*(x(3.*s+3.*v+i))*[-xcm;-ycm;-zcm]+...
              [x(i+3.*s+0.*v);x(i+3.*s+1.*v)]);
         end   
        
         if (j~=1) % Free equations 
         MinimizeThis = [xm(i,j);ym(i,j)]-(R*P*(x(3.*s+3.*v+i))*[x(j+0.*s);x(j+1.*s);x(j+2.*s)]+...
              [x(i+3.*s+0.*v);x(i+3.*s+1.*v)]);
         end
 
         
    FreeTiltIsoScaleCMvalues(zaehler) = MinimizeThis(1);
    zaehler = zaehler + 1;
    
    FreeTiltIsoScaleCMvalues(zaehler) = MinimizeThis(2);
    zaehler = zaehler + 1;
    
    clear MinimizeThis;
     
    end
end


if (a==0)  % Simplex Search Method
    
   FreeTiltIsoScaleCMvaluesSum = 0;
    
    for k=1:(zaehler-1)
        FreeTiltIsoScaleCMvaluesSum = FreeTiltIsoScaleCMvaluesSum + ...
                                                            (FreeTiltIsoScaleCMvalues(k))^2;
    end
    
    FreeTiltIsoScaleCMvalues = FreeTiltIsoScaleCMvaluesSum;

end

if (a==1)  % Conjugate gradients
    FreeTiltIsoScaleCMvalues = FreeTiltIsoScaleCMvalues;
end
%--------------------------------------------------------------------------
