function [Residuum,x3d,y3d,z3d,tx,ty,alpha,...
             isoscale,xscale,yscale,zscale] = Marker3dRigidBodyREF(handles)

%--------------------------------------------------------------------------
% MARKER3D Rigid Body REF
%--------------------------------------------------------------------------
disp(' ');
disp('------------------------------------------------------------------');
disp('Welcome to Marker3d Rigid Body REF');
disp('------------------------------------------------------------------');
disp(' ');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% DATA
%--------------------------------------------------------------------------
global v s xm ym theta Imdim

v = handles.v;
s = handles.s;
xm = handles.xm;
ym = handles.ym;
theta = handles.theta;
Imdim = handles.Imdim;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% CONTROL
%--------------------------------------------------------------------------
global o a

% Which origin
o = [1 handles.rm handles.xrm handles.yrm handles.zrm];

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
x3d0 = handles.x3dalg3d+(Imdim/2+1);
y3d0 = handles.y3dalg3d+(Imdim/2+1);
z3d0 = handles.z3dalg3d+(Imdim/2+1);
tx0 = handles.txalg3d;
ty0 = handles.tyalg3d;
alpha0 = handles.alphaalg3d(handles.rp);

StartValues = [x3d0 y3d0 z3d0 tx0 ty0 alpha0];
%--------------------------------------------------------------------------
clear handles;

%--------------------------------------------------------------------------
% OPTIMIZATION
%--------------------------------------------------------------------------
if (a==0)
RigidBodyREFresults = fminsearch(@RigidBodyREFfun,StartValues,options);
end
if (a==1)
RigidBodyREFresults = lsqnonlin(@RigidBodyREFfun,StartValues,lb,ub,options);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% RESULTS
%--------------------------------------------------------------------------

% Tilt axis azimuth
alpha(1:v) = RigidBodyREFresults(3.*s+2.*v+1);

% 3d Marker coordinates
for j=1:s
    mx(j) = RigidBodyREFresults(j+0.*s);
    my(j) = RigidBodyREFresults(j+1.*s);
    mz(j) = RigidBodyREFresults(j+2.*s);
end

x3d = mx-o(3);
y3d = my-o(4);
z3d = mz-o(5);

% Translations
for i=1:v
    tx(i) = RigidBodyREFresults(i+3.*s+0.*v);
    ty(i) = RigidBodyREFresults(i+3.*s+1.*v);
end

% Scale
isoscale(1:v) = 1;
xscale(1:v) = 1;
yscale(1:v) = 1;
zscale(1:v) = 1;

% Residuum
Residuum = 0;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% MODEL FUNCTION
%--------------------------------------------------------------------------
function RigidBodyREFvalues = RigidBodyREFfun(x)

% Constant real data
global v s xm ym theta Imdim
% Control variables
global o a

% Marker3d Rigid Body REF
zaehler = 1;

for j=1:s
    for i=1:v
       
    if isnan(xm(i,j)) == 1  % Missing data
        continue
    end

    % ROTATION AND PROJECTION OPERATOR
    P = [(sin(x(3.*s+2.*v+1)))^2*cos(theta(i))+(cos(x(3.*s+2.*v+1)))^2 ...
          sin(x(3.*s+2.*v+1))*cos(x(3.*s+2.*v+1))*(1-cos(theta(i))) ...
          sin(x(3.*s+2.*v+1))*sin(theta(i)); ...
          sin(x(3.*s+2.*v+1))*cos(x(3.*s+2.*v+1))*(1-cos(theta(i))) ...
         (cos(x(3.*s+2.*v+1)))^2*cos(theta(i))+(sin(x(3.*s+2.*v+1)))^2 ...
         -cos(x(3.*s+2.*v+1))*sin(theta(i))];
    
     % Optimization with IMAGE DIMENSION
         if (j==o(2)) % Constraint
         MinimizeThis = [xm(i,j);ym(i,j)]-...
                 (P*[o(3)-(Imdim/2+1);o(4)-(Imdim/2+1);o(5)-(Imdim/2+1)]+...
                 [x(i+3.*s+0.*v)+(Imdim/2+1);x(i+3.*s+1.*v)+(Imdim/2+1)]);
         end   
        
         if (j~=o(2)) % Free equations 
         MinimizeThis = [xm(i,j);ym(i,j)]-...
                 (P*[x(j+0.*s)-(Imdim/2+1);x(j+1.*s)-(Imdim/2+1);x(j+2.*s)-(Imdim/2+1)]+...
                 [x(i+3.*s+0.*v)+(Imdim/2+1);x(i+3.*s+1.*v)+(Imdim/2+1)]);
         end
 
         
    RigidBodyREFvalues(zaehler) = MinimizeThis(1);
    zaehler = zaehler + 1;
    
    RigidBodyREFvalues(zaehler) = MinimizeThis(2);
    zaehler = zaehler + 1;
    
    clear MinimizeThis;
     
    end
end


if (a==0)  % Simplex Search Method
    
    RigidBodyREFvaluesSum = 0;
    
    for k=1:(zaehler-1)
        RigidBodyREFvaluesSum = RigidBodyREFvaluesSum + (RigidBodyREFvalues(k))^2;
    end
    
    RigidBodyREFvalues = RigidBodyREFvaluesSum;

end

if (a==1)  % Conjugate gradients
    RigidBodyREFvalues = RigidBodyREFvalues;
end
%--------------------------------------------------------------------------
