function tom_equi_rotation(N)
global iii;
iii = 1;
N = 200;
dims = 4;
tic;
p = get_random_points_on_hyperspere(N, dims);
toc;
dist = compute_distance(p);

[mdist, mdistidx] = min(dist);

toc;
d2 = dist + dist';
d2(d2==0) = nan;

p2 = push_points(p);

idx = true(1, size(p,2));
for (i=1:dims)
    idx = idx & p(i,:) >= 0;
end;
%p = p(:, idx);

if (0)
for (i=1:size(p,2))
    dist = sum((repmat(p(:,i), [1,size(p,2)]) - p).^2, 1);
    if (sum(dist==0) ~= 1)
        error('point more times.');
    end;
    d_min(i) = min(dist(dist~=0));
    d_sum(i) = sum(dist(dist~=0));
end;    
end;

for (i=1:500)
    if (dims==3) 
        pp = p(:, p(1,:)>0 & p(2,:)>0 & p(3,:)>0);
        pp = p;
        figure(5); plot3(pp(1,:), pp(2,:), pp(3,:), 'g.'); axis equal;
        pp2 = p(:, p(3,:)>0);
        figure(4); plot(pp2(1,:), pp2(2,:), 'b.'); axis equal;
    elseif (0)
        pp4min = 0; pp4max=0.05; 
        pp = p(1:3,p(4,:) <= pp4max  & p(4,:)>=pp4min);
        figure(3); plot3(pp(1,:), pp(2,:), pp(3,:), 'g.'); axis equal;
        pause(0.01);
    end;
    
    tic;
    dists = sum(reshape(reshape(repmat(p, [N, 1]), [dims,N*N]) - repmat(p, [1,N]), [dims,N,N]) .^ 2, 1);
    sorted_dists = sort(dists, 2);
    dd = sort(sorted_dists(1,2,:));
    [dd(1), dd(end)]
    toc
    
    p = push_points(p);
end;
1
%figure(2); plot(sort(dist), 'b.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p2 = push_points(p)

[dims, N] = size(p);

p2 = nan(dims, N);
maxdist = 0;
for (i=1:N)
    vects = p - repmat(p(:,i), [1,N]);
    dist = sqrt(sum(vects.^2,1));
    delidx = find(~dist);
    vects(:,delidx) = [];
    dist(:,delidx) = [];
    maxdist = max(maxdist, min(dist));
    vects = vects ./ repmat(dist, [dims,1]);
    p2(1:dims, i) = sum(vects .* repmat(1./(abs(dist).^2), [dims,1]), 2);
end;


p3 = p - p2 / (N^4^(1/dims));
p4 = p3 ./ repmat(sqrt(sum(p3.^2, 1)), [dims,1]);

p2 = p4;

potential(p2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pot = potential(p)


for i=1:size(p,2)
    
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = compute_distance(p)

[dims, N] = size(p);
dist = zeros(N, N);
o = ones(N,1);
mo = -o;
for (i=1:N-1)
    %d = sqrt(sum((p - repmat(p(:,i), [1, N])) .^ 2, 1));
    d = acos(max(min(p(:,i+1:N)' * p(:,i), o(i+1:N)), mo(i+1:N)));
    dist(i, i:N) = [nan; d];
end;
dist(N,N) = nan;
dist = dist + dist';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = equi_unit_sphere_points(npts)

% see http://www.math.niu.edu/~rusin/known-math/95/sphere.faq

% Number of segments on the meridian
K = round(sqrt(npts * pi() / 4));



N = 1;
X(N) = 0;
Y(N) = 0;
Z(N) = 1;

for i=1:K-1
    M = round(2*K*sin(pi*(i/K)));    
    for j=0:M-1
        N = N+1;
        X(N) = sin(pi*(i/K))*cos(2*pi*(j/M));
        Y(N) = sin(pi*(i/K))*sin(2*pi*(j/M));
        Z(N) = cos(pi*(i/K));
    end;
end;
N = N + 1;
X(N) = 0;
Y(N) = 0;
Z(N) = -1;

pts = cat(1, X, Y, Z);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function returns an N by dims array, X, in which 
% each of the N rows has the dims Cartesian coordinates 
% of a random point uniformly-distributed over the 
% interior of an dims-dimensional hypersphere with 
% radius r and center at the origin.  The function 
% 'randn' is initially used to generate N sets of dims 
% random variables with independent multivariate 
% normal distribution, with mean 0 and variance 1.
% Then the incomplete gamma function, 'gammainc', 
% is used to map these points radially to fit in the 
% hypersphere of finite radius r with a uniform % spatial distribution.
% Roger Stafford - 12/23/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = get_random_points_inside_hyperspere(N,dims,radius)
 
X = randn(N,dims);
s2 = sum(X.^2,2);
X = X.*repmat(radius*(gammainc(s2/2,dims/2).^(1/dims))./sqrt(s2),1,dims);
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sphere = get_random_points_on_hyperspere(N, dims)

sphere = nan(dims, N);
nsphere = 0;
while (nsphere < N)
    s = randn(dims, N);
    dist = sqrt(sum(s.^2,1));
    idx = dist ~= 0;
    sidx = sum(idx);
    sphere(1:dims, nsphere+(1:sidx)) = s(1:dims, idx) ./ repmat(dist, [dims,1]);
    nsphere = nsphere + sidx;
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rot = quaternion2rot(q)

x = q(1);
y = q(2);
z = q(3);

x2 = x * x;
y2 = y * y;
z2 = z * z;
xy = x * y;
xz = x * z;
yz = y * z;
wx = w * x;
wy = w * y;
wz = w * z;

rot = [[ 1 - 2*(y2 + z2),     2*(xy - wz),     2*(xz + wy), 0 ]; ...
       [     2*(xy + wz), 1 - 2*(x2 + z2),     2*(yz - wx), 0 ]; ...
       [     2*(xz - wy),     2*(yz + wx), 1 - 2*(x2 + y2), 0 ]; ...
       [               0,               0,               0, 1 ]];




    
    