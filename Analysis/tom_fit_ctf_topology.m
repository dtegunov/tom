function [Fit Dz_det Dz_delta_det statistic x y z]=tom_fit_ctf_topology(img,Fit,EM,Search,region_size,ctf_size,step,topology)

%tom_fit_ctf_topology fits a defocus plane
%   
%
%
%PARAMETERS
%
%  INPUT
%   img                 micrograph
%   Fit, EM, Search     Parameter Structures from tom_fit_ctf_gui
%   region_size         region
%   ctf_size            size of ctf
%   step                step size
%   topology            'plane', only plane implemented, yet
%
%  OUTPUT
%  
%EXAMPLE
%     
% Example:
% [Fit Dz_det Dz_delta_det statistic x y z]==tom_fit_ctf_topology(img,Fit,EM,Search,[2048],[256],[256],'plane');
%
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by SN 04/07/10
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
%


% Example:
% [Dz_det Dz_delta_det statistic]=tom_fit_ctf_topology(img,Fit,EM,Search,[1024 1024],[256 256],[256]);
% 
if nargin<8
    topology='plane';
end;

Search_tmp=Search;

il=0;
iil=0;
idx=0;

for i=1:step:size(img,1)-region_size
    il=il+1;
    for ii=1:step:size(img,2)-region_size
        idx=idx+1;
        iil=iil+1;
        
        region=img(i:i+region_size-1,ii:ii+region_size-1);
        x(idx,1)=round((i+region_size(1)./2));
        y(idx,1)=round((ii+region_size(1)./2));
        
        ps=tom_calc_periodogram(double(region),ctf_size(1));
        ps=(log(fftshift(ps)));
        
        inner_rad_bg=8;
        outer_rad_bg=(ctf_size)/2-1;
                        
        [decay decay_image]=calc_decay(ps,Search.mask_inner_radius_bg,Search.mask_outer_radius_bg,16);
        background_corrected_ps=double(ps-decay_image);        
        
        mask_in_radius=Search.mask_inner_radius;
        mask_out_radius=Search.mask_outer_radius;
        ps_size=size(ps);
        mask_in = tom_spheremask(ones(ps_size),mask_in_radius,0,[ps_size(1)./2+1 ps_size(2)./2+1 1]);
        mask_out = tom_spheremask(ones(ps_size),mask_out_radius,0,[ps_size(1)./2+1 ps_size(2)./2+1 1]);
        mask=mask_out-mask_in;
        
        
        
        ps=background_corrected_ps.*mask;
%        tom_imagesc(ps); drawnow;
        warning off;
        [Fit]=tom_fit_ctf(ps,EM,Search_tmp);
        warning on;
%        Fit.Dz_det
        Dz_det(il,iil)=Fit.Dz_det;
        z(idx,1)=Fit.Dz_det;
        Dz_delta_det(il,iil)=Fit.Dz_delta_det;
        [a1 a2 a3 a4 ]=tom_dev(Fit.corr_all,'noinfo');
        statistic.cc_std(il,iil)=a4;
        statistic.cc_mean(il,iil)=a1;
        
        
    end;
    iil=0;
end;

% fit a 2d plane
[a z_fit residuals]=fit_plane_hesse(x,y,z);
% write in Fit structure
Fit.Plane.Pos_x=x;
Fit.Plane.Pos_y=y;
Fit.Plane.Dz_det=z;
Fit.Plane.Dz_fit=z_fit;
Fit.Plane.residuals=residuals;
Fit.Plane.Hesse_coeff=a;
Fit.Plane.Hesse_normal_form='0 = a(1) + a(2)*x+a(3)*y+a(4)*z , (a(2),a(3),a(4)) is the normal to the plane, abs(a(1)) its distance from (0,0,0)';

function [a z_fit residuals]=fit_plane_hesse(x,y,z)

N=length(x);
A=[ones(N,1),x,y,z];
[U,S,V]=svd(A,0);
sps=diag(S);
[SSmin,Imin]=min(sps);
a=V(:,Imin);
a=a/norm(a(2:4),2);
z_fit=-((a(1) + a(2)*x+a(3)*y)./a(4));
residuals=z_fit-z;

% the Hesse -normal - form of the plane is
%0 = a(1) + a(2)*x+a(3)*y+a(4)*z , (a(2),a(3),a(4)) is the normal to the plane
%and abs(a(1)) its distance from (0,0,0)


function fit_plane_PCA(X,vis)


maxlim = max(abs(X(:)))*1.1;
[coeff,score,roots] = princomp(X);
basis = coeff(:,1:2);
normal = coeff(:,3);
pctExplained = roots' ./ sum(roots);
[n,p] = size(X);
meanX = mean(X,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
residuals = X - Xfit;
error = abs((X - repmat(meanX,n,1))*normal);
sse = sum(error.^2);

if vis==1
    [xgrid,ygrid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),5), ...
        linspace(min(X(:,2)),max(X(:,2)),5));
    zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
    h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);
    
    hold on
    above = (X-repmat(meanX,n,1))*normal > 0;
    below = ~above;
    nabove = sum(above);
    X1 = [X(above,1) Xfit(above,1) nan*ones(nabove,1)];
    X2 = [X(above,2) Xfit(above,2) nan*ones(nabove,1)];
    X3 = [X(above,3) Xfit(above,3) nan*ones(nabove,1)];
    plot3(X1',X2',X3','-', X(above,1),X(above,2),X(above,3),'o', 'Color',[0 .7 0]);
    nbelow = sum(below);
    X1 = [X(below,1) Xfit(below,1) nan*ones(nbelow,1)];
    X2 = [X(below,2) Xfit(below,2) nan*ones(nbelow,1)];
    X3 = [X(below,3) Xfit(below,3) nan*ones(nbelow,1)];
    plot3(X1',X2',X3','-', X(below,1),X(below,2),X(below,3),'o', 'Color',[1 0 0]);
    
    hold off
    maxlim = max(abs(X(:)))*1.1;
    axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
    axis square
    view(-23.5,5);
end;


