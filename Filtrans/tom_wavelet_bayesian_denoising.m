function Id=tom_wavelet_bayesian_denoising(I, allowed_scale, SNR0, SNRF)
% Usage: Id=tom_wavelet_bayesian_denoising(I, allowed_scale, SNR0, SNRF)
%
% This function applies the a bayesian denoising procedure
% to all scales equal or greater than allowed_scale
%
% Typical use:
%    I=im2double(imread(original_image));
%    Id=xmipp_wavelet_bayesian_denoising(I,3,0.1,0.2);
%
% Reference:
% C.O.S.Sorzano, E. Ortiz, M. López, J. Rodrigo.
%    Improved Bayesian image denoising based on wavelets with
%    applications to Electron Microscopy.
%    Pattern Recognition, 39: 1205-1213 (2006) 

global C;
global P;
global Ncoefs;
global Neval;

allowed_scale=log2(size(I,1))-2-allowed_scale;

verbose=0;
if verbose
    disp('Computing Wavelet transform ...');
end
qmf = MakeONFilter('Daubechies',12);
WId=FWT2_PO(I,0,qmf);
powerI=sum(sum(I.^2));
Ncoefs_total=size(I,1)*size(I,2);

% Estimation of the power at each band ------------------------------------
if verbose
    disp('Evaluating power at each band ...');
end
max_scale=log(length(WId))/log(2);
scale=max_scale-1:-1:allowed_scale+1;
power=zeros(size(scale));
average=zeros(size(scale));
Ncoefs=zeros(size(scale));
for i=1:length(scale)
    sI=[xmipp_DWT_coefs(WId,scale(i),1);
        xmipp_DWT_coefs(WId,scale(i),2);
        xmipp_DWT_coefs(WId,scale(i),3)];
    power(i)=sum(sum(sI.^2));
    average(i)=mean(sI);
    Ncoefs(i)=length(sI);
end

% Evaluate the power of the unconsidered part of the image
sI=xmipp_DWT_coefs(WId,scale(length(scale)),0);
power_rest=sum(sum(sI.^2));
Ncoefs_rest=length(sI);

% Prepare the equation system
alpha=0.01;
C=[eye(length(scale)) eye(length(scale))];
P=(power./Ncoefs)';
A=zeros(2*(length(scale)-1),2*length(scale));
for i=2:length(scale)
    A(i-1,i-1)=1;
    A(i-1,i  )=-1;
    A(i-1+(length(scale)-1),i-1+(length(scale)))=1;
    A(i-1+(length(scale)-1),i  +(length(scale)))=-1;
end
A=[A; -eye(2*length(scale))];
b=zeros(size(A,1),1);

% Power constraints
aux0Ncoefs=SNR0*Ncoefs;
auxFNcoefs=SNRF*Ncoefs;
A=[A; -auxFNcoefs  Ncoefs;
       aux0Ncoefs -Ncoefs];
b=[b; 0; 0];

% Total power restriction
Aeq=[Ncoefs  Ncoefs];
beq=[powerI-power_rest];

if verbose
    disp('Decomposing power at each band ...');
end
[estimatedS,RESNORM,RESIDUAL,EXITFLAG]=...
    lsqlin(C,P,A,b,Aeq,beq,[],[],[],...
    optimset('LargeScale','off','Display','off'));

Neval=0;
estimatedS=fmincon(@xmipp_bayesian_goal, estimatedS,A,b,Aeq,beq,[],[],...
    @xmipp_bayesian_nonlinear_constraints,...
    optimset('LargeScale','off','GradObj','on','TolFun',1e-10,'Display','off'));

% Bayesian denoising ------------------------------------------------------
if verbose
    disp('Denoising ...');
end
for i=1:length(scale)
    N=estimatedS(i);
    S=estimatedS(i+length(scale));
    
    % Denoise
    WId=xmipp_bayesian_denoising_LL(WId, scale(i), 1, 0, S, N);
    WId=xmipp_bayesian_denoising_LL(WId, scale(i), 2, 0, S, N);
    WId=xmipp_bayesian_denoising_LL(WId, scale(i), 3, 0, S, N);
end

% Compute the inverse image
if verbose
    disp('Computing Inverse Wavelet transform ...');
end
Id=IWT2_PO(WId,0,qmf);

% Show information
if verbose
    disp(['Noise  power at each scale=' num2str(estimatedS(1:length(scale))')]);
    disp(['Signal power at each scale=' num2str(estimatedS((length(scale)+1:2*length(scale)))')]);
    disp(['Local SNR=' num2str([...
       estimatedS(length(scale)+1:2*length(scale))./estimatedS(1:length(scale))]')]);
end
tSNR=zeros(length(scale),1); % Cumulated (total) SNR
for i=1:length(scale)
   tSNR(i)=sum(Ncoefs(1:i).*estimatedS(length(scale)+1:length(scale)+i)')/...
           sum(Ncoefs(1:i).*estimatedS(1:i)');
end
if verbose
    disp(['Cumulated SNR=' num2str(tSNR')]);
end
end

% ---------------------------------------------------------------------
function [f,g,H]=xmipp_bayesian_goal(x)
    global C;
    global P;
    global Neval;
    
    % f=1/2*(Cx-P)'*(Cx-P)=1/2*x'C'Cx-P'Cx+P'P
    % g=C'Cx-C'P
    % H=C'C
    
	error=C*x-P;
	f=0.5*sum(error.*error);
    H=C'*C;
    g=H*x-C'*P;
    Neval=Neval+1;
end

% ---------------------------------------------------------------------
function [c,ceq]=xmipp_bayesian_nonlinear_constraints(x)
   global Ncoefs;
   ceq=[];
   
   scale=length(Ncoefs);
   tSNR=zeros(scale,1); % Cumulated (total) SNR
   lSNR=zeros(scale,1); % Local SNR
   c=zeros(2*scale-1,1);
   c(2*scale-1)=1;
   if sum(x(1:scale))>0
       for i=1:scale
           tSNR(i)=sum(Ncoefs(1:i).*x(scale+1:scale+i)')/...
                   sum(Ncoefs(1:i).*x(1:i)');
       end
       lSNR=x(scale+1:2*scale)./x(1:scale);
       c(1:scale-1)=lSNR(1:scale-1)-lSNR(2:scale);
       c(scale:2*scale-2)=tSNR(1:scale-1)-tSNR(2:scale);
       c(2*scale-1)=0;
   end
end

% ---------------------------------------------------------------------
function wc=xmipp_DWT_coefs(W,scale,orientation)

Nscale=2^scale;
[k0x,k0y]=quad2ix(scale,0,0,orientation);
[kFx,kFy]=quad2ix(scale,Nscale-1,Nscale-1,orientation);
wc=reshape(W(k0x:kFx,k0y:kFy),Nscale*Nscale,1);
end

% ---------------------------------------------------------------------
function WId=xmipp_bayesian_denoising_LL(WI, scale, orientation, mu, S, N)

WId=WI;
Nscale=2^scale;
[k0x,k0y]=quad2ix(scale,0,0,orientation);
[kFx,kFy]=quad2ix(scale,Nscale-1,Nscale-1,orientation);
y=WI(k0x:kFx,k0y:kFy);
WId(k0x:kFx,k0y:kFy)=(S/(S+N)).*exp(-1/2*((y-mu).^2/(S+N)))./...
     (exp(-1/2*((y-mu).^2/N))+...
      exp(-1/2*((y-mu).^2/(S+N))))...
   .*y;
end
