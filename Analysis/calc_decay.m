function [decay decay_img]=calc_decay(ps,mask_in,mask_out,sampling)

% [decay decay_img]=calc_decay(psb,8,120,32);

pol=tom_cart2polar(ps);
polc=pol(mask_in:mask_out,:);
idx=1;
for i=1:sampling./2:size(polc,2)-(sampling-1)
    s=sum(polc(:,i+(sampling-1)),2);

    % needs Curve Fitting TB, out, SN, 4.3.2010
    %    [curve, goodness] = fit( double([1:size(polc,1)]'),double(s), 'Poly4' );

    %    decay(idx,:)=feval(curve,[1:size(polc,1)]');
    [P,S,MU]=polyfit( double([1:size(polc,1)]'),double(s), 4 );
    decay(idx,:) = polyval(P,double([1:size(polc,1)]'),[],MU);
    idx=idx+1;

end;
idx=1;
decay_img=zeros(size(pol));
for i=1:sampling./2:size(pol,2)-(sampling-1)
    decay_img(mask_in:mask_out,i:i+sampling-1)=repmat(decay(idx,:)',[1 sampling]);    
    idx=idx+1;
end;
decay_img=tom_polar2cart(decay_img);