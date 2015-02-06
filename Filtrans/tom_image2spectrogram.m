function out=tom_image2spectrogram(im,window_size,offset)

% test=spectrogram(im(:,1) );
% 
% sz=size(test);
% 
% d=zeros(sz(1),sz(2),size(im,2));
% 
% for i=1:size(im,2)
%     d(:,:,i)=spectrogram(im(:,i)); %,10,10,1000,'yaxis');
% end;

sz=tom_cart2polar(zeros(window_size));
sz=size(sz);

out=zeros(size(im,1),size(im,2),sz(2));



for i=1:size(im,1)
    for ii=1:size(im,2)
        im_s=tom_move(im,[-(i-1) -(ii-1)]);
        figure(h1); tom_imagesc(im_s); drawnow;
        figure(h2); tom_imagesc(im_s(1:window_size(1),1:window_size(2))); drawnow;
        
        ps_tmp=tom_ps(im_s(1:window_size(1),1:window_size(2)) );
        figure(h3); tom_imagesc(ps_tmp); drawnow;
        ps_tmp=sum(tom_cart2polar(ps_tmp),1);
        out(i,ii,1:sz(2))=ps_tmp;
    end;
    disp(num2str(i));
end;

