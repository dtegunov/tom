function im = tom_rectangle(sz,x,y,width,height)

im = zeros(sz(1),sz(2),'single');

im(x:x+width,y:y+height) = 1;

im = im(1:sz(1),1:sz(2));
