function outdata = flood2(I,pos,floodmask)
% Flood fills image I from point (x,y) with color c.
LastFlood = false(size(I));
Flood = LastFlood;
Flood(pos(1),pos(2)) = 1;
Mask = floodmask;
FloodFilter = [0,1,0; 1,1,1; 0,1,0];
se = strel('disk',10);
nhood = getnhood(se);
while any(LastFlood(:) ~= Flood(:))
LastFlood = Flood;
Flood = (imdilate(Flood,nhood)) & Mask;
end
outdata=false(size(I));
outdata(Flood) = 1;

end

