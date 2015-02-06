function trlen = lentrk(tr);


s = size(tr);

ndat = s(2)

[b,u] = unique(tr(:,ndat));
ntracks = length(u)
u = [0;u];

for i=2:ntracks
    res(i-1,1) = tr(u(i),ndat-1) - tr(u(i-1)+1, ndat-1);
    res(i-1,2) = tr(u(i),ndat);
end

trlen=res;