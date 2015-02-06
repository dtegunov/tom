function imout = minneighbors(I, extent)

imout = I;

for y=1:size(I,2)
    for x=1:size(I,1)
        if I(x, y)
            minx = max(1,x-extent);
            maxx = min(size(I,1),x+extent);
            miny = max(1,y-extent);
            maxy = min(size(I,2),y+extent);
            
            sizex = maxx - minx + 1;
            sizey = maxy - miny + 1;
            
            if sum(sum(I(minx:maxx,miny:maxy))) < min(sizex, sizey)
                imout(x,y) = 0;
            end;
        end;
    end;
end;

end

