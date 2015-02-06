function tom_av2_remove_xray_stack(input,nstd)


if (tom_isemfile(input))
    header=tom_reademheader(input);
    for i=1:size(header.Header,3)
        im=tom_emreadc(input,'subregion',[1 1 i],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
        im.Value=xraycorrect(im.Value,nstd);
        tom_emwritec(input,im.Value,'subregion',[1 1 i],[header.Header.Size(1)-1 header.Header.Size(1)-1 0]);
    end;
    
    
else
    
    st=importdata(input);
    
    for i=1:length(st.textdata)
        try
            im=tom_spiderread(st.textdata{i});
            im.Value=xray_correct(im.Value,nstd);
            tom_spiderwrite(st.textdata{i},im.Value);
            if (mod(i,100)==0)
                disp(st.textdata{i});
            end;
        catch ME
            disp(st.textdata{i});
            disp( ME.stack);
            disp( ME.message);
        end;
    end;
    
end;



function I=xray_correct(J,nstd)

 s=2;
d=std2(J);
meanimg=mean(mean(J));
[x y]=find(J>(meanimg+nstd.*d) | J<(meanimg-nstd.*d)); % factor 10 added by SN
I=J;
if isempty(x)
    return;
else
    for i=1:size(x,1)
        if ((x(i)-s > 0) && (y(i)-s > 0) && (x(i)+s <= size(I,1)) ... % bug fixed FF
                && (y(i)+s <= size(I,2)) )
            I(x(i),y(i))=mask(I,x(i),y(i),d,meanimg,s);
        else
            I(x(i),y(i))=mean(mean(I));
        end;
    end
end
if size(x,1)>0
    disp(['values outside range found: ' num2str(size(x,1))]); %changed by SN
end;




%****** SubFunction mask ***********

function c=mask(A,xcoor,ycoor,dev,meanimg,s)

[xdim ydim]=size(A);
a=A(xcoor-s:xcoor+s,ycoor-s:ycoor+s);
t=find(a>(meanimg+dev) | a<(meanimg-dev));
if isempty(t)
    c=sum(sum(a))/((2*s+1)*(2*s+1));
elseif (size(a,1)*size(a,2) == size(t,1))%bug fixed FF
    c=0;
else
    a(t)=0;
    c=sum(sum(a))/(((2*s+1)*(2*s+1))-size(t,1));
end
