%%
function [img1 angleMatrix psr autoc] = tom_os3_bestPeak(img1,img2, angleMatrix, angl,psr,psrMap,autoc,autocorrMap,dim,options)
%replace peak maximum if necessary and its corresponding angle value
%tom_os3_bestPeak(maxi,peaks,ang,i,psr,psrMap,2)


%%
%2Do
%ersetze schleifen mit schneller vektor maximum operation aus matlab
%bibliothek

%%
    if(options.analysis.ccc == 1)
        optimum = @ccc2d;
    end

    if(options.analysis.psr == 1)
        optimum = @psr2d;
    end
    
    if(options.analysis.autocorr == 1)
        optimum = @aut2d;
    end
    
    if(options.analysis.ccc == 1 && options.analysis.psr == 1)
        optimum = @cccpsr2d;
    end;
    
    if(options.analysis.ccc == 1 && options.analysis.autocorr == 1)
        optimum = @cccaut2d;
    end;
    
    if(options.analysis.psr == 1 && options.analysis.autocorr == 1)
        optimum = @psraut2d;
    end;    
    
    if(options.analysis.ccc == 1 && options.analysis.psr == 1 && options.analysis.autocorr == 1)
        optimum = @cccpsraut2d;
    end;  
%%

mask = find(optimum(img1,img2,psr,psrMap,autoc,autocorrMap));

img1(mask) = img2(mask);
angleMatrix(mask) = angl(mask);
psr(mask) = psrMap(mask);
autoc(mask) = autocorrMap(mask);

%if(dim == 2)

    
    
%     for x = 1:(size(img1,1))
%         for y = 1:(size(img1,2))
% %if each value is smaller than the new, replace the values
%             if(optimum(img1(x,y) ,img2(x,y)  ,psr(x,y) ,psrMap(x,y) ,autoc(x,y) ,autocorrMap(x,y)))
% %           if(img1(x,y) < img2(x,y)  && psr(x,y) < psrMap(x,y) && autoc(x,y) < autocorrMap(x,y))
% %             if(img1(x,y) < img2(x,y)  && psr(x,y) < psrMap(x,y))
% %             if(img1(x,y) < img2(x,y))
% %             if(psr(x,y) < psrMap(x,y))
% %             if(autoc(x,y) < autocorrMap(x,y))
%                 img1(x,y) = img2(x,y);
%                 angleMatrix(x,y) = angl(x,y);
%                 psr(x,y) = psrMap(x,y);
%                 autoc(x,y) = autocorrMap(x,y);
%             end;
%         end;
%     end;

%%    
% else
%     for x = 1:(size(img1,1))
%         for y = 1:(size(img1,2))
%             for z=1:size(img1,3)
% %                 if((img1(x,y,z)) < (img2(x,y,z)))
% %                 if(psr(x,y,z)  < psrMap(x,y,z))  
% %                 if(autoc(x,y,z) < autocorrMap(x,y,z))
% %                 if(img1(x,y,z) < img2(x,y,z) && psr(x,y,z)  < psrMap(x,y,z))
% %                 if(img1(x,y,z) < img2(x,y,z) && autoc(x,y,z) < autocorrMap(x,y,z))    
%                 if(img1(x,y,z) < img2(x,y,z) && autoc(x,y,z) < autocorrMap(x,y,z) && psr(x,y,z)  < psrMap(x,y,z))  
%                     img1(x,y,z) = img2(x,y,z);
%                     angleMatrix(x,y,z) = angl(x,y,z);
%                     psr(x,y,z) = psrMap(x,y,z);
%                     autoc(x,y,z) = autocorrMap(x,y,z);
%                 end;
%             end;
%         end;
%     end;
%     
% end;

end
%%
function r =ccc2d(a,b,c,d,e,f)
    
    r = a<b;
    
end

%%
function r =psr2d(a,b,c,d,e,f)
    
    r = c<d;
    
end

%%
function r =aut2d(a,b,c,d,e,f)
    
    r = e<f;
    
end

%%
function r =cccpsr2d(a,b,c,d,e,f)
    
    r = a<b & c<d;
    
end

%%
function r =cccaut2d(a,b,c,d,e,f)
    
    r = a<b & e<f;
    
end

%%
function r =psraut2d(a,b,c,d,e,f)
    
    r = c<d & e<f;
    
end

%%
function r =cccpsraut2d(a,b,c,d,e,f)
    
    r = a<b & c<d & e<f;
    
end