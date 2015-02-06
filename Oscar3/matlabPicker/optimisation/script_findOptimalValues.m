function [res files pos mins] = script_findOptimalValues(path,startInterval,endInterval,steps,align2d,files)

if(ischar(align2d))
    load(align2d);
end;

if(strcmp(path(end),'/'))
    slash = '';
else
    slash ='/';
end;
%%
if(~exist('files'))
    
    files = dir([path slash '*mat']);
    files = {files.name};
    
    if(length(files) > 10)
        
        r = ceil(rand(10,1) * length(files));
        f = {};
        for i =1:length(r)
            f{i} = files{r(i)};
        end;
        files = f;
    end;
end;
%%
result = [];
d=1;
e=1;
f=1;

% list = {};
% 
% for a=startInterval(1):steps:endInterval(1)
%    e = 1; 
%     for b=startInterval(2):steps:endInterval(2)
%         f = 1;
% 	    for c=startInterval(3):steps:endInterval(3)
%             list{length(list)+1} = [d e f];
%             f = f+1;
%         end;
%         e=e+1;
%     end;
%     d=d+1;
% end;


for a=startInterval(1):steps:endInterval(1)
   e = 1; 
    for b=startInterval(2):steps:endInterval(2)
        f = 1;
	    for c=startInterval(3):steps:endInterval(3)

%             a = checkKoeff(list,[d e f]);
            
%             if(~isempty(a))
%                 res(d,e,f) = res(a(1),a(2),a(3));
%                 mins(d,e,f) = mi(a(1),a(2),a(3));
%                 
%                 
%                 
%                 continue;
%             end;

            [mea mi] = script_testCombination(files,[a b c],align2d,[path slash]);

            %speichere mittelwert in 3d volumen
            res(d,e,f) = mea;
            mins(d,e,f) = mi;
            disp(['item ' num2str(d) ' ' num2str(e) ' ' num2str(f) ' = ' num2str(a) ' ' num2str(b) ' ' num2str(c) ' :' num2str(mea)]);
            f=f+1;

        end;
        e=e+1;
    end;
    d=d+1;
end;
r.r = res;
r.startI = startInterval;
r.endI = endInterval;
r.step = steps;
pos = script_min(r);



% function checkKoeff(list,b)
% 
% for i=1:length(list)
%     
%     if(sum(b)/3 = )
%     
% end;