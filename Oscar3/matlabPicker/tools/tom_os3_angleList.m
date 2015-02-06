%%
% generates scanning array according to the angleIncrement and the start
% and stop angle for scanning
function res = tom_os3_angleList(options,dim,angle1)



%% dimension 2
if(dim == 2 )
    
    if(~exist('angle1'))
        angle1.start = options.correlation.angles.start;
        angle1.end = options.correlation.angles.end;
        angle1.inc = options.correlation.angles.increment;
    end;
    
    if(angle1.start < angle1.end)
        
        returnSize = floor((angle1.end(1) - angle1.start(1)) /angle1.inc(1))+1;
        
        res = zeros(returnSize,1);
       
        counter = 1;
        for i=angle1.start:angle1.inc:angle1.end
            
            res(counter) = mod(i,360);
            counter = counter +1;
            
        end;
        
    else
        
        returnSize = floor((360+angle1.end-angle1.start)/angle1.inc)+1;
        %%return size + 1 (angle1.start) counts too
        res = zeros(returnSize,1);
        
        counter = 1;
        for i=angle1.start+1:angle1.inc:angle1.end+360-1
            res(counter) = mod(i,360);
            counter = counter +1;
        end;
    end;
%% dimension 3
else
    
    
    tmp.start   = options.correlation.angles.start(1);
    tmp.inc     = options.correlation.angles.increment(1);
    tmp.end     = options.correlation.angles.end(1);
    
    phi         = tom_os3_angleList(options,2,tmp);

    tmp.start   = options.correlation.angles.start(2);
    tmp.inc     = options.correlation.angles.increment(2);
    tmp.end     = options.correlation.angles.end(2);
    
    psi         = tom_os3_angleList(options,2,tmp);

    tmp.start   = options.correlation.angles.start(3);
    tmp.inc     = options.correlation.angles.increment(3);
    tmp.end     = options.correlation.angles.end(3);
    
    theta       = tom_os3_angleList(options,2,tmp);
    
    %generate array for each phi psi theta combination
    %[phi psi theta]
    counter = 1;
    res = {};
    matrices = {};
    for ph = 1:length(phi)
        for ps = 1:length(psi)
            for th = 1:length(theta)
                
                angles = [phi(ph) psi(ps) theta(th)];
                
                angleMatrix = floor(tom_angles2rotmatrix(angles)*10000);
                
                found = false;
                i = 1;
                while(length(res) >0 && ~found && i <=length(res))
                    found = isequal(angleMatrix,matrices{i});
                    i=i+1;
                end;
                
                if(~found)
                    res{length(res)+1} = angles;
                    matrices{length(matrices)+1} = angleMatrix;                
                end;

            end;
        end;
    end;
    
    
    
    
end;



