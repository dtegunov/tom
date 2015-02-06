function globalOptimum = tom_os3_simulatedAnnealing(startI,endI,step,startPosition,annealing,results,align2d)

currentPosition.position = startPosition;
currentPosition.costs = tom_os3_costFunction(results,currentPosition.position,align2d);

lastBestPosition.position = [];
lastBestPosition.costs = 10000000000;

positionList={};

while(annealing > 0)
        
        
        
        randomPosition = startI + (round(rand(3,1)*100)/100)'.*(endI-startI);

        if(checkPosition(positionList,randomPosition))
            %area has already been visited
            continue;
        end;
        positionList = addToList(positionList,currentPosition.position);
        costs =tom_os3_costFunction(results,randomPosition,align2d);

        if(currentPosition.costs > costs)
            
            lastBestPosition.position = currentPosition.position;
            lastBestPosition.costs = currentPosition.costs;
            
            currentPosition.position = randomPosition;
            currentPosition.costs =  costs;
            disp(['new position ' mat2str(randomPosition) ' costs: ' mat2str(costs)]);
        else
            delta = abs(currentPosition.costs - costs);
            value = exp(-delta/annealing);
%             randomDelta = value * rand;
            randomV = rand;
            disp(['delta ' num2str(delta) ' value ' num2str(value) ' randomV ' num2str(randomV)]);
            if(randomV<=value)
                disp([ 'accept worse value :' num2str(costs) '>' num2str(currentPosition.costs)]);
                %replace if the optimum found so far is worse than the
                %currently found bad value
                if(lastBestPosition.costs >= currentPosition.costs)
                    lastBestPosition.position = currentPosition.position;
                    lastBestPosition.costs = currentPosition.costs;
                end;
                
                currentPosition.position = randomPosition;
                currentPosition.costs =  costs;
            end;
        end;

    annealing = annealing - step;
end;    

costs =tom_os3_costFunction(results,[1 0 0],align2d);
if(currentPosition.costs > costs)
            lastBestPosition.position = currentPosition.position;
            lastBestPosition.costs = currentPosition.costs;
            
            currentPosition.position = randomPosition;
            currentPosition.costs =  costs;
            disp(['new position ' mat2str(randomPosition) ' costs: ' mat2str(costs)]);
end;

costs =tom_os3_costFunction(results,[0 1 0],align2d);
if(currentPosition.costs > costs)
            lastBestPosition.position = currentPosition.position;
            lastBestPosition.costs = currentPosition.costs;
            
            currentPosition.position = randomPosition;
            currentPosition.costs =  costs;
            disp(['new position ' mat2str(randomPosition) ' costs: ' mat2str(costs)]);
end;

costs =tom_os3_costFunction(results,[0 0 1],align2d);
if(currentPosition.costs > costs)
            lastBestPosition.position = currentPosition.position;
            lastBestPosition.costs = currentPosition.costs;
            
            currentPosition.position = randomPosition;
            currentPosition.costs =  costs;
            disp(['new position ' mat2str(randomPosition) ' costs: ' mat2str(costs)]);
end;

% globalOptimum = currentPosition.position;
if(lastBestPosition.costs >= currentPosition.costs)
	globalOptimum.position = currentPosition.position;
	globalOptimum.costs =  currentPosition.costs;
else
	globalOptimum.position = lastBestPosition.position;
        globalOptimum.costs =  lastBestPosition.costs;
end;

%% 
function list = addToList(list,position)
    
    found = numel(find(position >= 1)) > 0;
    i=1;
    while(~found && i<=100)
        
        list{length(list)+1}=position * i;
        i = i+0.1;
    end;
    
    
    
function visited = checkPosition(positionList,randomPosition)

visited = false;
for i =1:length(positionList)
    
    p = positionList{i};
    
    if(numel(find(randomPosition == 0)) == 0)
        q = p / randomPosition;
    else
        q = ones(3,1);
        q(1) = 0;
    end;
    
    if(sqrt((p - randomPosition).^2) < 0.1 )
        visited = true;
        disp(['pruning ' mat2str(randomPosition)]);
        return;
    end;
end;
