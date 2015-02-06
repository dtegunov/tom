function align2d = tom_os3_pickList2Align2d(pickList,binning)


if(nargin == 1)
    binning = 1;
end;

align2d.dataset='';
align2d.filename ='';
align2d.position.x =0;
align2d.position.y =0;
align2d.class='default';
align2d.radius =0;
align2d.color =[0 0 0];
align2d.shift.x =0;
align2d.shift.y =0; 
align2d.angle=0;
align2d.isaligned=0;
align2d.ccc=0;
align2d.quality=0;
align2d.normed='none';



for i=1:length(pickList)
    pick = pickList{i};
    
    align2d(i).filename = pick.filename;
    align2d(i).position.x = pick.coordinates(1);
    align2d(i).position.y = pick.coordinates(2);
    align2d(i).radius = pick.templateSize(1)*2^(binning-1);

    align2d(i).color = [0 0 1];
    align2d(i).shift.x =0;
    align2d(i).shift.y =0; 
    if(pick.value > 1)
        pick.value = 0;
    end;
    align2d(i).ccc =  pick.value;
    align2d(i).isaligned=0;
    align2d(i).quality=0;
    align2d(i).normed='none';
    align2d(i).class='automatic';
    align2d(i).angle= pick.angle;
end;
