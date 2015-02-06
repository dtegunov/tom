function statistic=tom_corr_subregions_clean_statistic(statistic,nr_std)


%TOM_CORR_SUBREGIONS_CLEAN_STATISTIC refines the statistics generated by
%tom_corr_subregions by eliminating outliers, n-times the standard
%deviation off from the mean value.
%
%   statistic=tom_corr_subregions_clean_statistic(statistic,nr_std);
%
%
%PARAMETERS
%
%  INPUT
%   statistic           shift statistics generated by tom_corr_subregions.
%   nr_std              n*times the standard deviation off from the
%                       mean-value of shifts.
%                     
%  
%  OUTPUT
%   statistic           cleaned shift statistics.
%
%EXAMPLE
%
%   statistic=tom_corr_subregions_clean_statistic(statistic,3);
%
%REFERENCES
%
%SEE ALSO
%   TOM_CORR_SUBREGIONS, TOM_CORR
%
%   created by SN 30/09/08
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

    
hit=1;
nohit=1;
hits=statistic.shift_hit;
if (isfield(statistic,'shift_nohit'))
    nohits=statistic.shift_nohit;
else
    nohits=[];
end;

if (isfield(statistic,'shift_nohit')==0)
    statistic.pos_nohit=[];
end;


allhits=[hits; nohits];
statistic.shift_nohit=[];
statistic.shift_hit=[];
pos_hit=statistic.pos_hit;
statistic.pos_hit=[];
pos_nohit=statistic.pos_nohit;
statistic.pos_nohit=[];
allpos=[pos_hit; pos_nohit];

for i=1:size(allhits,1)
    devx_min=statistic.shift_hit_mean(1)-statistic.shift_hit_std(1).*nr_std;
    devx_max=statistic.shift_hit_mean(1)+statistic.shift_hit_std(1).*nr_std;

    devy_min=statistic.shift_hit_mean(2)-statistic.shift_hit_std(2).*nr_std;
    devy_max=statistic.shift_hit_mean(2)+statistic.shift_hit_std(2).*nr_std;

    
    if allhits(i,1)<devx_min ||  allhits(i,1)>devx_max ...
            || allhits(i,2)<devy_min || allhits(i,2)>devy_max 
                statistic.shift_nohit(nohit,:)=allhits(i,:);
                statistic.pos_nohit(nohit,:)=allpos(i,:);
                nohit=nohit+1;
    else
                statistic.shift_hit(hit,:)=allhits(i,:);
                statistic.pos_hit(hit,:)=allpos(i,:);
                hit=hit+1;
    end;
end;
statistic.hit=hit-1;
statistic.nohit=nohit-1;
statistic.shift_hit_mean(1)=mean(statistic.shift_hit(:,1));
statistic.shift_hit_mean(2)=mean(statistic.shift_hit(:,2));
statistic.shift_hit_std(1)=std(statistic.shift_hit(:,1));
statistic.shift_hit_std(2)=std(statistic.shift_hit(:,2));
