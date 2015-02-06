function sysinfo=tom_systeminfo(varargin)
%TOM_SYSTEMINFO creates ...
%
%   sysinfo=tom_systeminfo(varargin)
%
%   out=tom_systeminfo
%   out=tom_systeminfo('name1','name2','name3');
%
%PARAMETERS
%
%  INPUT
%   Give information about a computer. The result is a structure containing:
%   out.name           :name of the computer    
%   out.Accessible     :the computer is accesible = 1, else = 0
%   out.Type           :type of the computer. Should be 'Linux'
%   out.NrOfProcessors :Number of processors of the computer
%   out.NrOfUser       :Number of Users connected to the computer
%   out.Memory         :[Memory Total, Memory used, Memory free]
%   out.LoadAverage    :load averages for the past 1, 5, and 15 minutes
%   out.Uptime         :The current time, how long the system has been running
%
%   out=tom_systeminfo  get information of the default list of Baumeister computer
%    default list:
%    'cairo','dublin','haifa','jakarta','kyoto','lima','montreal','nairobi',
%    'tucson','xian'
%
%   out=tom_systeminfo('name1','name2','name3')  get information of the list of
%    computer: 'name1','name2','name3'
%
%   Note: This function is using ssh. If a password is required, this function
%    will not work (see manual of SSH)
%
%EXAMPLE
%   a=tom_systeminfo  get information of HP workstation in Baumeister group
%
%   a=tom_systeminfo('lima','cluster01','sally2')  get information of computer
%                                                  lima, cluster01 and sally2
%
%REFERENCES
%
%SEE ALSO
%   TOM_EMREAD, TOM_LAMBOOT
%
%   created by WDN 02/01/06
%   updated by WDN 20/02/06
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


Ordi=computer;
switch Ordi
    case {'GLNX86','GLNXA64'}
        sysinfo=struct('Name','','Accessible','','Type','','NrOfProcessors','',...
                       'NrOfUser','','Memory','','LoadAverage','','Uptime','');
        if isempty(varargin)
            %default list
            a=['cairo    ';'dublin   ';'haifa    ';'jakarta  ';'kyoto    ';,...
               'lima     ';'montreal ';'nairobi  ';'tucson   ';'xian     '];
            BaumeisterComputer=cellstr(a);
        else %use the list given by the user
            if iscellstr(varargin)                
                BaumeisterComputer=varargin';
            else
                BaumeisterComputer=varargin{1};
            end
        end
        disp('Wait a moment. Processing ....');
        for i=1:size(BaumeisterComputer,1)
            sysinfo(i).Name=BaumeisterComputer{i};
            [status,message]=unix(['ping -w 1 -c 1 ' BaumeisterComputer{i}]);
            if status~=1
                [sysinfo(i).Accessible,...
                    sysinfo(i).Type]=Nsys(BaumeisterComputer{i});
                if sysinfo(i).Accessible==1
                    sysinfo(i).NrOfProcessors=Nproc(BaumeisterComputer{i});
                    sysinfo(i).NrOfUser=Nuser(BaumeisterComputer{i});                    
                    sysinfo(i).Memory=Nmem(BaumeisterComputer{i});                    
                    [sysinfo(i).LoadAverage,...
                        sysinfo(i).Uptime]=Nave(BaumeisterComputer{i});
                end
            else
                sysinfo(i).Accessible=0;
               
            end
        end
    otherwise
        disp('This function works only under a linux system');
end


%----------------------------------
%--------  EXTRA FUNCTION  --------
%----------------------------------
function [acc nam]=Nsys(HostType)
[status,message]=unix(['ssh ' HostType ' uname -a']);%unix(['ssh ' HostType ' uname -a'],'-echo');
if status==0
    acc=1;
    if findstr('assword:',message)
        message(1:findstr('assword:',message)+8)='';     
        nam=fliplr(deblank(fliplr(deblank(message(1:findstr(message,HostType)-1)))));
    else
        nam=deblank(message(1:findstr(message,HostType)-1));
    end
else
    acc=0;
    nam='';
end
      
%----------------------------------
function nbup=Nproc(HostComp)
[status,message]=unix(['ssh ' HostComp ' cat /proc/cpuinfo |grep processor']);
nbup=str2num(message(size(message,2)-2:size(message,2)))+1;

%----------------------------------
function nbus=Nuser(HostUser)
[status,message]=unix(['ssh ' HostUser ' who -q |grep users']);
nbus=str2num(message(findstr(message,'=')+1:size(message,2)));

%----------------------------------
function mem=Nmem(HostMem)
%return Memory Total, Memory Used, Memory free, in Kb
[status,message]=unix(['ssh ' HostMem ' cat /proc/meminfo']);
message1=message(findstr(message,'MemTotal:')+9:findstr(message,'MemFree:')-1);
message1=deblank(fliplr(deblank(fliplr(message1))));
message2=message(findstr(message,'MemFree:')+8:findstr(message,'Buffers:')-1);
message2=deblank(fliplr(deblank(fliplr(message2))));
mem(1)=str2num(message1(1:findstr(message1,'kB')-1));
mem(3)=str2num(message2(1:findstr(message2,'kB')-1));
mem(2)=mem(1)-mem(3);

%----------------------------------
function [ave upt]=Nave(HostAve)
%return Average after 1mn, 5mn, 15mn
[status,message]=unix(['ssh ' HostAve ' uptime']);
message=(message(findstr(message,'up')-10:size(message,2)));
nbkoma=findstr(message,',');
sk=size(nbkoma,2);
upt=deblank(fliplr(deblank(fliplr(message(1:nbkoma(1)-1)))));
ave(1)=str2num(message(nbkoma(sk-1)-4:nbkoma(sk-1)-1));
ave(2)=str2num(message(nbkoma(sk)-4:nbkoma(sk)-1));
ave(3)=str2num(message(nbkoma(sk)+1:nbkoma(sk)+5));
