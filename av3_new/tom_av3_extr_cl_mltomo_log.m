function [st Align]=tom_av3_extr_cl_mltomo_log(f_log_ml,f_partlist_org,f_new_part_list,class_vect,supp_thr,peak_thr,iteration_number)
%  tom_av3_extr_cl_mltomo_log parse and filter mltomo.log file
%  
%     [st Align]=tom_av3_extr_cl_mltomo_log(f_log_ml,f_partlist_org,f_new_part_list,class_vect,supp_thr,peak_thr,iteration_number)
%  
%  PARAMETERS
%  
%    INPUT
%     log_ml_f            name of the mltomo log file
%     f_partlist_org      name of the org input particle list
%     f_new_part_list     filename for new particle_list
%     class_vect          vector containing the classes which should be
%                         extracted (starts counting at 0 ...chek log file)
%     supp_thr            (0)   threshold for centroid (template) probability                               
%                          0%    no particles are sorted out  
%                          100%  only fully assigned parts are used! 
%    
%    
%     peak_thr            (0)    threshold for translation probability   
%                          0%    no particles are sorted out  
%                          100%  only fully assigned parts are used!    
%
%    iteration_number    (last iter) iteration number used for extration
%                          default is the last iteration  
%
%    
%    OUTPUT
%     st                parsed log filst struct
%     Align             Align struct for tom/av3_new functions           
%  
%  EXAMPLE
%    
%   
%   %example using class number 0 and 1 ,supp_thr of 98% , peak_thr 0.9 (only almost conv partrs are used)
%   %iteration 11 was used for part extraction
%   [st,Align]=tom_av3_extr_cl_mltomo_log('mltomo.log','part_list_absFB.txt','sorted01.txt',[0 1],98,0.9,11);
% 
%   %plot conv of particles
%   figure; hist(st.supp(:),20);
%   
%   %genete a average of the Align
%   avg=tom_av3_average3(Align(end,:),'inverse',99,[0 1]);
%   figure; tom_dspcub(avg);
%
%   %sort particle with 
%   save('my_alg.mat','Align'); 
%   tom_av3_stackbrowser   %browse for my_alg.mat in tom_av3_stackbrowse
%   
%
%  NOTE:
% 
%  Use abs filename in f_partlist_org (...further processing)
%  
%  Align(run,i).CCC is used 2 store supp (centroid probability)
%  
%  RUN HAS TO BE FINISHED TO USE THE DEFAULT (LAST ITERATION)
%  if the run was user interrupted use iteration_number 
%  
%  the mltomo.log contains the entryCount in the particle/template-list not the particle/template number !!! 
%
%  To have all information about the particle use the format of tom_particles   
%
%  emheader of every image should contain: 
%
%  Header.comment: Position.x Position.y Position.z offset binning TomogramFilename
%  Header.Tiltangle: minimum tilt angle        (overloaded for storage)
%  Header.Header.Tiltaxis: maximum tilt angle  (overloaded for storage)
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_av3_average3,tom_av3_stackbrowser,tom_particles
%  
%     created by fb okt.2011
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (nargin < 5)
    supp_thr=0;
end;

if (nargin < 6)
    peak_thr=0;
end;

if (nargin < 7)
   call=['cat ' f_log_ml ' | grep -o ''' 'itr=....''' ' | tail -1 '];
   [a out]=unix(call);
   iteration_number=str2double(strtok(out,'iter='));
end;


show_param(f_log_ml,f_partlist_org,f_new_part_list,class_vect,supp_thr,peak_thr,iteration_number);

fid=fopen(f_partlist_org,'r');
i=1;
num_of_lines=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break; 
    end;
    num_of_lines=num_of_lines+1; 
    part_list.name{i}=tline;
    tline = fgetl(fid);
    if ~ischar(tline)
        break; 
    end;
    num_of_lines=num_of_lines+1;
    part_list.wedge{i}=tline;
    part_list.entry(i)=i-1;
    i=i+1;
end;
    



for i=1:iteration_number
    
    disp(['Using Iteration number: ' num2str(iteration_number)]);
    
    if (iteration_number<10)
        it_blank='  ';
    end;
    
    if (iteration_number>=10 && iteration_number<100)
        it_blank=' ';
    end;
    
    if (iteration_number>=100 )
        it_blank='';
    end;
    
    
    %get line number
    call=['cat ' f_log_ml ' | grep -n ''' '; itr=' it_blank num2str(iteration_number) ' : sigma''' ];
    [a out]=unix(call);
    lnr=str2double(strtok(out,':'));
    
    if isempty(out)
        iteration_number=iteration_number-1;
        disp(['No entries in logfile for: ' num2str(iteration_number+1) ]);
        continue;
    end;
    
    %get entry number
    call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  -o "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = ":" } ; {print $3} ''' ' | sed ''' 's/p//''' ];
    [a out]=unix(call);
    p_entry_nr=sscanf(out,'%d \n');
    
    if (isempty(out)==0)
        break;
    end;
    iteration_number=iteration_number-1;
    disp(['No entries in logfile for: ' num2str(iteration_number+1) ]);
    
end;

disp(' ');

%get template nr
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  -o "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = ":" } ; {print $2} ''' ' | sed ''' 's/t//''' ];
[a out]=unix(call);
p_templ=sscanf(out,'%d \n');

disp(['Run was using ' num2str(max(p_templ)+1) ' templates ']);

%get supp
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  -o "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = ":" } ; {print $4} ''' ' | awk ''' '{ print $2}'  ''''];
[a out]=unix(call);
supp=sscanf(out,'%f \n');

%get reference number
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = "=" } ; {print $3} '''  '| awk ' '''  BEGIN { FS = "[" } ; {print $1} '''];
[a out]=unix(call);
ref_num=sscanf(out,'%f \n');


%get angles
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = "[" } ; {print $2} ''' ];
[a out]=unix(call);
tmp=sscanf(out,'%f , %f , %f ]\n');
angles=reshape(tmp,[3 ,(length(tmp)/3)]);

%get shifts
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = "]" } ; {print $2} ''' ];
[a out]=unix(call);
tmp=strrep(strrep(out,'[',''),' ','');
tmp=sscanf(tmp,'%f,%f,%f\n');
shifts=reshape(tmp,[3 ,(length(tmp)/3)]);

%get peak val
call=['cat ' f_log_ml ' | sed ''' '1,' num2str(lnr) 'd' ''' ' ' | grep  "  ;' it_blank '.*:.*t.*p.*:.*%" | awk ' '''  BEGIN { FS = "=" } ; {print $NF} ''' ];
[a out]=unix(call);
peak_val=sscanf(out,'%f \n');




u_entry=unique(p_entry_nr);
for i=1:length(u_entry)
    st.entry_num(i)=u_entry(i);
    tmp_idx=find(p_entry_nr==u_entry(i)) ;
    [supp_val pos]=max(supp(tmp_idx)); %-1 class count starts at 0
    abs_idx=tmp_idx(pos);
    st.class(i)=pos-1;
    st.supp(i)=supp_val; 
    st.ref_num(i)=ref_num(abs_idx);
    st.angles(:,i)=angles(:,abs_idx)';
    st.shifts(:,i)=shifts(:,abs_idx)';
    st.peak_val(:,i)=peak_val(abs_idx);
    idx_part=find([part_list.entry==u_entry(i)]);
    if (isempty(idx_part))
        disp(['Cannot find entry: ' num2str(u_entry(i)) ' in ' f_partlist_org]);
        disp(['Length of: ' f_partlist_org ': ' num2str(num_of_lines) ' requsted ' num2str(((max(u_entry)+1)*2)) ]);
        error(['logfile and particle-list do not fit!']);
    end;
    st.part_name{i}=part_list.name{idx_part};
    st.part_wedge{i}=part_list.wedge{idx_part};
end;



fid=fopen(f_new_part_list,'wt');
zz=0;
used_idx=zeros(length(u_entry),1);
for i=1:length(u_entry)
  if ((isempty(find(ismember(class_vect,st.class(i))))==0)  && st.supp(i) >= supp_thr &&  st.peak_val(i) >= peak_thr )
    fprintf(fid,'%s\n',st.part_name{i});
    fprintf(fid,'%s\n',st.part_wedge{i});
    zz=zz+1;
    used_idx(zz)=i;
  end;
end;
tmp=used_idx(1:zz);
used_idx=tmp;
fclose(fid);



disp(['Length of: ' f_partlist_org  ': ' num2str(num_of_lines) ]);
disp(['Length of: ' f_new_part_list  ': ' num2str(zz.*2) ]);
disp(['Used parts: ' num2str(zz) ]);
disp(['Removed parts: ' num2str(length(u_entry) -zz) ]);


if (nargout > 1)
    
    disp('Building Aling struct: ');
    
    run=1;
    
    if (tom_isemfile(st.part_name{i})~=1)
        disp('Warning particles not readable');
        dummy_head=tom_emheader(st.part_name{i});
    end;
    
    for i=1:length(used_idx)
        
        
        f_name=st.part_name{used_idx(i)};
        
        try
            header = tom_reademheader(f_name);
        catch Me
            header.Header = dummy_head.Header;
        end;
            
        
        %Evaluate comment, format: Position.x Position.y Position.z offset binning TomogramFilename
        comment = char(header.Header.Comment');
        token = zeros(1,5);
        try
            if size(strtrim(comment),1) > 1
                for j = 1:5
                    [tok, comment] = strtok(comment);
                    token(j) = str2num(tok);
                end
            end
        end;
        
        
        Align(run,i).Filename = f_name;
        Align(run,i).Tomogram.Filename = strtrim(comment);
        Align(run,i).Tomogram.Header = header.Header;
        Align(run,i).Tomogram.Position.X = token(1); %Position of particle in Tomogram (values are unbinned)
        Align(run,i).Tomogram.Position.Y = token(2);
        Align(run,i).Tomogram.Position.Z = token(3);
        Align(run,i).Tomogram.Regfile = '';
        Align(run,i).Tomogram.Offset = token(4);     %Offset from Tomogram
        Align(run,i).Tomogram.Binning = token(5);    %Binning of Tomogram
        Align(run,i).Tomogram.AngleMin = header.Header.Tiltangle;
        Align(run,i).Tomogram.AngleMax = header.Header.Tiltaxis;
        Align(run,i).Shift.X = st.shifts(1,used_idx(i)); %Shift of particle, will be filled by tom_av3_extract_anglesshifts
        Align(run,i).Shift.Y = st.shifts(2,used_idx(i));
        Align(run,i).Shift.Z = st.shifts(3,used_idx(i));
        Align(run,i).Angle.Phi = st.angles(1,used_idx(i)); %Rotational angles of particle, will be filled by tom_av3_extract_anglesshifts
        Align(run,i).Angle.Psi = st.angles(2,used_idx(i));
        Align(run,i).Angle.Theta = st.angles(3,used_idx(i));
        Align(run,i).Angle.Rotmatrix = []; %Rotation matrix filled up with function tom_align_sum, not needed otherwise
        Align(run,i).CCC = st.supp(used_idx(i)); % cross correlation coefficient of particle, will be filled by tom_av3_extract_anglesshifts
        Align(run,i).PeakVal=st.peak_val(used_idx(i)); 
        Align(run,i).Class = st.class(used_idx(i));
        Align(run,i).ProjectionClass = st.ref_num(used_idx(i));
        Align(run,i).NormFlag = 0; %is particle phase normalized?
        Align(run,i).Filter = [0 0]; %is particle filtered with bandpass?
        
        
    end;
    
end;

function show_param(f_log_ml,f_partlist_org,f_new_part_list,class_vect,supp_thr,peak_thr,iteration_number)


disp(' ');
disp('======================================>');
disp(['f_log_ml: ' f_log_ml]);
disp(['f_partlist_org: ' f_partlist_org]);
disp(['f_new_part_list: ' f_new_part_list]);
disp(['class_vect: ' num2str(class_vect)]);
disp(['supp_thr: ' num2str(supp_thr)]);
disp(['peak_thr: ' num2str(peak_thr)]);
disp(['iteration_number: ' num2str(iteration_number)]);
disp('<======================================');
disp(' ');


