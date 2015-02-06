function Align = tom_av3_parseAlign_tomcorr3d(alignmentFile, particleFile, fileprefix, resultDir)
   %alignmentFile -> File with the alignment parameters e.g. hoffmann_mpsa0000_par_alignfile.txt
   %particleFile -> File with the particle file names e.g. hoffmann_mpsa0000_par_particles_list
   %fileprefix -> should be '' in case of absolute path and the path to the
   %in case of relative paths in the particle File.
   
    binning = 1;
    i = 1;
    fid = fopen(particleFile);
    line=1;
    part_no =0;
    while 1
        tline = fgetl(fid);
        if(tline==-1)
            break;
        end
        
        if(mod(line,4)==1)
            part_no = part_no+1;
            if regexp(tline,'\.\/')
                filenames{part_no} = [fileprefix tline(3:end)];
            else
                filenames{part_no} = tline;
            end
        end
        if(mod(line,4)==2)
            % extracting the wedgefilename if any is there
            wedge_tmp = regexp(tline,' *','split');
            wedge_tmp = wedge_tmp(length(wedge_tmp));
            wedge_tmp = regexprep(wedge_tmp{1}, '"', '');
            regexp(wedge_tmp,'\.\/','split');
            if regexp(wedge_tmp,'\.\/')
                wedge_info{part_no} = [fileprefix wedge_tmp(3:end)];
            else
                wedge_info{part_no} = wedge_tmp;
            end
        end
        
        line = line +1;
    end    
    fid = fopen(alignmentFile);
    while 1
        tline = fgetl(fid);
        
        if(tline==-1)
            break;
        end
        if(regexp(tline,'^#.'))
            % Comment
        elseif (regexp(tline,'^binning'))
            r = regexp(tline,'=','split');
            binning = str2double(r{2});
        elseif (tline) % prune empty lines
            run = 1;
            values = regexp(tline,' *','split');
            %Align(run,i).Filename = filenames{i};
            Align(run,i).Filename = filenames{str2num(values{3})+1}; %change fb to cope with split particle-list !! (fsc ole ole)
            Align(run,i).WedgeInfo = wedge_info{i};
            Align(run,i).Tomogram.Filename = '';
            Align(run,i).Tomogram.Header = '';
            Align(run,i).Tomogram.Position.X = 0; %Position of particle in Tomogram (values are unbinned)
            Align(run,i).Tomogram.Position.Y = 0;
            Align(run,i).Tomogram.Position.Z = 0;
            Align(run,i).Tomogram.Regfile = '';
            Align(run,i).Tomogram.Offset = '';        %Offset from Tomogram
            Align(run,i).Tomogram.Binning = binning;  %Binning of Tomogram
            Align(run,i).Tomogram.AngleMin = ''; 
            Align(run,i).Tomogram.AngleMax = '';
            Align(run,i).Shift.X = str2double(values{7}); %Shift of particle, will be filled by tom_av3_extract_anglesshifts
            Align(run,i).Shift.Y = str2double(values{8});
            Align(run,i).Shift.Z = str2double(values{9});
            Align(run,i).Angle.Phi = str2double(values{4}); %Rotational angles of particle, will be filled by tom_av3_extract_anglesshifts
            Align(run,i).Angle.Psi = str2double(values{5});
            Align(run,i).Angle.Theta = str2double(values{6});
            Align(run,i).Angle.Rotmatrix = []; %Rotation matrix filled up with function tom_align_sum, not needed otherwise
            Align(run,i).CCC = str2double(values{10}); % cross correlation coefficient of particle, will be filled by tom_av3_extract_anglesshifts
            Align(run,i).Class = 1;
            Align(run,i).ProjectionClass = 0;
            Align(run,i).NormFlag = 0; %is particle phase normalized?
            Align(run,i).Filter = [0 0]; %is particle filtered with bandpass?
            i=i+1;
        end 
    end
    
    save([resultDir 'Align.mat'], 'Align');
end