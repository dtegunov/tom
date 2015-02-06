function tom_av3_align_particle_stack(alignFile,resultDir,mask,wedgeInFile)
load(alignFile);

particle = tom_emread( Align(1).Filename);particle = particle.Value;
m = size(particle,1);
avg = zeros(m,m,m);
compound_wedge = zeros(m,m,m);
wedge = tom_wedge(zeros(m,m,m), 30);


for i = 1:length(Align) 
    filename  = Align(i).Filename;
    shift     = Align(i).Shift;
    rot       = Align(i).Angle;
    Align(i).Shift.X = 0;
    Align(i).Shift.Y = 0;
    Align(i).Shift.Z = 0;
    Align(i).Angle.Phi = 0;
    Align(i).Angle.Psi = 0;
    Align(i).Angle.Theta = 0;
    % load particle
    particle  = tom_emread(filename);particle = particle.Value;
    
    % save aligned particle
    rotPart = double(tom_rotate(tom_shift(particle,[shift.X shift.Y shift.Z]),[rot.Phi rot.Psi, rot.Theta]));
    splitF = regexp(filename,'/','split');
    resultFilename = splitF(length(splitF));
    tom_emwrite([resultDir resultFilename{1}],rotPart.*mask );
    Align(i).Filename = [resultDir resultFilename{1}];

    % load wedge if a wedge file is given (import if wedge is not aligned the
    % z axis)
    if(wedgeInFile==1)
        wedge = tom_emread([ Align(i).WedgeInfo]);
        wedge = wedge.Value;
    end
    
    % save aligned wedge
    rotwedge   = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[rot.Phi rot.Psi, rot.Theta])),0.5,1,'z'),0,0.5);
    tom_emwrite([resultDir 'wedge_' resultFilename{1} ],rotwedge );
        
    avg = avg + tom_norm(rotPart,'mean0+1std');
    compound_wedge = compound_wedge + rotwedge;
    
    disp(['Particle no ' Align(i).Filename ' added to average']);
end;
tom_emwrite([resultDir 'avg.em'],avg );
tom_emwrite([resultDir 'compound_wedge.em'], compound_wedge );

avg_corr = tom_av3_correct_average(avg, compound_wedge);
tom_emwrite([resultDir 'avg_corr.em'],avg_corr );
disp(['Averaging finished - ' num2str(i) ' particles averaged ... '  ]);

save('Align.mat', 'Align');
