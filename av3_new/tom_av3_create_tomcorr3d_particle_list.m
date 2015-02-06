function sum =tom_av3_create_tomcorr3d_particle_list(align,quality_threshold, resultFile)
   load(align);
   sum = 1;
   fid =fopen(resultFile,'w');
   for i=1:length(Align)
      % if(Align(i).quality>=quality_threshold)
            splitF = regexp(Align(i).Filename,'/','split');
            resultFilename = splitF(length(splitF));

            fprintf(fid, ['/fs/pool/pool-nickell/ml_processing/data/experimental/MspA_hoffmann/particles/particles_96_ctfcor_PV09/' resultFilename{1} '\n']); %particles
       %     fprintf(fid, [Align(i).Filename '\n']);
            fprintf(fid, ['emfile /fs/pool/pool-nickell/ml_processing/data/experimental/MspA_hoffmann/particles/particles_96_ctfcor_PV09_align_new/wedge_' resultFilename{1} '\n']);
            fprintf(fid, ['allpass\n']);
            fprintf(fid, ['\n']);
            sum = sum +1;
       %end
   end
