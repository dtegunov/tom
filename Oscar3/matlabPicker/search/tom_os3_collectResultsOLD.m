function [peakVolume anglesVolume psrVolume autocVolume] = tom_os3_collectResults(results,volumeSize)
%MY_COLLECTRESULTS creates ...
%
%   tom_os3_collectResults()
%
%PARAMETERS
%
%  INPUT
%   filelabel           ...
%  
%  OUTPUT
%   data		...
%
%EXAMPLE
%   tom_os3_collectResults(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

%% init variables
    
    peakVolume     = zeros(volumeSize);
    anglesVolume   = peakVolume;
    psrVolume      = peakVolume;
    autocVolume    = peakVolume;
    numberResults  = length(results.peaks);
    done           = zeros(numberResults,1);
%%
    for resultIterator = 1:numberResults

        %avoid working on the same result set more than once
        if(~done(resultIterator))
            
            job = results.jobs{resultIterator};
            
            peakSubVolume  = results.peaks{resultIterator};
            angleSubVolume = results.angles{resultIterator};
            psrSubVolume   = results.psr{resultIterator};
            autocSubVolume = results.autocorr{resultIterator};
            
            %find another result set which with the same coordinates
            findIterator = 1;
            while(findIterator <= numberResults)
            
                if(~done(findIterator) && findIterator ~= resultIterator)
                    j = results.jobs{findIterator};
                    if(isequal(job.coordinates,j.coordinates))
                        
                        %apply the maximum operation to both volumes
                        [peakSubVolume angleSubVolume psrSubVolume autocSubVolume] = tom_os3_bestPeak(peakSubVolume,results.peaks{findIterator}, ...
                                                                                                 angleSubVolume,results.angles{findIterator}, ...
                                                                                                 psrSubVolume,results.psr{findIterator}, ...
                                                                                                 autocSubVolume,results.autocorr{findIterator}, ...
                                                                                                 job.dimension);
                        
                        done(findIterator) = 1;
                    end;
                end;
                findIterator = findIterator + 1;
            end;
            
            %use a b c d as placeholders for the long expressions
            
            a = peakVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6));
            b = anglesVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6));
            c = psrVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6));
            d = autocVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6));
            
            [a b c d] = tom_os3_bestPeak(a , peakSubVolume, b , angleSubVolume, c , psrSubVolume , d , autocSubVolume,job.dimension);
            
            peakVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6)) = a;
            anglesVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6)) = b;
            psrVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6)) = c;
            autocVolume(job.coordinates(1):job.coordinates(2),job.coordinates(3):job.coordinates(4),job.coordinates(5):job.coordinates(6)) = d;
            
            
            done(resultIterator) = 1;
        end;
    end;


    res.peaks   = peakVolume;
    res.angles  = anglesVolume;
    res.psr     = psrVolume;
    res.autoc   = autocVolume;