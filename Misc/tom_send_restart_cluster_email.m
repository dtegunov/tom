function tom_send_restart_cluster_email(nr)
%TOM_SEND_RESTART_CLUSTER_EMAIL creates ...
%
%   tom_send_restart_cluster_email(nr)
%
%PARAMETERS
%
%  INPUT
%   nr                  Nachricht
%  
%  OUTPUT
%
%EXAMPLE
%   ... = tom_send_restart_cluster_email(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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


sendmail('rz-linux@biochem.mpg.de',['Bitte cluster ' num2str(nr) ' mit Matlab worker neu starten. Danke.'],['Bitte cluster und Matlab worker neu starten. Danke']);

%sendmail('rz-linux@biochem.mpg.de','Bitte cluster06 und Matlab worker neu starten. Danke.',['Bitte cluster06 und Matlab worker neu starten. Danke']);

%test

%sendmail('nickell@biochem.mpg.de','Bitte cluster und Matlab worker neu starten. Danke.',['Bitte cluster und Matlab worker neu starten. Danke']);
