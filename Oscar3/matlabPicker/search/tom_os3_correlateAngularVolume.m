function ac = tom_os3_correlateAngularVolume(angularCorrelationVolume, patternFingerprint)

sizeAC = size(angularCorrelationVolume);

ac = zeros(sizeAC(1),sizeAC(2));

patternSTD = std(patternFingerprint(:));

for x=1:sizeAC(1)
    for y=1:sizeAC(2)

        volumeSTD = std(angularCorrelationVolume(x,y,:));

%          if(patternSTD*0.7 < volumeSTD &&  volumeSTD < patternSTD *1.3)
            [a ac(x,y)] = tom_peak(tom_corr(patternFingerprint,squeeze(angularCorrelationVolume(x,y,:))','norm'));
%          end;
                
    end;
end;