function rotCount = tom_av3_create_rot_file(zRotIncrement, tiltIncrement ,tiltThreshold, resultFile)
phi = 0;
psi = 0;
the = 0;
rotCount=1;
while phi < 90
    while psi <= tiltThreshold
        while the <= tiltThreshold
            ang(1,rotCount) = phi/180*pi;
            ang(2,rotCount) = psi/180*pi;
            ang(3,rotCount) = the/180*pi;
            rotCount = rotCount+1;
            the = the +tiltIncrement;
        end
        the = 0;
        psi = psi +tiltIncrement;
    end
    psi = 0;
    phi = phi+zRotIncrement;
end
tom_emwrite(resultFile,ang);