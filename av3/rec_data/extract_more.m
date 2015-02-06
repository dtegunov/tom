function motl = extract_refined_new(motl)
% motl = extract_refined_new(motl)
%   function extracts particles form tomograms - the clue is that the shift
%   in columns 11, 12 and 13 are taken into account
%   07/25/03 FF
%   last change 11/05/04 FF




itomo=0;
for ind = 1:size(motl,2)
    itomotmp = motl(5,ind);
    if (itomotmp > itomo)
        itomo = itomotmp;
        switch itomo
        case 1
            filename = '';
        case 2
            filename = '../../DN2/DN2_REC_3';
        case 3
            filename = '../../DN3/DN3_VOL_1';
        case 4
            filename = '../../DN4/DN4_VOL_PV_1';
        case 5
            filename = '';
        case 6
            filename = '../../DN6/DN6_VOL_5';
        case 8
            filename = '../../DN8/VDN8_1';
        case 9
            filename = '../../DN9/VDN9_1';
        case 10
            filename = '../../DN10/DN10_VOL_1';
        case 11
            filename = '../../../DN11/DN11_VOL_1';
        case 12
            filename = '../../DN12/DN12_VOL_1';
        case 13
            filename = '../../../DN13/DN13_VOL_1';
        case 15
            filename = '../../../DN15/DN15_VOL_2';
        end;
        vol = tom_emread(filename);
    end;
        xcent = motl(8,ind);
        ycent = motl(9,ind);
        zcent = motl(10,ind);
        xshift = motl(11,ind);motl(11,ind)=0;motl(14,ind)=0;
        yshift = motl(12,ind);motl(12,ind)=0;motl(15,ind)=0;
        zshift = motl(13,ind);motl(13,ind)=0;motl(16,ind)=0;
        psi = pi/180*motl(18,ind);
        theta = pi/180*motl(19,ind);
        % Koordinaten im ungedrehten Koordinatensystem koennen durch
        % Rotation erzeugt werden - die shift sind im (zurueck!) gedrehten
        % Koordinatensystem bestimmt, d.h. die Koordinaten sind bzgl des um
        % -theta, -psi rotierten Partikels bestimmt. Die ungedrehten shift
        % koennen deshalb durch Rotation bestimmt werden.
        matr = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0;0 0 1];
        matr = matr*[1 0 0 ; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
        rshift = [xshift yshift zshift];
        rshift = round(matr*rshift');
        xcent = xcent + rshift(1);
        ycent = ycent + rshift(2);
        zcent = zcent + rshift(3);
        %part = vol.Value(xcent-32:xcent+31,ycent-32:ycent+31,zcent-32:zcent+31);
        part = tom_red(vol.Value,[xcent-32 ycent-32 zcent-32], [64 64 64]);
        motl(8,ind)=xcent;motl(9,ind)=ycent;motl(10,ind)=zcent;
        [meanv, maxv, minv, stdv, varv] = tom_dev(part,'noinfo');
        part = (part-meanv)/stdv;
        tom_dspcub(tom_limit(double(tom_rotate(part,[-psi*180/pi,0,-theta*180/pi])),-2,2));
        indx = motl(4,ind);
        name = ['../PAR_EQ_' num2str(indx) '.em'];
        tom_emwrite(name, part);
        disp([name ' written']);
end;
