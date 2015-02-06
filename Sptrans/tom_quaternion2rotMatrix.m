%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%http://vamos.sourceforge.net/matrixfaq.htm#Q47
%Assuming that a quaternion has been created in the form:
% 
%     Q = |X Y Z W|
% 
%   Then the quaternion can then be converted into a 4x4 rotation
%   matrix using the following expression:
% 
% 
%         |       2     2                                |
%         | 1 - 2Y  - 2Z    2XY - 2ZW      2XZ + 2YW     |
%         |                                              |
%         |                       2     2                |
%     M = | 2XY + 2ZW       1 - 2X  - 2Z   2YZ - 2XW     |
%         |                                              |
%         |                                      2     2 |
%         | 2XZ - 2YW       2YZ + 2XW      1 - 2X  - 2Y  |
%         |                                              |
% 
% 
%   If a 4x4 matrix is required, then the bottom row and right-most column
%   may be added.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rot = rot3D_quaternion2rotMatrix(q)

    x = q(1);
    y = q(2);
    z = q(3);
    w = q(4);

    x2 = x * x;
    y2 = y * y;
    z2 = z * z;
    xy = x * y;
    xz = x * z;
    yz = y * z;
    wx = w * x;
    wy = w * y;
    wz = w * z;

    rot = [[ 1 - 2*(y2 + z2),     2*(xy - wz),     2*(xz + wy), 0 ]; ...
           [     2*(xy + wz), 1 - 2*(x2 + z2),     2*(yz - wx), 0 ]; ...
           [     2*(xz - wy),     2*(yz + wx), 1 - 2*(x2 + y2), 0 ]; ...
           [               0,               0,               0, 1 ]];
end