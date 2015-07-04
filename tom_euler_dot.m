function [ cosine ] = tom_euler_dot( angles1, angles2 )

m1 = tom_mat_euler(angles1);
m2 = tom_mat_euler(angles2);

t1 = m1 * [0; 0; 1];
t2 = m2 * [0; 0; 1];

cosine = sum(t1.*t2);

end

