% 1 water per volume el.
% 16+2 g/mol -> 1000 * 1/18 mol wasser/l
n =1000/(18*(100*1000*1000)^3)*6.023*10^23 %H20 pro nm^3
pixsize = 1/(n^(1/3))*10; %seitenlaenge in A des cubes, in dem 1 H20
