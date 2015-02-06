% This script will compile all the C files
files=dir('*.c');
for i=1:length(files)
    mex(files(i).name,'-v');
end
