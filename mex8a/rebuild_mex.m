% Compiles the mex functions from the current directory...
%
% Parameters:
%  tag: Name of the target to compile. Can be either a string or a cell
%       array of strings. These are the names of the functions.
%  flags: Cell array of strings with the compiler flags. Directly passed to
%       mex.
%  libname: name of the library libtomc, to link agains.
%       A filename ending with .a is a static library (means static
%       linkage). File extension .so is a shared library, means dynamic
%       linkage.
%       The directory in which libname can be found is used for the 
%       -rpath option for the linker (in case of dynamic linking).
%       In this case libname should be an absolute filename because
%       otherwise the absolute name is obtained using fileattrib, which
%       resolves symbolic links.
function rebuild_mex(tag, flags, libname)

clear functions

mfilepath = fileparts(mfilename('fullpath'));


% Set the include path.
include = {['-I' fullfile(mfilepath, '../../include')]};

% get the name of the lapack-flag (-llapack).
% if (isunix())
% 	lapacklib = '-lmwlapack';
% elseif (ispc())
% 	lapacklib = 'C:\Program Files\MATLAB\R2007b\extern\lib\win64\microsoft\libmwlapack.lib';
% 	if (~exist(lapacklib,'file'))
% 		lapacklib = 'C:\Program Files\MATLAB\R2007b\extern\lib\win32\microsoft\libmwlapack.lib';
% 	end;
% 	if (~exist(lapacklib,'file'))
% 		error('lapack-library from matlab not found.');
% 	end;
% end;

if (~exist('tag','var'))
    tag = {};
elseif (ischar(tag))
    tag = {tag};
end;

params = { '-v', '-largeArrayDims' };

if (exist('flags','var'))
    params{end+(1:length(flags))} = flags{:};
end;

% Prepare the library name...
if (isunix()) 
    if (~exist('libname','var'))
        libname = '../../.libs/libtomc.a';
    end;
    if (~exist(libname, 'file'))
        error(['The library ' libname ' was not found. did you compile the source?']);
    end;
    if (length(libname)>2 && strcmp(libname(end-1:end), '.a'))
        libparams = {libname};
    elseif (length(libname)>3 && strcmp(libname(end-2:end), '.so'))
        s = libname;
        if (libname(1) ~= '/')
            [s1, s] = fileattrib(libname);
            clear s1;
            s = s.Name;        
        end;
        [d1,d2,d3,d4] = fileparts(s);
        if (length(d2)<4 || ~strcmp(d2(1:3), 'lib'))
            error(['a shared library must have the name lib????.so (' libname ').']);
        end;
        libparams = {['-Wl,-rpath,' d1], ['-L' d1], ['-l' d2(4:end)] };
        clear d1 d2 d3 d4;
    else
        error(['unrecognized library name ' libname]);
    end;
    fftw_libs = {'-lfftw3', '-lfftw3f'};
else
    libparams = {'../tom/io/io.c'};
    fftw_libs = {fullfile(matlabroot, 'bin', 'win32', 'libfftw3.lib'), fullfile(matlabroot, 'bin', 'win32', 'libfftw3f.lib')};
end;


% if (use_shared)
%     mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emread.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
% else
%     mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emread.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
% end;

if (isempty(tag) || any(strcmp(tag, 'tom_mex_emread')))
    mex(params{:}, fullfile(mfilepath, 'tom_mex_emread.c'), include{:}, libparams{:});
    %WINDOWS: mex tom_mex_emread.c -I../../include ../tom/io/io.c -v
end
if (isempty(tag) || any(strcmp(tag, 'tom_mex_emwrite'))) 
    mex(params{:}, fullfile(mfilepath, 'tom_mex_emwrite.c'), include{:}, libparams{:});
    %WINDOWS: mex tom_mex_emwrite.c -I../../include ../tom/io/io.c -v
end
if (isempty(tag) || any(strcmp(tag, 'tom_mex_emwrite_new'))) 
    mex(params{:}, fullfile(mfilepath, 'tom_mex_emwrite_new.c'), include{:}, libparams{:});
    %WINDOWS: mex tom_mex_emwrite_new.c -I../../include ../tom/io/io.c -v
end;
if (isempty(tag) || any(strcmp(tag, 'tom_mex_emwrite_append_stack'))) 
    mex(params{:}, fullfile(mfilepath, 'tom_mex_emwrite_append_stack.c'), include{:}, libparams{:});
    %WINDOWS: mex tom_mex_emwrite_append_stack.c -I../../include ../tom/io/io.c -v
end;
if (isempty(tag) || any(strcmp(tag, 'tom_mex_emwrite_paste'))) 
    mex(params{:}, fullfile(mfilepath, 'tom_mex_emwrite_paste.c'), include{:}, libparams{:});
    %WINDOWS: mex tom_mex_emwrite_paste.c -I../../include ../tom/io/io.c -v
end;
if (isempty(tag) || any(strcmp(tag, 'tom_mex_rotate'))) 
    mex(params{:}, fullfile(mfilepath, 'tom_mex_rotate.cpp'), include{:}, libparams{:}, fftw_libs{:});
end;






% 
% if (idx_compile(i))
%     if (use_shared)
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_corr3d.cpp'), include{:}, libs{:}, '-Wl,-rpath,/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-L/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-ltomc', '-lfftw3 -lfftw3f');
%     else
%         libs = {'core__fftw_plan.o'
%                 'core__interpol.o'
%                 'core__io.o'
%                 'core__rotate.o'
%                 'core__transform.o'
%                 'core__transform_cpp.o'
%                 'core__volume.o'
%                 'core__volume_container.o'
%                 'core__volume_fcn.o'
%                 'core__wedge.o'
%                 'corr__corr.o'
%                 'corr__correlation_handler.o'};        
%         for (ii=1:length(libs)) libs{ii} = fullfile(mfilepath, ['/.libs/' libs{ii}]); end;
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_corr3d.cpp'), libs{:}, '-lfftw3 -lfftw3f', include{:});
%     end;
%     if (~isempty(dst_dir))
%         system(['cd ' mfilepath '; cp .libs/tom_mex_corr3d.mexa64 src/mex/tom_mex_corr3d.m src/mex/tom_av3_corr3d.m /fs/pool/pool-bmsan-apps/tom_dev/av3_new/; cd /fs/pool/pool-bmsan-apps/tom_dev/av3_new/; chmod 777 tom_mex_corr3d.mexa64 tom_mex_corr3d.m tom_av3_corr3d.m']);
%     end;
% end;
% i = i+1;
% 
% 
% if (idx_compile(i))
%     if (use_shared)
%         error('TODO...');
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_parse_av4.cpp'), include{:}, libs{:}, '-Wl,-rpath,/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-L/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-ltomc', '-lfftw3 -lfftw3f');
%     else
%         libs = {    'boost__md5.o'
%                     'core__interpol.o'
%                     'core__io.o'
%                     'core__transform.o'
%                     'core__transform_cpp.o'
%                     'core__volume.o'
%                     'core__volume_fcn.o'
%                     'core__wedge.o'
%                     'corr__config_files.o'
%                     'corr__filename_generators.o'
%                     'helper__snippets.o'
%                };        
%         for (i=1:length(libs)) libs{i} = fullfile(mfilepath, ['/.libs/' libs{i}]); end;
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_parse_av4.cpp'), libs{:}, '-lfftw3 -lfftw3f', include{:});
%     end;
% end;
% i = i+1;
% 
% 
% 
% if (0)
%     if (idx_compile(2))
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_proj3d.c'), fullfile(mfilepath, '/src/tom/core/transform.c'), fullfile(mfilepath, '/src/tom/core/interpol.c'), include{:});
%     end;
%     if (idx_compile(3))
%     end;
%     if (idx_compile(6))
%         libs = {'core__fftw_plan.o'
%                 'core__interpol.o'
%                 'core__io.o'
%                 'core__rotate.o'
%                 'core__transform.o'
%                 'core__transform_cpp.o' 
%                 'core__volume.o'
%                 'core__volume_container.o'
%                 'core__volume_fcn.o'
%                 'core__wedge.o'
%                 'corr__corr.o'
%                 'corr__correlation_handler.o'};        
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_rotate2.cpp'), include{:}, libs{:}, '-lfftw3 -lfftw3f');
%         if (isunix())
%             system(['cd ' mfilepath '; cp .libs/tom_mex_rotate2.mexa64 /fs/sandy01/lv03/pool/bmsan/apps/tom_dev/Sptrans/tom_rotatec2.mexa64; chmod 777 /fs/sandy01/lv03/pool/bmsan/apps/tom_dev/Sptrans/tom_rotatec2.mexa64']);
%         end;
%     end;
%     if (idx_compile(7))
%         mex(params{:}, fullfile(mfilepath, '/src/mex/pam_mysql.cpp'), include{:}, '-I/usr/include/mysql', '-lmysqlclient');
%     end;
%     if (idx_compile(8))
%         mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
%     end;
% end;
% 


