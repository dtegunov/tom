function rebuild_mex(flags, idx_compile, do_64bit)

clear functions

mfilepath = fileparts(mfilename('fullpath'));

% Change to the libs directory where to save the binaries.
oldpath = pwd();
cd(fullfile(mfilepath, '.libs'));

% Set the include path.
include = {['-I' fullfile(mfilepath, '/include')]};

% get the name of the lapack-flag (-llapack).
if (isunix())
	lapacklib = '-lmwlapack';
elseif (ispc())
	lapacklib = 'C:\Program Files\MATLAB\R2007b\extern\lib\win64\microsoft\libmwlapack.lib';
	if (~exist(lapacklib,'file'))
		lapacklib = 'C:\Program Files\MATLAB\R2007b\extern\lib\win32\microsoft\libmwlapack.lib';
	end;
	if (~exist(lapacklib,'file'))
		error('lapack-library from matlab not found.');
	end;
end;


params = { '-v', '-largeArrayDims' };
%params{end+1} = '-g';

if (exist('flags','var'))
    params{end+(1:length(flags))} = flags{:};
end;

idx_compileb = false(1,8);
idx_compileb(idx_compile) = true;
idx_compile = idx_compileb;


if (ispc())
    dst_dir = 'W:\tom_dev\mex8a';
    dst_dir = '';
else
    dst_dir = '/fs/pool/pool-bmsan-apps/tom_dev/mex8a/';
end;

if (ispc())
    if (~exist('do_64bit', 'var') || do_64bit)
        mex_extension = '.mexw64';
    else
        mex_extension = '.mexw32';
    end;
else
    if (~exist('do_64bit', 'var') || do_64bit)
        mex_extension = '.mexa64';
    else
        mex_extension = '.mexa32';
    end;
end;

if (0 && isunix())
    if (system(['cd ' mfilepath '; bash -c ''NDEBUG=yes make shared''']))
        error('make not successfull');
    end;
    copyfile(fullfile(mfilepath, '.libs/libtomc.so*'), dst_dir, 'f');
    fileattrib(fullfile(dst_dir, 'libtomc.so*'), '+w');
end;

use_shared = 0 && isunix();
do_copy_mfiles = 0 && isunix();
dst_dir = '';

i = 1;
try

    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emread.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
        else
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emread.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;
        
        if (~isempty(dst_dir))
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_emread' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_emread' mex_extension]), '+w');
            if (do_copy_mfiles)
                copyfile(fullfile(mfilepath, 'src/mex/tom_mex_emread.m'), dst_dir, 'f');
                copyfile(fullfile(mfilepath, 'src/mex/tom_emreadc3.m'), dst_dir, 'f');
                copyfile(fullfile(mfilepath, 'src/mex/tom_emheader_data_2_struct.m'), dst_dir, 'f');
                fileattrib(fullfile(dst_dir, 'tom_mex_emread.m'), '+w');
                fileattrib(fullfile(dst_dir, 'tom_emreadc3.m'), '+w');
                fileattrib(fullfile(dst_dir, 'tom_emheader_data_2_struct.m'), '+w');
            end;
        end;
    end;
    
    i = i + 1;
    
    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
        else
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;
        
        if (~isempty(dst_dir))
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_emwrite' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_emwrite' mex_extension]), '+w');
            if (do_copy_mfiles)
                copyfile(fullfile(mfilepath, 'src/mex/tom_mex_emwrite.m'), dst_dir, 'f');
                copyfile(fullfile(mfilepath, 'src/mex/tom_emwritec3.m'), dst_dir, 'f');
                copyfile(fullfile(mfilepath, 'src/mex/tom_emheader_struct_2_data.m'), dst_dir, 'f');
                fileattrib(fullfile(dst_dir, 'tom_mex_emwrite.m'), '+w');
                fileattrib(fullfile(dst_dir, 'tom_emwritec3.m'), '+w');
                fileattrib(fullfile(dst_dir, 'tom_emheader_struct_2_data.m'), '+w');
            end;
        end;
    end;
    i = i + 1;

    
    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_append_stack.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
        else
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_append_stack.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;
        
        if (~isempty(dst_dir))
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_emwrite_append_stack' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_emwrite_append_stack' mex_extension]), '+w');
            if (do_copy_mfiles)
                copyfile(fullfile(mfilepath, 'src/mex/tom_mex_emwrite_append_stack.m'), dst_dir, 'f');
                fileattrib(fullfile(dst_dir, 'tom_mex_emwrite_append_stack.m'), '+w');
            end;
        end;
    end;
    i = i + 1;
    
    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_paste.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
        else
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_paste.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;

        if (~isempty(dst_dir))
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_emwrite_paste' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_emwrite_paste' mex_extension]), '+w');
            if (do_copy_mfiles)
                copyfile(fullfile(mfilepath, 'src/mex/tom_mex_emwrite_paste.m'), dst_dir, 'f');
                fileattrib(fullfile(dst_dir, 'tom_mex_emwrite_paste.m'), '+w');
            end;
        end;
    end;
    i = i + 1;

    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_new.c'), include{:}, ['-Wl,-rpath,' dst_dir], ['-L' dst_dir], '-ltomc');
        else
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite_new.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;

        if (~isempty(dst_dir))
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_emwrite_new' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_emwrite_new' mex_extension]), '+w');
            if (do_copy_mfiles)
                copyfile(fullfile(mfilepath, 'src/mex/tom_mex_emwrite_new.m'), dst_dir, 'f');
                fileattrib(fullfile(dst_dir, 'tom_mex_emwrite_new.m'), '+w');
            end;
        end;
    end;
    i = i + 1;
    

    if (idx_compile(i))
        if (use_shared)
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_corr3d.cpp'), include{:}, libs{:}, '-Wl,-rpath,/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-L/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-ltomc', '-lfftw3 -lfftw3f');
        else
            libs = {'core__fftw_plan.o'
                    'core__interpol.o'
                    'core__io.o'
                    'core__rotate.o'
                    'core__transform.o'
                    'core__transform_cpp.o'
                    'core__volume.o'
                    'core__volume_container.o'
                    'core__volume_fcn.o'
                    'core__wedge.o'
                    'corr__corr.o'
                    'corr__correlation_handler.o'};        
            for (ii=1:length(libs)) libs{ii} = fullfile(mfilepath, ['/.libs/' libs{ii}]); end;
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_corr3d.cpp'), libs{:}, '-lfftw3 -lfftw3f', include{:});
        end;
        if (~isempty(dst_dir))
            system(['cd ' mfilepath '; cp .libs/tom_mex_corr3d.mexa64 src/mex/tom_mex_corr3d.m src/mex/tom_av3_corr3d.m /fs/pool/pool-bmsan-apps/tom_dev/av3_new/; cd /fs/pool/pool-bmsan-apps/tom_dev/av3_new/; chmod 777 tom_mex_corr3d.mexa64 tom_mex_corr3d.m tom_av3_corr3d.m']);
        end;
    end;
    i = i+1;
    
    
    if (idx_compile(i))
        if (use_shared)
            error('TODO...');
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_parse_av4.cpp'), include{:}, libs{:}, '-Wl,-rpath,/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-L/fs/pool/pool-bmsan-apps/tom_dev/mex8a', '-ltomc', '-lfftw3 -lfftw3f');
        else
            libs = {    'boost__md5.o'
                        'core__interpol.o'
                        'core__io.o'
                        'core__transform.o'
                        'core__transform_cpp.o'
                        'core__volume.o'
                        'core__volume_fcn.o'
                        'core__wedge.o'
                        'corr__config_files.o'
                        'corr__filename_generators.o'
                        'helper__snippets.o'
                   };        
            for (i=1:length(libs)) libs{i} = fullfile(mfilepath, ['/.libs/' libs{i}]); end;
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_parse_av4.cpp'), libs{:}, '-lfftw3 -lfftw3f', include{:});
        end;
    end;
    i = i+1;

    
    
    if (0)
        if (idx_compile(2))
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_proj3d.c'), fullfile(mfilepath, '/src/tom/core/transform.c'), fullfile(mfilepath, '/src/tom/core/interpol.c'), include{:}, lapacklib);
        end;
        if (idx_compile(3))
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_rotate.c'), fullfile(mfilepath, '/src/tom/core/rotate.c'), fullfile(mfilepath, '/src/tom/core/interpol.c'), include{:});
            copyfile(fullfile(mfilepath, ['.libs/tom_mex_rotate' mex_extension]), dst_dir, 'f');
            fileattrib(fullfile(dst_dir, ['tom_mex_rotate' mex_extension]), '+w');
        end;
        if (idx_compile(6))
            libs = {'core__fftw_plan.o'
                    'core__interpol.o'
                    'core__io.o'
                    'core__rotate.o'
                    'core__transform.o'
                    'core__transform_cpp.o' 
                    'core__volume.o'
                    'core__volume_container.o'
                    'core__volume_fcn.o'
                    'core__wedge.o'
                    'corr__corr.o'
                    'corr__correlation_handler.o'};        
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_rotate2.cpp'), include{:}, libs{:}, '-lfftw3 -lfftw3f');
            if (isunix())
                system(['cd ' mfilepath '; cp .libs/tom_mex_rotate2.mexa64 /fs/sandy01/lv03/pool/bmsan/apps/tom_dev/Sptrans/tom_rotatec2.mexa64; chmod 777 /fs/sandy01/lv03/pool/bmsan/apps/tom_dev/Sptrans/tom_rotatec2.mexa64']);
            end;
        end;
        if (idx_compile(7))
            mex(params{:}, fullfile(mfilepath, '/src/mex/pam_mysql.cpp'), include{:}, '-I/usr/include/mysql', '-lmysqlclient');
        end;
        if (idx_compile(8))
            mex(params{:}, fullfile(mfilepath, '/src/mex/tom_mex_emwrite.c'), fullfile(mfilepath, '/src/tom/core/io.c'), include{:});
        end;
    end;
catch
    % Back to the old directory
    cd(oldpath);
    rethrow(lasterror);
end;
cd(oldpath);


