#!/bin/bash



# compile boost under AIX, IBM XL Compiler:
# http://www-1.ibm.com/support/docview.wss?rs=2239&context=SSJT9L&dc=DA420&dc=DA410&dc=DA440&dc=DA430&uid=swg27006911&loc=en_US&cs=utf-8&lang=en

BOOSTVRSN=1_34_0

PREFIX=~/usr/
SRCDIR="$PREFIX/src/"
PWD_=`pwd`
ARCHIVEDIR="$PWD_/../"

# Create the destination directory and extract the library.
mkdir -p "$SRCDIR"
cd "$SRCDIR"
rm -rf "./boost_$BOOSTVRSN/"
gunzip -d -c "$ARCHIVEDIR/boost_$BOOSTVRSN.tar.gz" | tar -xf -


# Apply the patch.
cp "$PWD_/boost_modfile_$BOOSTVRSN.txt" "./boost_$BOOSTVRSN.patch"
patch -p0 < "./boost_$BOOSTVRSN.patch"





sh ./configure "--prefix=$PREFIX" --with-libraries=filesystem


#change lines in ./tools/build/v2/tools/vacpp.jam to:
#           flags vacpp VA_C_COMPILER  <threading>single : xlc -DBOOST_NO_IS_ABSTRACT ;
#           flags vacpp VA_C_COMPILER  <threading>multi : xlc_r -DBOOST_NO_IS_ABSTRACT ;
#           flags vacpp VA_CXX_COMPILER <threading>single : xlC -DBOOST_NO_IS_ABSTRACT ;
#           flags vacpp VA_CXX_COMPILER <threading>multi : xlC_r -DBOOST_NO_IS_ABSTRACT ;
#           ar -q $(ARFLAGS:E="") "$(<)" "$(>)"

#change lines 732 in /boost/filesystem/operations.hpp
#         #if 0
#         static const boost::filesystem::basic_directory_iterator<Path> end_itr;
#         unsigned long count = 1;
#         if ( !boost::filesystem::is_symlink( ph ) // don't recurse symbolic links
#           && boost::filesystem::is_directory( ph ) )
#         {
#           for ( boost::filesystem::basic_directory_iterator<Path> itr( ph );
#                 itr != end_itr; ++itr )
#           {
#             count += remove_all_aux( itr->path() );
#           }
#         }
#         boost::filesystem::remove( ph );
#         return count;
#         #else
#         throw std::exception(__FILE__ ": there was an error during compilation. This function is now errouneous :)");
#         #endif

./tools/jam/src/bin.aixppc/bjam --user-config=user-config.jam --with-filesystem address-model=64 -n




