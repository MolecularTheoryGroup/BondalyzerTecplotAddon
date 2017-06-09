#!/bin/sh

#
#  Distribution looks like:
#
#  README.txt
#  doc/
#      HTML/
#           tecint.css
#           tecint.html (plus gif files ande tecint.css)
#           (a bunch of .gif files)
#      PDF/
#          tecint.pdf
#  include/
#          tecint.h
#  lib/
#      <platform>/ (all supported Unix platforms, plus win32 for windows)
#                 libtecint.a

if test `whoami` != "build" ; then
  echo "This script must be run as build."
  exit 1
fi

if test ! -x /global/scripts/SetupTecEnv ; then
  echo "Cannot execute /global/scripts/SetupTecEnv"
  exit 1
fi

. /global/scripts/SetupTecEnv

if test -z "$TECPLOTBUILDSANDBOX" ; then
  echo "Please set TECPLOTBUILDSANDBOX."
  exit 1
fi

# Tecplot build scripts are used to get the
# current model list and build machine names.

RAWDISTDIR=/rawdist/tecint
TECDEV=$TECPLOTBUILDSANDBOX/v11.0-4
BUILDUTILDIR=$TECDEV/autobuild/unixutil
AllPlatforms=`awk '{print $1}' $BUILDUTILDIR/MasterPlatList | sed 's/^u[a-z.0-9]*//'`

if test ! -d $BUILDUTILDIR ; then
  echo "Directory $BUILDUTILDIR does not exist."
  exit 1
fi

if test ! -f $BUILDUTILDIR/MasterPlatList ; then
  echo "File $BUILDUTILDIR/MasterPlatList does not exist."
  exit 1
fi

# Create a subdirectory
MakeDir()
{
  echo MakeDir $1

  if [ ! -d $1 ] ; then
    mkdir -p $1
    if [ $? -ne 0 ] ; then
      exit $?
    fi
    chmod 755 $1
  fi
}

# Copy a single file to the named location
ListAndCopyFile()
{
  echo Copy `/bin/ls -l $1`
  /bin/cp -pf $1 $2
  if [ $? -ne 0 ] ; then
    exit $?
  fi
  chmod 644 $2
}

# Wipe out everything there
rm -r $RAWDISTDIR

# Make the basic directories
for dir in $RAWDISTDIR/doc/PDF $RAWDISTDIR/doc/HTML $RAWDISTDIR/include $RAWDISTDIR/lib ; do
  MakeDir $dir
done

for platdir in $AllPlatforms ; do
  dir=$RAWDISTDIR/lib/$platdir
  MakeDir $dir
done

MakeDir $RAWDISTDIR/lib/win32

# Copy the documentation, README.txt and include file
ListAndCopyFile doc/PDF/tecint.pdf $RAWDISTDIR/doc/PDF/tecint.pdf
for file in doc/HTML/* ; do
  ListAndCopyFile $file $RAWDISTDIR/$file
done
ListAndCopyFile README.txt $RAWDISTDIR/README.txt

# Copy the include file
ListAndCopyFile lib/tecint.h $RAWDISTDIR/include/tecint.h

# Copy the libraries
for platdir in $AllPlatforms ; do
  file=lib/$platdir/release/libtecint.a
  ListAndCopyFile $file $RAWDISTDIR/lib/$platdir/`basename $file`
done

if test -f lib/Release/tecint.lib ; then
  ListAndCopyFile lib/Release/tecint.lib $RAWDISTDIR/lib/win32/tecint.lib
  echo
  echo "*****************************************************************************"
  echo "For the distribution, tar and compress $RAWDISTDIR, and copy to"
  echo "the distribution directory."
  echo "*****************************************************************************"
  echo
else
  echo
  echo "*****************************************************************************"
  echo "Remember to copy the WIN32 Release build to $RAWDISTDIR/lib/win32/tecint.lib,"
  echo "then tar and compress $RAWDISTDIR, and copy to the distribution directory."
  echo "*****************************************************************************"
  echo
fi
