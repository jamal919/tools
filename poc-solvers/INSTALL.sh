#!/bin/bash

LOG_NAME=.INSTALL.log
DIR_NAME=cmade

print_help(){
cat <<EOF
USE
  ./INSTALL.sh [OPTIONS] [ -DVAR1=val1 [ -DVAR2=val2 .... ] ]

DESCRIPTION
  Compile PaStiX, if necessary, and the poc-solvers.
  The build directory is $DIR_NAME by default.
  The ouput is saved to <build directory>/$LOG_NAME

OPTIONS
  -h,--help  print this help
  -d,--debug  set debug mode
  -D,--separate-debug  set debug mode with separate cmade-dbg build directory
  -1,--sequential : set sequential mode. Implies --pastix-type none
  -S,--separate-sequential : set sequential mode (see above) with separate cmade-nompi build directory.
  --redo : clean-up output dir and, if any, zdhips
  -s,--scotch-dir followed by the directory where to search for SCOTCH. Default: /usr
  -p,--pastix-type : followed by auto, distribution, custom (default and recommended) or none
  --zdhips  compile zdhips
EOF
}

#*******************************************************************************
# parse options

# backup arguments
args=("$@")

# set defaults
SCOTCH_HOME=/usr
CMAKE_VAR_ARGS=
PASTIX_COMPILE_TARGET=
PASTIX_TYPE=custom
REDO=
ZDHIPS_MAKE_OPTS=
ZDHIPS=

testIfFollowed(){
  [ "$2" -gt 1 ] && return
  echo "*** option $1 not followed by anything ***"
  print_help
  exit -1
  }

while [ $# -gt 0 ]; do
  case "$1" in
  -h|--help)
    print_help
    exit
    ;;
  -D|--separate-debug)
    DIR_NAME=cmade-dbg
    ;&
  -d|--debug)
    CMAKE_VAR_ARGS="$CMAKE_VAR_ARGS -DCMAKE_BUILD_TYPE=Debug"
    PASTIX_COMPILE_TARGET=debug
    ZDHIPS_MAKE_OPTS="OPTFLAGS=-g"
    shift
    ;;
  -S|--separate-sequential)
    DIR_NAME=cmade-nompi
    ;&
  -1|--sequential)
    # From the output of echo `sed -re 's/set\(([^ ]+) YES .* enable.* MPI.*/-D\1=NO/;t;d' src/CMakeLists.txt`
    CMAKE_VAR_ARGS="$CMAKE_VAR_ARGS -DHIPS=NO -DPARPACK=NO -DMUMPS=NO"
    PASTIX_TYPE=none
    shift
    ;;
  --redo)
    REDO=yes
    shift
    ;;
  -s|--scotch-dir)
    testIfFollowed "$1" "$#"
    SCOTCH_HOME="$2"
    shift 2
    ;;
  -p|--pastix-type)
    testIfFollowed "$1" "$#"
    PASTIX_TYPE="$2"
    shift 2
    ;;
  --zdhips)
    ZDHIPS=yes
    shift
    ;;
  -D*)
    CMAKE_VAR_ARGS="$CMAKE_VAR_ARGS $1"
    shift
    ;;
    esac
  done

#*******************************************************************************
set -x

ncpu=`perl -e 'while(<>){s/^processor.*: *// or next;$n=$_}print $n+1' /proc/cpuinfo` #number of CPU

# backup path of this file
THIS_PATH=`readlink -f $0`

# record version to log
hg -R `dirname $THIS_PATH` parent --template "Mercurial revision {rev}:{node|short} of {date|isodate}\n"

# go to a directory parent to this file
cd `dirname $0`/..

#*******************************************************************************
# OUTPUTS

# directory
OUT_DIR="poc-solvers/$DIR_NAME"
mkdir -p $OUT_DIR

# log file MUST START WITH .
LOG_PATH=$OUT_DIR/$LOG_NAME
if [ "`readlink -f /proc/$$/fd/0`" == "`readlink -f /proc/$$/fd/1`" ];then
  echo "> $THIS_PATH ${args[*]}" >$LOG_PATH
  $THIS_PATH ${args[@]} |& tee -a $LOG_PATH
  exit
  fi

# optional clean-up
if [ "$REDO" ];then
  # make sure ...
  true >"$OUT_DIR"/deleteThis
  # ... this completes to at least one file
  rm -r "$OUT_DIR"/*
  [ "$ZDHIPS" ] &&
    make -C zdhips clean
  fi

#*******************************************************************************
# set SCOTCH_INC and SCOTCH_LIB as necessary

findAndSelect(){
  f=${@/#/-o -name }
  WORDS=(`find $SCOTCH_HOME -false $f`)
  
  case ${#WORDS[@]} in
  0)
    set +x
    echo "$f" NOT FOUND. CHECK YOUR INSTALLATION OF SCOTCH.
    exit 1
    ;;
  1)
    WORD="$WORDS"
    ;;
  *)
    select WORD in "${WORDS[@]}";do break;done
    ;;
    esac
  
  WORD="`dirname $WORD`"
  }

SCOTCH_LIB=
SCOTCH_INC=

findScotch(){
  
  [ "x$SCOTCH_LIB" != x -a "x$SCOTCH_INC" != x ] &&
    return
  
  findAndSelect libscotch.{so,a}
  SCOTCH_LIB="$WORD"

  findAndSelect scotch.h
  SCOTCH_INC="$WORD"
  }

#*******************************************************************************
case "$PASTIX_TYPE" in
a*)
  PASTIX_TYPE=auto
  ;;
d*)
  PASTIX_TYPE=distribution
  ;;
c*)
  PASTIX_TYPE=custom
  ;;
n*)
  PASTIX_TYPE=none
  ;;
*)
  set +x
  echo "*** ARGUMENT ERROR: PASTIX_TYPE=$PASTIX_TYPE ***"
  print_help
  exit -1
  ;;
  esac

if [ "$PASTIX_TYPE" = auto ];then
  if which pastix-conf;then
    PASTIX_TYPE=distribution
  else
    PASTIX_TYPE=custom
    fi
  fi

case "$PASTIX_TYPE" in
distribution)
  true "You have PaStiX"
  ;;
custom)
  # See https://gforge.inria.fr/frs/?group_id=186
  version=5.2.2.22;remodeDirName=35070
  baseName=pastix_$version
  if [ \! -f $baseName/install/pastix-conf ];then # if you do not have PASTIX
    true "======================================================="
    true "INSTALLING PaStiX $version"
    if [ \! -f $baseName.tar.bz2 ];then
      wget https://gforge.inria.fr/frs/download.php/file/$remodeDirName/$baseName.tar.bz2 ||
        exit -1
      fi
    
    if [ \! -f $baseName/src/utils/src/pastix-conf.sh ];then
      tar xf $baseName.tar.bz2 # extract it
      # patch it for gcc 7
      sed -i~ -re 's/^( *if \(!?)(isnan.+)(\) \{ *$)/\11 \/*\2*\/\3/' pastix_5.2.2.22/src/sopalin/src/variable_csc.c
      fi
    # doc in $baseName/INSTALL.txt and $baseName/doc/refcard/refcard.pdf
    
    if [ \! -f $baseName-config.in ];then # If you do not have $baseName-config.in
      # create it and save it out of the PASTIX directory.
      # It will be usefull for when you want to build PASTIX again.
      
      findScotch
      
      # see the comments in the example config file
      sed -r \
        -e "s|(^ *SCOTCH_HOME \?= *).*|\1$SCOTCH_HOME|"\
        -e "s|(^ *SCOTCH_INC \?= *).*|\1$SCOTCH_INC|" \
        -e "s|(^ *SCOTCH_LIB \?= *).*|\1$SCOTCH_LIB|" \
        -e "161,162 s|^#||" \
        -e "164,168 s|^|#|" \
        $baseName/src/config/LINUX-GNU.in >$baseName-config.in
      
      grep -E '^ *SCOTCH_' $baseName-config.in
      set +x
      read -n 1 -p 'IF YOU WANT TO CHECK OR EDIT $baseName-config.in, PRESS ANYTHING BUT Enter TO EXIT:'
      if [ "$REPLY" ]; then
        echo
        exit
        fi
      fi
    
    set -x
    cp -pv $baseName-config.in $baseName/src/config.in #copy the right config file
    make -C $baseName/src all drivers -j$ncpu $PASTIX_COMPILE_TARGET || {
      true make sure you have hwloc properly installed
      exit -1
      }
    fi
  
  PATH=$PWD/$baseName/install:$PATH
  export CMAKE_INCLUDE_PATH=$PWD/$baseName/install${CMAKE_INCLUDE_PATH:+:$CMAKE_INCLUDE_PATH}
  export LD_LIBRARY_PATH=$PWD/$baseName/install${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  ;;
none)
  CMAKE_VAR_ARGS="$CMAKE_VAR_ARGS -DPASTIX=NO"
  ;;
*)
  set +x
  echo "PROGRAMMING ERROR: PASTIX_TYPE=$PASTIX_TYPE"
  exit 8
  esac

[ "$PASTIX_TYPE" != none ] && pastix-conf

#*******************************************************************************
if [ "$ZDHIPS" ];then
  
  dhips=zdhips/dhips
  zhips=zdhips/zhips
  
  if [ \! -f $dhips/LIB/libhips.a -o \! -f $zhips/LIB/libhips.a -o \! -f $dhips/LIB/dhips.h -o \! -f $zhips/LIB/zhips.h ];then
    findScotch
    mkdir -p zdhips
    ln -s ../poc-solvers/Makefile.zdhips zdhips/Makefile
    make -C zdhips $ZDHIPS_MAKE_OPTS SCOTCHINCDIR=$SCOTCH_INC SCOTCHLIBDIR=$SCOTCH_LIB -j$ncpu
    fi
  
  fi

#*******************************************************************************
if [ \! -f $OUT_DIR/lib/poc-solvers.LDFLAGS ];then # if you do not have the poc-solvers
  # INSTALL the poc-solvers
  cd $OUT_DIR
  CMAKE_LIBRARY_PATH=$LD_LIBRARY_PATH cmake $CMAKE_VAR_ARGS ../src &&
    make -j$ncpu install ||
    cat <<-EOF
	OK, it failed, but make sure you have LD_LIBRARY_PATH and CMAKE_INCLUDE_PATH environment variables set properly, as most people rarely do. Currently, we have :
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH
	CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH
	EOF
  fi
