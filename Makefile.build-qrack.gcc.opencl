CWD = $(shell pwd)
MACH = $(shell uname -p)
QRACK_TARGET = $(MACH)

ENABLE_CUDA := 0

ifeq ($(MACH),aarch64)
  QRACK_TARGET = AArch64
endif

DEBUG_BUILD := 0

# Disable Sphinx.
HAVE_SPHINX := 0

CC = /usr/bin/gcc
CXX = /usr/bin/g++
CMAKE = /usr/bin/cmake
GMAKE = /usr/bin/gmake
GMAKE_NUMJOBS = -j4

SRCDIR = qrack
TOPDIR = $(CWD)
TOPSRCDIR = $(TOPDIR)/$(SRCDIR)
TOPBUILDDIR = $(TOPDIR)/qrack-build
LD_LIBRARY_PATH =
CUDA_INSTALL = /usr/local/cuda

PYTHON = /usr/bin/python3.11
PYTHON_EXECUTABLE = $(PYTHON)
PYTHON_LIBRARY = /usr/lib64/libpython3.11.so
PYTHON_LIBRARIES = /usr/lib64/libpython3.11.so
PYTHON_INCLUDE_DIRS = /usr/include/python3.11
SWIG = /usr/bin/swig
SWIG_EXECUTABLE = $(SWIG)

OFLAG = -O3
GFLAG = -g
ENABLE_EXPENSIVE_CHECKS = 0

ifeq ($(DEBUG_BUILD),1)
  OFLAG = -O0
  GFLAG = -g3
  ENABLE_EXPENSIVE_CHECKS = 1
  TOPBUILDDIR = $(TOPDIR)/qrack-build-debug
endif

# Use this prefix when creating internal binaries.
CMAKE_PREFIX = /usr/local
CMAKE_INSTALL_BINDIR = $(CMAKE_PREFIX)/bin
CMAKE_INSTALL_LIBDIR = $(CMAKE_PREFIX)/lib64
CMAKE_INSTALL_LIBEXECDIR = $(CMAKE_PREFIX)/libexec
CMAKE_INSTALL_INCLUDEDIR = $(CMAKE_PREFIX)/include
CMAKE_INSTALL_DATADIR = $(CMAKE_PREFIX)/share
CMAKE_INSTALL_DATAROOTDIR = $(CMAKE_PREFIX)/share

CFLAGS = $(GFLAG) $(OFLAG) -pthread -std=c99 -fno-strict-aliasing
CFLAGS += -Wall -Wcast-align -Wno-long-long -Woverflow
CFLAGS += -Wno-unused-command-line-argument
CFLAGS += -Wno-documentation-unknown-command
CFLAGS += -Wstack-protector -ffunction-sections -fdata-sections
CFLAGS += -fkeep-static-consts -fstack-protector-all
CFLAGS += -fvisibility=default

ifeq ($(DEBUG_BUILD),1)
CFLAGS += -fno-omit-frame-pointer
endif

CFLAGS += -finput-charset=UTF-8 -fPIC
CFLAGS += -fuse-ld=gold
CFLAGS += -Wl,-rpath -Wl,$(CMAKE_INSTALL_LIBDIR)

ifeq ($(DEBUG_BUILD),1)
CFLAGS += -Wl,-O0
else
CFLAGS += -Wl,-O3
endif

CXXFLAGS = $(GFLAG) $(OFLAG) -pthread -std=c++17 -fno-strict-aliasing
CXXFLAGS += -fexceptions -frtti -fstack-protector-all
CXXFLAGS += -Wall -Wcast-align -Wno-long-long -Woverflow
CXXFLAGS += -Wno-redundant-move -finput-charset=UTF-8 -fPIC
CXXFLAGS += -Wno-unused-command-line-argument
CXXFLAGS += -Wstack-protector -fkeep-static-consts
CXXFLAGS += -fdata-sections -ffunction-sections
CXXFLAGS += -fvisibility=default

ifeq ($(DEBUG_BUILD),1)
CXXFLAGS += -fno-omit-frame-pointer
endif

CXXFLAGS += -fuse-ld=gold
CXXFLAGS += -Wl,-rpath -Wl,$(CMAKE_INSTALL_LIBDIR)

ifeq ($(DEBUG_BUILD),1)
CXXFLAGS += -Wl,-O0
else
CXXFLAGS += -Wl,-O3
endif

CPPFLAGS = -D_REENTRANT -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
CPPFLAGS += -D__STDC_LIMIT_MACROS
CPPFLAGS += -D_GNU_SOURCE -D_XOPEN_SOURCE=700
CPPFLAGS += -DCL_HPP_MINIMUM_OPENCL_VERSION=110
CPPFLAGS += -DCL_HPP_TARGET_OPENCL_VERSION=300
CPPFLAGS += -I$(CUDA_INSTALL)/include
CPPFLAGS += -I$(CMAKE_PREFIX)/include

LDFLAGS = -L$(TOPBUILDDIR)/lib
LDFLAGS += -L$(TOPBUILDDIR)/lib64
LDFLAGS += -L$(CMAKE_INSTALL_LIBDIR)
LDFLAGS += -Wl,-rpath -Wl,$(CMAKE_INSTALL_LIBDIR)

ifeq ($(ENABLE_CUDA),1)
LDFLAGS += -L$(CUDA_INSTALL)/lib64
LDFLAGS += -Wl,-rpath -Wl,$(CUDA_INSTALL)/lib64
endif

ifeq ($(DEBUG_BUILD),1)
LDFLAGS += -Wl,-O0
else
LDFLAGS += -Wl,-O3
endif

LDFLAGS += -fuse-ld=gold

ifeq ($(DEBUG_BUILD),1)
  CPPFLAGS += -D_DEBUG
else
  CPPFLAGS += -DNDEBUG
endif

ifeq ($(ENABLE_EXPENSIVE_CHECKS),1)
  CPPFLAGS += -D_GLIBCXX_DEBUG -DXDEBUG
endif

LIBFFI_INCDIR = /usr/include

CMAKE_OPTIONS = -DCMAKE_C_COMPILER:FILEPATH=$(CC)
CMAKE_OPTIONS += -DCMAKE_CXX_COMPILER:FILEPATH=$(CXX)
CMAKE_OPTIONS += -DCMAKE_C_CFLAGS:STRING="$(CPPFLAGS) $(CFLAGS)"
CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS:STRING="$(CPPFLAGS) $(CXXFLAGS)"
CMAKE_OPTIONS += -DCMAKE_EXE_LINKER_FLAGS:STRING="$(LDFLAGS) -fPIE"
CMAKE_OPTIONS += -DCMAKE_SHARED_LINKER_FLAGS:STRING="$(LDFLAGS) -fPIC"
CMAKE_OPTIONS += -DCMAKE_INSTALL_PREFIX:FILEPATH=$(CMAKE_PREFIX)
CMAKE_OPTIONS += -DCMAKE_AR:FILEPATH=/usr/bin/ar

CMAKE_OPTIONS += -DCMAKE_INSTALL_BINDIR:STRING="$(CMAKE_INSTALL_BINDIR)"
CMAKE_OPTIONS += -DCMAKE_INSTALL_LIBDIR:STRING="$(CMAKE_INSTALL_LIBDIR)"
CMAKE_OPTIONS += -DCMAKE_INSTALL_LIBEXECDIR:STRING="$(CMAKE_INSTALL_LIBEXECDIR)"
CMAKE_OPTIONS += -DCMAKE_INSTALL_INCLUDEDIR:STRING="$(CMAKE_INSTALL_INCLUDEDIR)"
CMAKE_OPTIONS += -DCMAKE_INSTALL_DATADIR:STRING="$(CMAKE_INSTALL_DATADIR)"
CMAKE_OPTIONS += -DCMAKE_INSTALL_DATAROOTDIR:STRING="$(CMAKE_INSTALL_DATAROOTDIR)"
CMAKE_OPTIONS += -DCMAKE_SUPPRESS_REGENERATION:BOOL=ON
CMAKE_OPTIONS += -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

CMAKE_OPTIONS += -DCMAKE_MAKE_PROGRAM:STRING="/usr/bin/gmake"
CMAKE_OPTIONS += -DCMAKE_ASM_COMPILER:STRING="$(CC)"
CMAKE_OPTIONS += -DCMAKE_BUILD_RPATH="$(TOPBUILDDIR)/lib\;$(TOPBUILDDIR)/lib64"
CMAKE_OPTIONS += -DCMAKE_INSTALL_RPATH="$(CMAKE_INSTALL_LIBDIR)"

ifeq ($(DEBUG_BUILD),1)
  CMAKE_OPTIONS += -DCMAKE_BUILD_TYPE:STRING=Debug
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS_DEBUG:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS_DEBUG:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_ASM_FLAGS_DEBUG:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_ASM_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
else
  CMAKE_OPTIONS += -DCMAKE_BUILD_TYPE:STRING=Release
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS_DEBUG:STRING="$(OFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS_DEBUG:STRING="$(OFLAG)"
  CMAKE_OPTIONS += -DCMAKE_ASM_FLAGS_DEBUG:STRING="$(OFLAG)"
  CMAKE_OPTIONS += -DCMAKE_ASM_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_C_FLAGS_RELEASE:STRING="$(OFLAG) $(GFLAG)"
  CMAKE_OPTIONS += -DCMAKE_CXX_FLAGS_RELEASE:STRING="$(OFLAG) $(GFLAG)"
endif

CMAKE_OPTIONS += -DFFI_INCLUDE_DIR:STRING=$(LIBFFI_INCDIR)
CMAKE_OPTIONS += -DFFI_LIBRARY_DIR:STRING=$(LIBFFI_LIBDIR)

ifeq ($(HAVE_SPHINX),1)
  CMAKE_OPTIONS += -DSPHINX_EXECUTABLE:STRING="/usr/bin/sphinx-build"
  CMAKE_OPTIONS += -DSPHINX_OUTPUT_HTML:BOOL=ON
  CMAKE_OPTIONS += -DSPHINX_OUTPUT_MAN:BOOL=ON
  CMAKE_OPTIONS += -DSPHINX_WARNINGS_AS_ERRORS:BOOL=OFF
endif

# Benchmarks
CMAKE_OPTIONS += -DHAVE_POSIX_REGEX=true
CMAKE_OPTIONS += -DHAVE_STEADY_CLOCK=true

ifeq ($(ENABLE_CUDA),1)
CMAKE_OPTIONS += -DENABLE_CUDA:BOOL=ON
endif

CMAKE_OPTIONS += -DENABLE_PTHREAD:BOOL=ON
CMAKE_OPTIONS += -DENABLE_QUNIT_CPU_PARALLEL:BOOL=ON
CMAKE_OPTIONS += -DENABLE_EMIT_LLVM:BOOL=OFF
CMAKE_OPTIONS += -DENABLE_ROT_API:BOOL=ON
CMAKE_OPTIONS += -DENABLE_REG_GATES:BOOL=ON
CMAKE_OPTIONS += -DENABLE_VM6502Q_DEBUG:BOOL=ON
CMAKE_OPTIONS += -DENABLE_BCD:BOOL=ON
CMAKE_OPTIONS += -DCPP_STD=17
CMAKE_OPTIONS += -DQRACK_CPP_STD_OPT:STRING=-std=c++17

# Python
CMAKE_OPTIONS += -DPYTHON:FILEPATH="$(PYTHON)"
CMAKE_OPTIONS += -DPYTHON_EXECUTABLE:FILEPATH="$(PYTHON)"
CMAKE_OPTIONS += -DPYTHON_LIBRARY:FILEPATH="$(PYTHON_LIBRARY)"
CMAKE_OPTIONS += -DPYTHON_LIBRARIES:FILEPATH="$(PYTHON_LIBRARIES)"
CMAKE_OPTIONS += -DPYTHON_INCLUDE_DIRS:FILEPATH="$(PYTHON_INCLUDE_DIRS)"
CMAKE_OPTIONS += -DSWIG_EXECUTABLE:FILEPATH="$(SWIG_EXECUTABLE)"

# Path
PATH :=

ifeq ($(ENABLE_CUDA),1)
PATH = $(CUDA_INSTALL)/bin:/usr/bin:/usr/sbin:/bin:/sbin:/usr/local/bin
else
PATH = /usr/bin:/usr/sbin:/bin:/sbin:/usr/local/bin
endif

LLVM_BUILD_ENV = CC="$(CC)"
LLVM_BUILD_ENV += CXX="$(CXX)"
LLVM_BUILD_ENV += CFLAGS="$(CFLAGS)"
LLVM_BUILD_ENV += CXXFLAGS="$(CXXFLAGS)"
LLVM_BUILD_ENV += LDFLAGS="$(LDFLAGS)"
LLVM_BUILD_ENV += LD_OPTIONS="$(LD_OPTIONS)"
LLVM_BUILD_ENV += PATH="$(PATH)"
LLVM_BUILD_ENV += PYTHON="$(PYTHON)"
LLVM_BUILD_ENV += PYTHON_LIBRARY="$(PYTHON_LIBRARY)"
LLVM_BUILD_ENV += LANG="C"
LLVM_BUILD_ENV += LC_ALL="C"

configure:
	( cd $(TOPSRCDIR) ; \
	  mkdir -p $(TOPBUILDDIR) ; \
	  if [ ! -f $(TOPBUILDDIR)/.configured ] ; \
	    then echo "Configuring $(QRACK_TARGET) QRACK with $(CMAKE) $(CMAKE_OPTIONS)" ; \
	    cd $(TOPBUILDDIR) ; \
	    /usr/bin/env - $(LLVM_BUILD_ENV) $(CMAKE) $(CMAKE_OPTIONS) $(TOPSRCDIR) ; \
	    if [ $$? -eq 0 ] ; \
	    then \
	    	cd $(TOPBUILDDIR) ; \
	    	touch $(TOPBUILDDIR)/.configured ; \
	    else \
		echo "Failure configuring $(QRACK_TARGET) LLVM." ; \
		exit 1; \
	     fi \
	  else \
	    	echo "LLVM has already been configured." ; \
	  fi )

confclean:
	  rm -f $(TOPBUILDDIR)/.configured

build: configure
	( if [ ! -d $(TOPBUIDDIR) ] || [ ! -f $(TOPBUILDDIR)/.configured ] ; \
	    then echo "LLVM has not been configured" ; \
	  else \
	    cd $(TOPBUILDDIR) ; \
	    $(GMAKE) $(GMAKE_NUMJOBS) ; \
	    touch $(TOPBUILDDIR)/.built ; \
	  fi )

clean:
	( if [ -f $(TOPBUIDLDIR) ] ; \
	    then cd $(TOPBUILDDIR) ; \
	    $(GMAKE) clean ; \
	    rm -f $(TOPBUILDDIR)/.built ; \
	  fi )

clobber:
	  rm -rf $(TOPBUILDDIR)

all: build


