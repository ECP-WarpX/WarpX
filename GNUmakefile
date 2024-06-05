AMREX_HOME  ?= ../amrex
PICSAR_HOME ?= ../picsar
OPENBC_HOME ?= ../openbc_poisson

#DEBUG = FALSE
DEBUG	= TRUE

WARN_ALL = TRUE
#WARN_ERROR=TRUE

#DIM = 1
#DIM = 2
DIM = 3

QED	       = TRUE
#QED_TABLE_GEN = TRUE

COMP = gcc
#COMP = intel
#COMP = pgi

TINY_PROFILE   = TRUE
#PROFILE       = TRUE
#COMM_PROFILE  = TRUE
#TRACE_PROFILE = TRUE

USE_OMP   = TRUE
USE_GPU   = FALSE

EBASE     = main

USE_PYTHON_MAIN = FALSE

USE_SENSEI_INSITU = FALSE
USE_ASCENT_INSITU = FALSE
USE_OPENPMD = TRUE

WarpxBinDir = Bin

USE_FFT = FALSE
USE_RZ = FALSE

USE_EB = FALSE

WARPX_HOME := .
include $(WARPX_HOME)/Source/Make.WarpX





VPATH_LOCATIONS += $(HEFFTE_HOME)/include
INCLUDE_LOCATIONS += $(HEFFTE_HOME)/include
LIBRARY_LOCATIONS += $(HEFFTE_HOME)/lib

libraries += -lheffte

ifeq ($(USE_CUDA),TRUE)
  libraries += -lcufft
else ifeq ($(USE_HIP),TRUE)
  # Use rocFFT.  ROC_PATH is defined in amrex
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocfft/include
  LIBRARY_LOCATIONS += $(ROC_PATH)/rocfft/lib
  LIBRARIES += -L$(ROC_PATH)/rocfft/lib -lrocfft
else
  libraries += -lfftw3_mpi -lfftw3f -lfftw3
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
