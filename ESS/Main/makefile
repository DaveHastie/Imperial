# This file is part of ESS++.
#      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
#                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
#                    David Hastie (d.hastie@imperial.ac.uk)  
#      2010
#      
# ESS++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESS++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ESS++.  If not, see <http://www.gnu.org/licenses/>.
#

CUDA			= 1
HPC				= 1

LIBDIR			= ../Routines/
LIBDIR2			= ../Classes/

ifeq ($(HPC),0)
	GOTODIR		= /home/dhastie/LinAlg/GotoBLAS2
endif

ifeq ($(HPC),0)
	CC	 		= gcc
	CXX 		= g++
else
	CC			= icc	
	CXX			= icpc
endif

IFLAGSTMP		= -I$(LIBDIR)
ifeq ($(HPC),0)
	LFLAGSTMP		= -L$(LIBDIR) -L$(GOTODIR)
else
	LFLAGSTMP		= -L$(LIBDIR)
endif
LIBSTMP			= -lgsl

ifeq ($(CUDA),0)
ifeq ($(HPC),0)
	IFLAGS		= $(IFLAGSTMP)
	LFLAGS 		= $(LFLAGSTMP)
	LIBS 		= $(LIBSTMP) -lgslcblas
else
	IFLAGS		= $(IFLAGSTMP) -I$(GSL_HOME)/include/
	LFLAGS		= $(LFLAGSTMP) -L$(GSL_HOME)/lib/ -L$(MKL_HOME)/lib/64/
	LIBS		= $(LIBSTMP) -lmkl_lapack -lmkl -lguide -lpthread
endif
else
ifeq ($(HPC),0)
	IFLAGS		= $(IFLAGSTMP) -I/usr/local/cuda/include -I/usr/local/cula/include
	LFLAGS		= $(LFLAGSTMP) -L/usr/local/cuda/lib64 -L/usr/local/cula/lib64
	LIBS		= $(LIBSTMP) -lgoto2 -lpthread -llapack -lcuda -lcublas -lcudart -lcula
else
	IFLAGS		= $(IFLAGSTMP) -I$(GSL_HOME)/include/ -I$(CUDA_HOME)/include/ -I$(CULA_ROOT)/include/ 
	LFLAGS		= $(LFLAGSTMP) -L$(GSL_HOME)/lib/ -L$(MKL_HOME)/lib/64/ -L$(CUDA_HOME)/lib64/ -L$(CULA_ROOT)/lib64/ 
	LIBS		= $(LIBSTMP) -lmkl_lapack -lmkl -lguide -lpthread -lcuda -lcublas -lcudart -lcula
endif
endif

CXXFLAGS		= 
ifeq ($(HPC),0)
	CPPFLAGS		= -W -Wall --pedantic -O3  $(IFLAGS) -D_CUDA_=$(CUDA)
	#CPPFLAGS		= -W -Wall -g -pg -O0 --pedantic $(IFLAGS) -D_CUDA_=$(CUDA)
	#CPPFLAGS		= -W -Wall -O4 -g $(IFLAGS) -D_CUDA_=$(CUDA)
	# Use the compile option -pg for profiling.
else
	CPPFLAGS		= -Wall -O3  $(IFLAGS) -D_CUDA_=$(CUDA)
endif


all:		ESS .depend 

clean:		
		rm -f *.o
		rm -f .depend
		rm -f ESS
		rm -f *~
		rm -f ../Routines/*.o
		rm -f ../Routines/*~
		rm -f ../Classes/*.o
		rm -f ../Classes/*~
		rm -f ../.depend
dep .depend:
	$(CXX) $(CPPFLAGS) -MM ESS.cc >>.depend

infiles1	= ESS.o  $(LIBDIR)dyn_name.o $(LIBDIR)matrix_handling.o  $(LIBDIR)post_processing.o $(LIBDIR)rand.o $(LIBDIR2)String_Matrices.o $(LIBDIR2)Int_Matrices.o $(LIBDIR2)Double_Matrices.o  $(LIBDIR2)Double_Matrices_cont.o  $(LIBDIR)regression.o $(LIBDIR)cond_post.o $(LIBDIR)xml_file_read.o $(LIBDIR2)Prior_param.o $(LIBDIR)moves.o $(LIBDIR2)Temperatures.o $(LIBDIR2)g_AdMH.o $(LIBDIR2)Move_monitor.o $(LIBDIR2)DR.o $(LIBDIR2)CM.o $(LIBDIR2)Int_Matrices_var_dim.o

ESS:	$(infiles1)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o ESS $(infiles1) $(LIBS)

