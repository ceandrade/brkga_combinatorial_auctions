################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on : Jul 09, 2011 by andrade
# Last update: Jan 20, 2014 by andrade
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################

###############################################################################
# User Defines
###############################################################################

# Set debug mode
#USER_DEFINES += -DDEBUG -DFULLDEBUG

###############################################################################
# User Flags
###############################################################################

ifeq ($(shell uname -s),Darwin)
    OS = macosx
else
    OS = linux
endif

## Compiler flags for debugging
#USER_FLAGS = -g3 -fexceptions  -fno-omit-frame-pointer \
	#-fno-optimize-sibling-calls -fno-inline

# Compiler flags for performace
USER_FLAGS += -Wall -Wextra -O3 -funroll-loops -ftracer \
	-fpeel-loops -fprefetch-loop-arrays -fomit-frame-pointer -fopenmp

###############################################################################
# Build options
###############################################################################

EXE = brkga-ca

###############################################################################
# The user source files
###############################################################################

# Include dir
USER_INCDIRS = \
	./include \
	./brkga \
	./lib

# Source directories
SRC_DIRS = \
	./problem_defs/ \
	./decoders/ \
	./cplex/ \
	./	

# Object files
OBJS = \
	./problem_defs/combinatorial_auction.o \
	./decoders/combinatorial_auction_decoder.o \
	./cplex/branch_callback.o \
	./cplex/heuristic_callback.o \
	./cplex/chromosome_generator.o \
	./cplex/take_duals.o \
	./main.o
	#./main_tunning.o

###############################################################################
# Lib and include definitions
###############################################################################

##############################
# Cplex and Concert settings
##############################

LIBFORMAT = static_pic

ifeq ($(OS), macosx)
	SYSTEM = x86-64_osx
	ILOG_DIR = ~/ILOG/CPLEX_Studio1262
else
	SYSTEM = x86-64_linux
	ILOG_DIR = ~/ILOG/CPLEX_Studio1262
endif

CPLEXDIR = $(ILOG_DIR)/cplex
CONCERTDIR = $(ILOG_DIR)/concert

CPLEXINCDIRS = \
	$(CPLEXDIR)/include \
	$(CONCERTDIR)/include

CPLEXLIBDIRS = \
	$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \
	$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CPLEXLIBS = -lilocplex -lcplex -lconcert -lm
#
# For Mac OSX
ifeq ($(OS), macosx)
	CPLEXLIBS += -framework CoreFoundation -framework IOKit
endif

# x86-64
CPLEXFLAGS = -fPIC -pthread -fexceptions -DIL_STD -DNDEBUG

##############################
# Consolidate paths
##############################

# Consolidate include paths
USER_INCDIRS += \
	$(CPLEXINCDIRS)
	
# Consolidate include paths
USER_LIBDIRS += \
	$(CPLEXLIBDIRS)

# Libraries necessary to link.
LIBS = \
	$(CPLEXLIBS)
	
# Compiler flags
USER_FLAGS += $(CPLEXFLAGS)
	

###############################################################################
# Compiler defs
###############################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = $(USER_FLAGS)

# Necessary Include dirs
# Put -I in front of dirs
INCLUDES = `for i in $(USER_INCDIRS); do echo $$i | sed -e 's/^/-I/'; done`

# Necessary library dirs
# Put -L in front of dirs
LIBDIRS = `for i in $(USER_LIBDIRS); do echo $$i | sed -e 's/^/-L/'; done`

###############################################################################
# Build Rules
###############################################################################

all: $(EXE)

.PHONY: all strip clean doc docclean depclean
.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	@echo "--> Linking objects... "
	$(CXX) $(CXXFLAGS) $(LINKER_FLAGS) $(OBJS) -o $@ $(LIBDIRS) $(LIBS) 
	@echo

.cpp.o:
	@echo "--> Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(USER_DEFINES) -c $< -o $@
	@echo

strip: $(EXE)
	@echo "-->Stripping $<..."
	strip $<
	@echo

doc:
	doxygen Doxyfile

clean:
	@echo "--> Cleaning compiled..."
	rm -rf $(EXE) $(OBJS) $(CLIQUER_OBJS)
	rm -rf brkga-ca.dSYM
	
depclean: docclean clean
	rm -rf `for i in $(SRC_DIRS); do echo $$i*~; done` 
	rm -rf Debug
	
docclean:
	@echo "--> Cleaning doc..."
	rm -rf doc

