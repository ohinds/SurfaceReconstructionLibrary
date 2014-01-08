############################################################################
# Makefile for the CVCNS reconFromSlices application
# Oliver Hinds <oph@bu.edu> 2004-02-04
# 
#######################################################################

################################ VARS ################################

# project name
export PROJECT = reconFromSlices

# whether to compile with debug, optimize flags
export DEBUG = 0
export OPTIM = 1
export STRIP = 1
export PROF = 0
export MEMLEAK = 0

export OS = linux

# directories 
export BIN_DIR = $(PWD)
export BIN_DEST_DIR = /usr/local/bin

################################ APPS ################################

export RM = /bin/rm -v
export ECHO = /bin/echo
export CC = /usr/bin/gcc

# debug flag
ifeq ($(DEBUG),1)
	DEBUG_FLAG = -g
endif

# profile flag
ifeq ($(PROF),1)
	PROF_FLAG = -pg
endif

# optimize flag
ifeq ($(OPTIM),1)
	OPTIM_FLAG = -O
endif

# strip flag
ifeq ($(STRIP),1)
	STRIP_FLAG = -s
endif

# memleak catch flag
ifeq ($(MEMLEAK),1)
	MEMLEAK_FLAG = -DMEMLEAK
endif

# params for static compilation of libs
ifeq ($(OS),win)
	SR_LIB_TYPE = static
endif
ifeq ($(OS),mac)
	SR_LIB_TYPE = static
endif
ifeq ($(SR_LIB_TYPE),static)
	LIBSR_LIB_DIR = /usr/local/lib
	LIBSR = $(LIBSR_LIB_DIR)/libsr.a
else
	LIBSR = -lsr
endif

ifeq ($(OS),win)
	VP_LIB_TYPE = static
endif
ifeq ($(OS),mac)
	VP_LIB_TYPE = static
endif
ifeq ($(VP_LIB_TYPE),static)
	LIBVP_LIB_DIR = /usr/local/lib
	LIBVP = $(LIBVP_LIB_DIR)/libvp.a
else
	LIBVP = -lvp
endif

# flags for the compiler and linker
export CFLAGS = -Werror -Wall $(MEMLEAK_FLAG) $(PROF_FLAG) $(DEBUG_FLAG) $(OPTIM_FLAG) $(STRIP_FLAG) 
export LDFLAGS = $(PROF_FLAG) $(LIBSR)  $(LIBVP) -lm -lgsl -lgslcblas -lz -lglut

SRC_FILES = $(wildcard *.c) $(wildcard *.h)

############################## TARGETS ###############################
all: $(PROJECT)

debug: 
	$(MAKE) DEBUG=1 OPTIM=0 STRIP=0 $(PROJECT)

$(PROJECT): $(SRC_FILES)
	gcc -c $(CFLAGS) $@.c
	gcc  -o $(BIN_DIR)/$@ $@.o $(LDFLAGS) $(CFLAGS) 

install: $(PROJECT)
	@$(ECHO) 'make: installing $(PROJECT) for $(OS)...'
	cp $(BIN_DIR)/$(PROJECT) $(BIN_DEST_DIR)
	@$(ECHO) '############################################'
	@$(ECHO) 'make: installed $(PROJECT) successfully!'
	@$(ECHO) '############################################'

############################### CLEAN ################################

clean:
	@$(ECHO) 'make: removing object and autosave files'
	$(RM) -f *.o *~ $(PROJECT)

######################################################################
### $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/reconFromSlices/makefile,v $
### Local Variables:
### mode: makefile
### fill-column: 76
### comment-column: 0
### End: 
