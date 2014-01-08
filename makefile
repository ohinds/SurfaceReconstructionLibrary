############################################################################
# Makefile for the CVCNS SurfaceReconstruction Library
# http://eslab.bu.edu
# Copyright 2006 Oliver Hinds <oph@bu.edu> 
#######################################################################

# project name
export PROJECT = libsr
export MAJOR_VER = 1
export MINOR_VER = 1
export RELEASE_VER = 0

# operating system, mac, linux and win are supported
export OS = linux

export BITS = 64

# whether to build a static or shared library
export LIB_TYPE = shared

# whether to compile with debug, optimize flags
export DEBUG = 0
export OPTIM = 1
export STRIP = 1
export PROF = 0
export MEMLEAK = 0

# directories 
export LIB_DEST_DIR = /usr/local/lib/
export HDR_DEST_DIR = /usr/local/include/
export LIB_DIR = /tmp/
export SRC_DIR = $(PWD)/src/
export OBJ_DIR = $(PWD)/obj/
export BIN_DIR = $(PWD)/bin/

################################ APPS ################################

export RM = /bin/rm -v
export ECHO = /bin/echo
export CC = /usr/bin/gcc
export AR = /usr/bin/ar
export INSTALL=sudo /usr/bin/install
export ROOTLN = sudo ln
export LDCONFIG=/sbin/ldconfig

# debug flag
ifeq ($(DEBUG),1)
	DEBUG_FLAG = -g -DSR_DEBUG_FLAG
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

# flags for the compiler and linker
VP_LIB_TYPE = shared
ifeq ($(OS),win)
	VP_LIB_TYPE = static
endif
ifeq ($(OS),mac)
	VP_LIB_TYPE = static
endif
ifeq ($(VP_LIB_TYPE),static)
	LIBVP_LIB_DIR = /usr/local/lib
	LIBVP = $(LIBVP_LIB_DIR)/libvp.a
endif

ifeq ($(OS),linux)
	export fPIC = -fPIC
endif


export CINCL =  -I$(SRC_DIR)
export CFLAGS = $(fPIC) -Wall -Werror $(MEMLEAK_FLAG) $(PROF_FLAG) $(DEBUG_FLAG) $(OPTIM_FLAG) $(STRIP_FLAG) $(CINCL)
export LDFLAGS = $(PROF_FLAG) $(LIBVP)

# bits of platform
ifeq ($(BITS),64)
	CFLAGS += -DBITS64
else
	CFLAGS += -DBITS32
endif

# differences between mac and *nix
ifeq ($(OS),mac)
	LIB_TYPE = static
	CFLAGS += -DMAC -I$(FINK_ROOT)/sw/include
	LDFLAGS += -lobjc -L$(FINK_ROOT)/sw/lib -ljpeg
endif

# specific flags for windows, these override earlier defs of cflags and ldflags
ifeq ($(OS),win)
	LIB_TYPE = static
        export CFLAGS = -Wall -Werror -I/usr/include $(PROF_FLAG) $(DEBUG_FLAG) $(OPTIM_FLAG) $(STRIP_FLAG)
	export LDFLAGS = -lm -lgsl $(LIBVP)
endif

############################## SUFFIXES ##############################

## if any command to generate target fails, the target is deleted.
# (see http://www.gnu.org/manual/make/html_chapter/make_5.html#SEC48)
.DELETE_ON_ERROR:

.SUFFIXES:
.SUFFIXES:  .o .c

# suffix rule for subsidiary source files
# (see http://www.gnu.org/manual/make/html_chapter/make_10.html#SEC111)
$(OBJ_DIR)/%.o: %.c %.h 
	@$(ECHO) '[make: building $@]'
	$(CC) $(CFLAGS) -o $@ -c $< 

HDR_FILES = $(wildcard *.h)
SRC_FILES = $(wildcard ./*.c)
TMP_FILES = $(patsubst ./%,$(OBJ_DIR)/%,$(SRC_FILES))
OBJ_FILES = $(TMP_FILES:.c=.o) 

default: $(PROJECT)
debug:	 
	$(MAKE) DEBUG=1 OPTIM=0 STRIP=0 $(PROJECT)
profile:	 
	$(MAKE) DEBUG=1 PROF=1 OPTIM=0 STRIP=0 $(PROJECT)
$(PROJECT): $(OBJ_FILES)
	@$(ECHO) 'make: building $@ for $(OS)...'
	cd $(SRC_DIR) && $(MAKE)
ifeq ($(LIB_TYPE),shared)
	$(CC)  -shared  -Wl,-soname,$(PROJECT).so.$(MAJOR_VER) $(OBJ_DIR)/*.o -o $(LIB_DIR)/$(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LDFLAGS) 
	@$(ECHO) '############################################'
	@$(ECHO) 'make: built [$@.so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER)] successfully!'
	@$(ECHO) '############################################'
else  # static
	$(AR) rcs $(LIB_DIR)/$(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(OBJ_DIR)/*.o
	@$(ECHO) '############################################'
	@$(ECHO) 'make: built [$@.a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER)] successfully!'
	@$(ECHO) '############################################'
endif


############################### INSTALL ################################

install: $(PROJECT)
	$(INSTALL) $(SRC_DIR)/*.h $(HDR_DEST_DIR)
	$(INSTALL) $(SRC_DIR)/*.extern $(HDR_DEST_DIR)
ifeq ($(LIB_TYPE),shared)
	@$(ECHO) 'make: installing $(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) for $(OS)...'
	$(INSTALL) $(LIB_DIR)/$(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)
	$(ROOTLN) -sf $(LIB_DEST_DIR)/$(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)/$(PROJECT).so.$(MAJOR_VER)
	$(ROOTLN) -sf $(LIB_DEST_DIR)/$(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)/$(PROJECT).so
ifeq ($(OS),linux)
	$(LDCONFIG) -n $(LIB_DEST_DIR)
endif
	@$(ECHO) '############################################'
	@$(ECHO) 'make: installed [$(PROJECT).so.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER)] successfully!'
	@$(ECHO) '############################################'
else # static library target
	@$(ECHO) 'make: installing $(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) for $(OS)...'
	($INSTALL) $(LIB_DIR)/$(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)
	$(ROOTLN) -sf $(LIB_DEST_DIR)/$(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)/$(PROJECT).a.$(MAJOR_VER)
	$(ROOTLN) -sf $(LIB_DEST_DIR)/$(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER) $(LIB_DEST_DIR)/$(PROJECT).a
	@$(ECHO) '############################################'
	@$(ECHO) 'make: installed [$(PROJECT).a.$(MAJOR_VER).$(MINOR_VER).$(RELEASE_VER)] successfully!'
	@$(ECHO) '############################################'
endif

############################### CLEAN ################################

clean:
	@$(ECHO) 'make: removing object and autosave files'
	-$(RM) -f $(OBJ_DIR)/*.o *~ $(SRC_DIR)/*~


######################################################################
### $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/makefile,v $
### Local Variables:
### mode: makefile
### fill-column: 76
### comment-column: 0
### End: 
