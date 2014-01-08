############################################################################
# Makefile for the src dir of libsr
# Oliver Hinds <oph@bu.edu> 2005-06-22
#######################################################################

## if any command to generate target flibvps, the target is deleted.
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

default: $(OBJ_FILES)

############################### CLEAN ################################

clean:
	@$(ECHO) 'make: removing object and autosave files'
	-$(RM) -f *.o *~

######################################################################
### $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/makefile,v $
### Local Variables:
### mode: makefile
### fill-column: 76
### comment-column: 0
### End: 
