CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../struct
OPTM =  -g


VPATH = ../01_data_structures:../02_parsers:../03_output:\
	../04_geometric_match:../05_preproc:../06_main:\
	../07_utils


PLATFORM = $(shell uname)


LIBS =  -lm -lc  -L/usr/lib  -llapack  
INCLUDE = -I ../00_include


DEFINES =   # -DDMALLOC # -DDEBUG
SOURCES.F = 
SOURCES.f = 

SOURCES.c =  struct.c struct_aa_code.c struct_compare_reduced_reps.c struct_align_backbone.c  \
	     struct_complement_match.c  struct_descr_out.c struct_digest.c  struct_geometry.c \
	     struct_linalg.c  struct_map.c struct_map2rotation.c       \
	     struct_pdb_out.c  struct_quaternion.c struct_read_pdb.c   \
	     struct_read_cmd_file.c struct_representation.c            \
	     struct_recursive_output.c struct_scoring_fn.c  struct_sse2descr.c \
	     struct_structure2sse.c \
	     struct_table.c  struct_table_input.c struct_input.c struct_utils.c



AUX =  makefile   README \
       struct_comparison.h  struct_geometry.h  struct.h  struct_pdb.h  \
       struct_utils.h


SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)



.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm  $(OBJECTS) 

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


