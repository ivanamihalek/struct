CCOMPL =  gcc 
CFLAGS =  -Wall  -O
#CFLAGS += -fopenmp


NAME = ../struct


VPATH = ../01_data_structures:../02_parsers:../03_output:\
	../04_geometric_match:../05_preproc:../06_main:\
	../07_utils


LIBS    =  -lm -lc  -L/usr/lib  -llapack #  -lprofiler
INCLUDE =  -I ../00_include


DEFINES =   # -DOMP  -DDEBUG  -DDMALLOC 

SOURCES.c =  struct.c struct_aa_code.c  struct_alignment_out.c                                   \
	     struct_map_reduced_reps.c struct_align_backbone.c  struct_hungarian.c               \
	     struct_complement_match.c  struct_descr_out.c struct_digest.c  struct_geometry.c    \
	     struct_input_structure.c struct_linalg.c  struct_map.c struct_map2rotation.c        \
	     struct_optimize_bb_almt.c                                                           \
	     struct_pdb_out.c  struct_process_input.c  struct_quaternion.c struct_read_pdb.c     \
	     struct_read_cmd_file.c struct_representation.c struct_results_out.c                 \
	     struct_map_out.c struct_scoring_fn.c  struct_sse2descr.c                            \
	     struct_structure2sse.c struct_table.c                                               \
	     struct_table_input.c  struct_utils.c 



AUX =  makefile   README \
       struct_comparison.h  struct_geometry.h  struct.h  struct_pdb.h  \
       struct_utils.h


SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(CFLAGS)  $(DEFINES) -o $@  $(OBJECTS)  $(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL) $(CFLAGS) $(DEFINES)  -c $<  $(INCLUDE)

clean:
	rm  $(OBJECTS) 



