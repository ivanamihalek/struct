CCOMPL =  gcc 
CFLAGS =  -Wall  -O
#CFLAGS += -fopenmp
NVCC = nvcc

NAME = ../struct


VPATH = ../01_data_structures:../02_parsers:../03_output:\
	../04_geometric_match:../05_preproc:../06_main:\
	../07_utils

CUDA_FLAGS = -I/opt/cuda-sdk/C/common/inc/   -arch sm_20

# need  -Wl,-no_pie
LIBS    =  -lm -lc  -L/usr/lib  -llapack   #-lprofiler  -Wl,-no_pie
INCLUDE =  -I ../00_include


DEFINES =   # -DGPU  # -DOMP  -DDEBUG  -DDMALLOC 

SOURCES.c =  struct.c struct_aa_code.c  struct_alignment_out.c  struct_curvature.c	         \
	     struct_map_reduced_reps.c  struct_align_backbone.c  struct_hungarian.c              \
	     struct_direction_match.c  struct_descr_out.c struct_digest.c  struct_geometry.c    \
	     struct_input_structure.c struct_linalg.c  struct_map.c struct_map2rotation.c        \
	     struct_optimize_bb_almt.c  struct_opt_quat.c                                        \
	     struct_pdb_out.c  struct_process_input.c  struct_quaternion.c struct_read_pdb.c     \
	     struct_read_cmd_file.c struct_representation.c struct_results_out.c                 \
	     struct_map_out.c struct_scoring_fn.c  struct_sse2descr.c                            \
	     struct_structure2sse.c struct_table.c    struct_triplets.c                          \
	     struct_table_input.c  struct_utils.c struct_lmmin.c 	 

SOURCES.cu = 


AUX =  makefile   README \
       struct_binheap.h  struct_curve_tors.h  struct_geometry.h  struct.h  struct_hungarian.h	\
       struct_lmmin.h	struct_pdb.h  struct_structure2sse.h  struct_triplets.h  struct_utils.h


SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS = $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)

GPUOBJS = $(SOURCES.cu:.cu=.o)


#$(NAME): $(OBJECTS) $(GPUOBJS)
#	$(NVCC) -o $@  $(OBJECTS)  $(GPUOBJS) $(LIBS) 

$(NAME): $(OBJECTS) 
	$(CCOMPL) -o $@  $(OBJECTS)  $(LIBS) 

#.cu.o:  $(SOURCES.cu) doesn't work for some reason
#	$(NVCC) $(CUDA_FLAGS) -c $<

struct_exhaustive_cuda_pieces.o: struct_exhaustive_cuda_pieces.cu
	$(NVCC) $(CUDA_FLAGS) -c $<  $(INCLUDE)
.c.o:
	$(CCOMPL) $(CFLAGS) $(DEFINES)  -c $<  $(INCLUDE)

clean:
	rm  $(OBJECTS)  $(GPUOBJS)



