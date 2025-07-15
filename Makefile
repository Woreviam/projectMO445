FLAGS = -fPIC -std=gnu11 -Wall -Wno-unused-result -pedantic

BIN = $(NEWIFT_DIR)/bin

LIBSVM_DIR = $(NEWIFT_DIR)/externals/libsvm
LIBNIFTI_DIR= $(NEWIFT_DIR)/externals/libnifti
LIBJPEG_DIR= $(NEWIFT_DIR)/externals/libjpeg
LIBPNG_DIR= $(NEWIFT_DIR)/externals/libpng
TSNE_DIR= $(NEWIFT_DIR)/externals/tsne
ZLIB=$(NEWIFT_DIR)/externals/zlib
#LIBSDL2_DIR=$(NEWIFT_DIR)/externals/libsdl2_minimal

LIBIFT_INC   = -I $(NEWIFT_DIR)/include
LIBIFT_LD	= -L $(NEWIFT_DIR)/lib -lift
LIBSVM_INC   = -I $(LIBSVM_DIR)/include 
LIBSVM_LD	= -L $(LIBSVM_DIR)/lib -lstdc++
LIBCBLAS_INC = -I /usr/local/opt/openblas/include
LIBCBLAS_LD  = -L /usr/local/opt/openblas/lib -L /lib64/atlas-sse3 -L /usr/lib/atlas-base -llapack -lblas -lcblas

LIBNIFTI_INC = -I $(LIBNIFTI_DIR)/include
LIBJPEG_INC = -I $(LIBJPEG_DIR)/include
TSNE_INC = -I $(TSNE_DIR)/include
ZLIB_INC = -I $(ZLIB)/include
#LIBSDL2_INC = -I $(LIBSDL2_DIR)/include
#LIBSDL2_LD = -L $(LIBSDL2_DIR) -lSDL2

EXTERNALS_LD = -lm

OS=$(shell uname -s)

# For Mac OS, add a special flag to enable the openmp. Install: brew install libomp
# Otherwise, add the openmp flag commonly 
ifeq ($(OS), Darwin)
	EXTERNALS_LD += -Xpreprocessor -fopenmp -lomp # -lomp is required to link the program, lib, and openmp
else
	EXTERNALS_LD += -fopenmp
endif

INCLUDES = $(LIBIFT_INC) $(LIBSVM_INC) $(LIBCBLAS_INC) $(LIBNIFTI_INC) $(LIBJPEG_INC) $(TSNE_INC) $(ZLIB_INC)
LIBS	 = $(LIBIFT_LD) $(LIBSVM_LD) $(LIBCBLAS_LD) $(EXTERNALS_LD)



# PYTHON_LIB = $(shell python-config 2>/dev/null)

# including the Python Dev Library
# ifdef PYTHON_LIB
#	 INCLUDES += $(shell python-config --includes)
#	LIBS += $(shell python-config --libs)
# endif


ifeq ($(IFT_DEBUG), 1)
	export FLAGS += -pg -g -fsanitize=address -fsanitize=leak -DIFT_DEBUG=1
else
	export FLAGS += -O3
endif


#FOR GPU support uncomment the next line (in all the Makefiles that will use GPU)
#IFT_GPU=1

#CUDA path in case IFT_GPU is enabled
export CUDA_DIR1=/usr/local/cuda
export CUDA_DIR2=/opt/cuda

ifeq ($(IFT_GPU), 1)
	export FLAGS += -DIFT_GPU=1
	INCLUDES += -I $(CUDA_DIR1)/include
	INCLUDES += -I $(CUDA_DIR2)/include
	LIBS	 += -L $(CUDA_DIR1)/lib64 -L $(CUDA_DIR2)/lib64 -lcublas -lcudart
    MAKE_CMD=$(MAKE) IFT_GPU=1
else
    MAKE_CMD=$(MAKE)
endif


libift:
	cd $(NEWIFT_DIR); $(MAKE_CMD); cd -;


$@.c: $@.c
	$(MAKE) libift
	$(CC) $(FLAGS) $@.c -o $(BIN)/$@ $(INCLUDES) $(LIBS)


clean:
	rm -rf $(BIN)/*;
