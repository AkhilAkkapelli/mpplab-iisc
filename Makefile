# Fortran compiler
FORT_COMP = gfortran

# Source and executable
SRC = MultiPrecisionInteger.f90 TestHighPrecisionInteger.f90
ifeq ($(OS), Windows_NT)
	EXE := mpl.exe
	FFLAGS := -D_WIN64
else
	EXE := mpl
	FFLAGS := -D__linux__
endif

# Compiler Flags
# CFLAGS = -i8 -fp-model precise -O3 -mcmodel=large -qopenmp -qmkl-ilp64  -I/opt/intel/oneapi/mkl/2024.2/include/ -I/opt/intel/oneapi/mkl/2024.2/include/intel64/ilp64/

# Libraries
# LIB = /opt/intel/oneapi/mkl/2024.2/lib/libmkl_blas95_ilp64.a /opt/intel/oneapi/mkl/2024.2/lib/libmkl_lapack95_ilp64.a

# Build target
all: $(EXE)

$(EXE): $(SRC)
	@$(FORT_COMP) $(FFLAGS) -o $(EXE) $(SRC) -cpp

clean:
	@rm $(EXE)


