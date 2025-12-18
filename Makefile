# Fortran compiler
FORT_COMP = gfortran

# Source and executable
SRC = MultiPrecisionInteger.f90 TestHighPrecisionInteger.f90
EXE = mpl

# Compiler Flags
# CFLAGS = -i8 -fp-model precise -O3 -mcmodel=large -qopenmp -qmkl-ilp64  -I/opt/intel/oneapi/mkl/2024.2/include/ -I/opt/intel/oneapi/mkl/2024.2/include/intel64/ilp64/

# Libraries
# LIB = /opt/intel/oneapi/mkl/2024.2/lib/libmkl_blas95_ilp64.a /opt/intel/oneapi/mkl/2024.2/lib/libmkl_lapack95_ilp64.a

# Build target
all: $(EXE)

$(EXE): $(SRC)
	$(FORT_COMP) -o $(EXE) $(SRC) -cpp

clean:
	rm -f $(EXE)

