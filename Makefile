CC      = gcc
MPI     = mpicc
CFLAGS += -Wall -Wextra -g -O4 -fopenmp
LDLIBS += -lm -lrt
OFLAGS += -O3 -mavx
MPI_INCLUDE_DIR=-I/usr/include/openmpi-x86_64

BIN_DIR=build

.PHONY: all
all: $(BIN_DIR)/stencil-basic $(BIN_DIR)/stencil-omp $(BIN_DIR)/stencil-mpi

$(BIN_DIR)/stencil-basic: src/stencil_basic.c
        $(CC) $(CFLAGS) $(OFLAGS) $(LDLIBS) $< -o $@

$(BIN_DIR)/stencil-omp: src/stencil_omp.c
        $(CC) $(CFLAGS) $(OFLAGS) $(LDLIBS) $< -o $@

$(BIN_DIR)/stencil-mpi: src/stencil_mpi.c
        $(MPI) $(CFLAGS) $(OFLAGS) $(LDLIBS) $(MPI_INCLUDE_DIR) -o $@ $<

clean:
        rm -f $(BIN_DIR)/stencil-basic $(BIN_DIR)/stencil-omp $(BIN_DIR)/stencil-mpi $(BIN_DIR)/stencil-all


