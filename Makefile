CC      = gcc
CFLAGS += -Wall -Wextra -g -O4 -fopenmp
LDLIBS += -lm -lrt

BIN_DIR=build

.PHONY: all
all: $(BIN_DIR)/stencil-basic $(BIN_DIR)/stencil-omp

$(BIN_DIR)/stencil-basic: src/stencil_basic.c
	$(CC) $(CFLAGS) $(LDLIBS) $< -o $@

$(BIN_DIR)/stencil-omp: src/stencil_omp.c
	$(CC) $(CFLAGS) $(LDLIBS) $< -o $@

clean:
	rm -f stencil-baseline stencil-openmp
