CC = mpicc

INC = ./include
SRC = ./src
OBJ = ./obj
BIN = ./bin
EXE = main

SRCS = $(wildcard $(SRC)/*.c)
OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))

GMRES_DIR = $(PETSC_DIR)/src/ksp/ksp/impls/gmres/

CFLAGS = -g -Wall -I$(INC) -I$(PETSC_DIR)/include -I$(PETSC_DIR)$(PETSC_ARCH)/include -I$(GMRES_DIR) -I$(PCG_DIR)/include -I$(METIS_DIR)/include
LIBS_PATH = -L$(PETSC_DIR)$(PETSC_ARCH)/lib -L$(PCG_DIR)/src -L$(METIS_DIR)/lib
LIBS = -lm -lpetsc -llapack -lpcg_random -lmetis

all: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$(EXE) $(OBJS) $(LIBS_PATH) $(LIBS)

$(OBJS): $(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -c $? -o $@

.PHONY: clean

clean: all $(OBJS)
	rm -f $(OBJS)
