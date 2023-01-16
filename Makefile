CC = mpicc

INC = ./include
SRC = ./src
OBJ = ./obj
BIN = ./bin
EXE = main

SRCS = $(wildcard $(SRC)/*.c)
OBJS = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))

GMRES_DIR = $(PETSC_DIR)/src/ksp/ksp/impls/gmres/
MATSBAIJ_DIR = $(PETSC_DIR)/src/mat/impls/sbaij/seq/

CFLAGS = -g -Wall -I$(INC) -I$(PETSC_DIR)/include -I$(PETSC_DIR)$(PETSC_ARCH)/include -I$(GMRES_DIR) -I$(MATSBAIJ_DIR) -I$(PCG_DIR)/include -I$(METIS_DIR)/include
CFLAGS_SERVER = $(CFLAGS) -I$(GK_DIR)/include

LIBS_PATH = -L$(PETSC_DIR)$(PETSC_ARCH)/lib -L$(PCG_DIR)/src -L$(METIS_DIR)/lib
LIBS_PATH_SERVER = $(LIBS_PATH) -L$(GK_DIR)/lib

LIBS = -lm -lpetsc -llapack -lpcg_random -lmetis
LIBS_SERVER = $(LIBS) -lGKlib

all: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$(EXE) $(OBJS) $(LIBS_PATH) $(LIBS)
	
server: CFLAGS=CFLAGS_SERVER

server: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$(EXE) $(OBJS) $(LIBS_PATH_SERVER) $(LIBS_SERVER)

$(OBJS): $(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -c $? -o $@

.PHONY: clean

clean: all $(OBJS)
	rm -f $(OBJS)
