###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda

##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=g++
CC_FLAGS=
CC_LIBS=

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc
NVCC_FLAGS= -arch=sm_61 -O3
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

# Test file diretory:
TEST_DIR = extra

##########################################################

## Make variables ##

# Target executable name:
EXE = runHamming

# Object files:
OBJS = $(OBJ_DIR)/main.o # $(OBJ_DIR)/cuda_kernel.o

# Test files:
TEST1 = $(TEST_DIR)/test.fa $(TEST_DIR)/test.fa
TEST2 = $(TEST_DIR)/test1k.fa $(TEST_DIR)/test10k.fa
TEST3 = $(TEST_DIR)/test1k.fa $(TEST_DIR)/test100k.fa

# Result file
RESULT = results.tsv

# Log file
LOG = $(TEST_DIR)/run.log

NVPROF = nvprof
NVPROF_FLAGS = -f --print-gpu-trace --device-buffer-size on
NVPROF_OUT = $(TEST_DIR)/runHamming-analysis.nvprof
NVPROF_LOG = $(TEST_DIR)/nvprof.log

##########################################################

## Compile ##

# # Link c++ and CUDA compiled object files to target executable:
# $(EXE) : $(OBJS)
# 	$(CC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# # Compile main .cpp file to object files:
# $(OBJ_DIR)/%.o : %.cpp
# 	$(CC) $(CC_FLAGS) -c $< -o $@

# # Compile C++ source files to object files:
# $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp include/%.h
# 	$(CC) $(CC_FLAGS) -c $< -o $@

# # Compile CUDA source files to object files:
# $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu $(INC_DIR)/%.cuh
# 	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

# Make options

# Link c++ and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(NVCC) $(NVCC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : %.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

all: $(EXE)

.PHONY : all

# Clean objects in object directory.
clean:
	@rm -f $(OBJS) $(EXE)

# Clean objects and test outputs.
cleanall:
	@rm -f $(OBJS) $(EXE) $(RESULT) $(LOG) $(NVPROF_OUT) $(NVPROF_LOG)

.PHONY : clean cleanall

# Test options
testall: test test10k test100k

# Make test #1:
test: $(EXE) $(TEST1)
	@./$(EXE) $(TEST1) > $(LOG)
	@echo "Test PASSED" || echo "Test FAILED"

# Make test #2:
test10k: $(EXE) $(TEST2)
	@./$(EXE) $(TEST2) > $(LOG)
	@echo "Test 10k PASSED" || echo "Test 10k FAILED"

# Make test #3:
test100k: $(EXE) $(TEST3)
	@./$(EXE) $(TEST3) > $(LOG)
	@echo "Test 100k PASSED" || echo "Test 100k FAILED"

.PHONY : testall test test10k test100k

# Make test #2:
profile: $(EXE) $(TEST2)
	@$(NVPROF) $(NVPROF_FLAGS) -o $(NVPROF_OUT) ./$(EXE) $(TEST2) > $(NVPROF_LOG)
	@echo "Profiling COMPLETED" || echo "Profiling FAILED"

.PHONY : profile
