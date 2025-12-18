PROJECT := connected_components_mpi

CC := mpicc
GCC := gcc

CFLAGS  := -Wall -Wextra -Wpedantic -std=c11 -O3 -Isrc -fopenmp
LDFLAGS := -fopenmp
LDLIBS  := -lm

# Directories
SRC_DIR := src
BIN_DIR := bin
OBJ_DIR := obj

# Source files
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

TARGET := $(BIN_DIR)/$(PROJECT)

# Pretty Output
ECHO := /bin/echo -e
COLOR_RESET := \033[0m
COLOR_GREEN := \033[1;32m
COLOR_YELLOW := \033[1;33m
COLOR_BLUE := \033[1;34m
COLOR_MAGENTA := \033[1;35m
COLOR_CYAN := \033[1;36m

.PHONY: all
all: $(TARGET)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(TARGET): $(OBJS) | $(BIN_DIR)
	@$(ECHO) "$(COLOR_GREEN)Linking (MPI+OpenMP):$(COLOR_RESET) $@"
	@$(CC) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)
	@$(ECHO) "$(COLOR_CYAN)Build complete!$(COLOR_RESET)"
	@$(ECHO) "$(COLOR_CYAN)Run with: srun ... $(TARGET) -n trials input.bin$(COLOR_RESET)"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@$(ECHO) "$(COLOR_BLUE)Compiling:$(COLOR_RESET) $<"
	@$(CC) $(CFLAGS) -c $< -o $@

.PHONY: converter
converter: | $(BIN_DIR)
	@$(ECHO) "$(COLOR_BLUE)Compiling:$(COLOR_RESET) mtx_to_bin"
	@$(GCC) -Wall -Wextra -O3 -fopenmp src/converter/mtx_to_bin.c -o $(BIN_DIR)/mtx_to_bin -lm
	@$(ECHO) "$(COLOR_GREEN)Converter built: $(BIN_DIR)/mtx_to_bin$(COLOR_RESET)"

.PHONY: clean
clean:
	@$(ECHO) "$(COLOR_YELLOW)Cleaning...$(COLOR_RESET)"
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@$(ECHO) "$(COLOR_GREEN)✓ Clean complete$(COLOR_RESET)"

.PHONY: rebuild
rebuild: clean all

.PHONY: help
help:
	@echo "Hybrid MPI+OpenMP Connected Components"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build MPI+OpenMP version (default)"
	@echo "  converter - Build mtx_to_bin utility"
	@echo "  clean     - Remove build artifacts"
	@echo "  rebuild   - Clean and rebuild"
	@echo "  help      - Show this message"
	@echo ""
	@echo "Usage:"
	@echo "  SLURM: sbatch run_slurm.sh"
	@echo ""
	@echo "Environment:"
	@echo "  OMP_NUM_THREADS - OpenMP threads per rank"
	@echo "  OMP_PROC_BIND   - Thread affinity"

