PROJECT := connected_components_mpi

export CC := gcc

CFLAGS := -Wall -Wextra -Wpedantic -std=c11 -O3 -march=native -Isrc/core -Isrc/algorithms -Isrc/utils -fopenmp
LDFLAGS := -fopenmp -lm

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
	@$(ECHO) "$(COLOR_GREEN)Linking:$(COLOR_RESET) $@"
	@$(CC) $(LDFLAGS) $(OBJS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@$(ECHO) "$(COLOR_BLUE)Compiling:$(COLOR_RESET) $<"
	@$(CC) $(CFLAGS) -c $< -o $@

converter: | $(BIN_DIR)
	@$(ECHO) "$(COLOR_BLUE)Compiling:$(COLOR_RESET) mtx_to_bin"
	@gcc $(CFLAGS) $(LDFLAGS) $(SRC_DIR)/converter/mtx_to_bin.c -o $(BIN_DIR)/mtx_to_bin


.PHONY: clean
clean:
	@$(ECHO) "$(COLOR_YELLOW)Cleaning...$(COLOR_RESET)"
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@$(ECHO) "$(COLOR_GREEN)✓ Clean complete$(COLOR_RESET)"

.PHONY: rebuild
rebuild: clean all

