# References

# Compiler
CC = g++

# Compiler flags
CFLAGS = -std=c++11 -Wall -g3

# Source file and target executable
SRC = main.cpp
BUILD_DIR = build
TARGET = $(BUILD_DIR)/main.out

all: $(TARGET)

$(TARGET): $(SRC)
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -rf $(BUILD_DIR)
