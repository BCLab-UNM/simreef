# ============================================
# SubstrateGenerator Makefile
# ============================================

# Compiler and flags
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -march=native -ffast-math
LDFLAGS  :=

# Target and source
TARGET   := reef_spa_substrate
SRC      := reef_spa_bitmap.cpp

# Default build rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Clean up build artifacts
clean:
	rm -f $(TARGET) *.o

# Rebuild from scratch
rebuild: clean all

# Run a quick test (optional sanity check)
test: $(TARGET)
	./$(TARGET) test.bmp 256 256 1.0 0.0 3.0 0.25 0.25 0.25 0.25 42 100 30 4 0.5 2.0 0.3 0.002
	file test.bmp

.PHONY: all clean rebuild test

