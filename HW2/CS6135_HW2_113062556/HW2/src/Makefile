# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -std=c++17 -O2 -fopenmp

# Boost library path
BOOST_INCLUDE = /usr/local/include/boost
BOOST_LIB = /usr/local/lib
LDFLAGS = -L$(BOOST_LIB)

# Target executables and directory
TARGET_DIR = ../bin
TARGET = $(TARGET_DIR)/hw2
TARGET_PAR = $(TARGET_DIR)/hw2_parallel

# Source files and object files for main.cpp
SRC = main.cpp
OBJ = $(SRC:.cpp=.o)

# Source files and object files for main_p.cpp (parallel version)
SRC_PAR = main_p.cpp
OBJ_PAR = $(SRC_PAR:.cpp=.o)

# Default rule: build both executables
all: $(TARGET) $(TARGET_PAR)

# Compile main.cpp to object file
$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) -c $< -o $@

# Link object files to create hw2
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) $(OBJ) -o $@ $(LDFLAGS)

# Compile main_p.cpp to object file
$(OBJ_PAR): $(SRC_PAR)
	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) -c $< -o $@

# Link object files to create hw2_parallel
$(TARGET_PAR): $(OBJ_PAR)
	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) $(OBJ_PAR) -o $@ $(LDFLAGS)

# Clean up object files and executables
clean:
	rm -f $(OBJ) $(OBJ_PAR) $(TARGET) $(TARGET_PAR)
