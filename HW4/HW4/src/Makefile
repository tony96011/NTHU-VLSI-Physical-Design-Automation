CXX = g++
CXXFLAGS = -Wall -std=c++17 -O2

# Boost library path
BOOST_INCLUDE = /usr/local/include/boost
BOOST_LIB = /usr/local/lib
LDFLAGS = -L$(BOOST_LIB)

# Target directory and executable
TARGET_DIR = ../bin
TARGET = $(TARGET_DIR)/hw4

SRCS = main.cpp Floorplan.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

# Linking
# $(TARGET): $(OBJS)
# 	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) -o $@ $(OBJS) $(LDFLAGS)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile each cpp
# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) -I$(BOOST_INCLUDE) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

