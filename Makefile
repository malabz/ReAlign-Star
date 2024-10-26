# Define compiler and compile options
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

# Define targets and dependencies
TARGET = realign_star
SRCS = src/main.cpp src/Fasta.cpp
OBJS = $(SRCS:.cpp=.o)
INSTALL_DIR = $(HOME)/.realign_star/bin

# Detect user's shell
SHELL_NAME := $(shell basename $$SHELL)

all: $(TARGET)

# Link object files
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile the source code
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Clean
clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all $(TARGET) clean
