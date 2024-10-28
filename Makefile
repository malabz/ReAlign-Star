# Define compiler and compile options
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

# Define targets and dependencies
TARGET = realign_star
SRCS = src/main.cpp src/Fasta.cpp
OBJS = $(SRCS:.cpp=.o)
INSTALL_DIR = $(HOME)/.realign_star/bin
JAR_FILE = profileAlignment.jar

# Detect user's shell
SHELL_NAME := $(shell basename $$SHELL)

all: $(TARGET) install_jar

# Link object files
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile the source code
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Install profileAlignment.jar in the INSTALL_DIR
install_jar:
	@mkdir -p $(INSTALL_DIR)
	@cp src/$(JAR_FILE) $(INSTALL_DIR)/$(JAR_FILE)
	@echo "Installed $(JAR_FILE) to $(INSTALL_DIR)"

# Clean
clean:
	rm -f $(TARGET) $(OBJS)
	rm -f $(INSTALL_DIR)/$(JAR_FILE)

.PHONY: all $(TARGET) clean install_jar
