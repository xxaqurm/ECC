CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++17
TARGET = main.out
BUILD_DIR = build
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:src/%.cpp=$(BUILD_DIR)/%.o)

all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean