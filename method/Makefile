CXX = nvcc
CXXFLAGS = -Iinclude -Xcompiler "-g" -lz
LDFLAGS = -lboost_program_options
TARGET = yamet

# Source files
CPP_SOURCES = code/*.cpp

# Final target executable
$(TARGET): $(CUDA_SOURCES) $(CPP_SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CPP_SOURCES) -o $(TARGET)

# Clean up
clean:
	rm -f $(TARGET)