CXX = nvcc
CXXFLAGS = -Iinclude -Xcompiler "-g" -lz
LDFLAGS = -lcudart -lboost_program_options
TARGET = yamet

# Source files
CUDA_SOURCES = code/*.cu
CPP_SOURCES = code/*.cpp

# Final target executable
$(TARGET): $(CUDA_SOURCES) $(CPP_SOURCES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CUDA_SOURCES) $(CPP_SOURCES) -o $(TARGET)

# Clean up
clean:
	rm -f $(TARGET)