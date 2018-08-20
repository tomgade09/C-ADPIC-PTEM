# Compiler binaries and settings
CC         := g++
NVCC       := nvcc
CXXFLAGS   := -std=c++17 -fPIC -pedantic -O2 -fopenmp -Wall -Wno-unused-variable
NVCFLAGS   := -std=c++14 -rdc=true -O2 -gencode=arch=compute_30,code=\"sm_30,compute_30\" -x cu -m64 -cudart static -Xcompiler "-fPIC" -Xlinker "-shared -fopenmp"

# Build-essential directories and defines
CUDAINC    := /usr/local/cuda/include
CUDALIB    := /usr/local/cuda/lib64
LIBS       := -L$(CUDALIB)
INCS       := -I$(CUDAINC) -I./include
SRC        := ./src
OUT        := ./lib
BUILD      := $(OUT)/build
TARGET     := $(OUT)/geoplasmasim.so

# Lists of files
SOURCES    := $(shell find $(SRC) -name "*.cpp")# | sed 's/^\.\/src\/SimAttributes\/SimAttributes\.cpp//')
OBJECTS    := $(patsubst $(SRC)/%.cpp,$(BUILD)/%.o,$(SOURCES))
CUSOURCES  := $(shell find $(SRC) -name "*.cu")
CUOBJECTS  := $(patsubst $(SRC)/%.cu,$(BUILD)/%.obj,$(CUSOURCES))

#Default rule
$(TARGET): $(OBJECTS) $(CUOBJECTS)
	$(NVCC) -shared -lcudadevrt -Xcompiler "-fopenmp" -Xlinker "-fopenmp" -o $(TARGET) $(OBJECTS) $(CUOBJECTS)


# File rules
$(BUILD)/Simulation/Simulation.o: $(SRC)/Simulation/Simulation.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/Simulation/Simulation.cpp -o $(BUILD)/Simulation/Simulation.o

$(BUILD)/Simulation/iterateSimCPU.o: $(SRC)/Simulation/iterateSimCPU.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/Simulation/iterateSimCPU.cpp -o $(BUILD)/Simulation/iterateSimCPU.o

$(BUILD)/Satellite/Satellite.o: $(SRC)/Satellite/Satellite.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/Satellite/Satellite.cpp -o $(BUILD)/Satellite/Satellite.o

$(BUILD)/Simulation/previousSimulation.o: $(SRC)/Simulation/previousSimulation.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/Simulation/previousSimulation.cpp -o $(BUILD)/Simulation/previousSimulation.o

$(BUILD)/API/utilsAPI.o: $(SRC)/API/utilsAPI.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/API/utilsAPI.cpp -o $(BUILD)/API/utilsAPI.o

$(BUILD)/API/LogFileAPI.o: $(SRC)/API/LogFileAPI.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/API/LogFileAPI.cpp -o $(BUILD)/API/LogFileAPI.o

$(BUILD)/API/SimulationAPI.o: $(SRC)/API/SimulationAPI.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/API/SimulationAPI.cpp -o $(BUILD)/API/SimulationAPI.o

$(BUILD)/LogFile/LogFile.o: $(SRC)/LogFile/LogFile.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/LogFile/LogFile.cpp -o $(BUILD)/LogFile/LogFile.o

$(BUILD)/SimAttributes/SimAttributes.o: $(SRC)/SimAttributes/SimAttributes.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/SimAttributes/SimAttributes.cpp -o $(BUILD)/SimAttributes/SimAttributes.o

$(BUILD)/utils/fileIO.o: $(SRC)/utils/fileIO.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/fileIO.cpp -o $(BUILD)/utils/fileIO.o

$(BUILD)/utils/random.o: $(SRC)/utils/random.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/random.cpp -o $(BUILD)/utils/random.o

$(BUILD)/utils/writeIOclasses.o: $(SRC)/utils/writeIOclasses.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/writeIOclasses.cpp -o $(BUILD)/utils/writeIOclasses.o

$(BUILD)/utils/PPHelperClasses.o: $(SRC)/utils/PPHelperClasses.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/PPHelperClasses.cpp -o $(BUILD)/utils/PPHelperClasses.o

$(BUILD)/utils/readIOclasses.o: $(SRC)/utils/readIOclasses.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/readIOclasses.cpp -o $(BUILD)/utils/readIOclasses.o

$(BUILD)/utils/string.o: $(SRC)/utils/string.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/string.cpp -o $(BUILD)/utils/string.o

$(BUILD)/utils/postprocess.o: $(SRC)/utils/postprocess.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/postprocess.cpp -o $(BUILD)/utils/postprocess.o

$(BUILD)/utils/numerical.o: $(SRC)/utils/numerical.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/utils/numerical.cpp -o $(BUILD)/utils/numerical.o

$(BUILD)/Particle/Particle.o: $(SRC)/Particle/Particle.cpp
	$(CC) $(CXXFLAGS) -c $(INCS) $(LIBS) $(SRC)/Particle/Particle.cpp -o $(BUILD)/Particle/Particle.o

#CUDA Source Files
$(BUILD)/Simulation/simulationphysics.obj: $(SRC)/Simulation/simulationphysics.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/Simulation/simulationphysics.cu -o $(BUILD)/Simulation/simulationphysics.obj

$(BUILD)/Satellite/Satellite.obj: $(SRC)/Satellite/Satellite.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/Satellite/Satellite.cu -o $(BUILD)/Satellite/Satellite.obj

$(BUILD)/EField/QSPS.obj: $(SRC)/EField/QSPS.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/EField/QSPS.cu -o $(BUILD)/EField/QSPS.obj

$(BUILD)/EField/EField.obj: $(SRC)/EField/EField.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/EField/EField.cu -o $(BUILD)/EField/EField.obj

$(BUILD)/Particle/Particle.obj: $(SRC)/Particle/Particle.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/Particle/Particle.cu -o $(BUILD)/Particle/Particle.obj

$(BUILD)/BField/DipoleB.obj: $(SRC)/BField/DipoleB.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/BField/DipoleB.cu -o $(BUILD)/BField/DipoleB.obj

$(BUILD)/BField/DipoleBLUT.obj: $(SRC)/BField/DipoleBLUT.cu
	$(NVCC) $(NVCFLAGS) -c $(INCS) $(LIBS) $(SRC)/BField/DipoleBLUT.cu -o $(BUILD)/BField/DipoleBLUT.obj


.PHONY: clean
clean:
	find ./lib -name "*.o" -type f -delete
	find ./lib -name "*.obj" -type f -delete
	rm $(TARGET)

