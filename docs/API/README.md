# API


![API Model](./API.jpg)


### What is it?
---
**`[functionName]API`**

The C API provides access to a number of class member and standalone functions in the library.  Functions are exported as "extern C", and as such, are able to be accessed by Python ctypes, as well as other callers requiring C linkage.  The above image gives a general representation of how API references to class member functions work.  Naming convention is generally: `[functionName]API`.

*Note: In addition, the functions within `namespace utils` are also exported, but with C++ linkage.  These need to be called from C++ or through some other method that allows for this linkage.  Generally, these functions shouldn't need to be accessed through Python, hence the C++ linkage.  They are also not listed here, but documented in the [utils](./../utils/README.md) docs.*


### Use
---
Generally speaking, API functions take in a pointer to an instance of a C++ class (if it's an API wrapper for a class member function, which is *not* `create`) as the first argument, followed by the function arguments, if any.  If it is a class creation API wrapper, it takes the arguments required for that class's constructor and returns a pointer to the class instance created.  If it is a class destruction API wrapper, it takes in the pointer to the class and returns nothing.


### Simulation API Functions
```
//Access Functions
double      getSimulationTimeAPI(Simulation* simulation);
double      getDtAPI(Simulation* simulation);
double      getSimMinAPI(Simulation* simulation);
double      getSimMaxAPI(Simulation* simulation);
int         getNumberOfParticleTypesAPI(Simulation* simulation);
int         getNumberOfParticlesAPI(Simulation* simulation, int partInd);
int         getNumberOfAttributesAPI(Simulation* simulation, int partInd);
const char* getParticleNameAPI(Simulation* simulation, int partInd);
const char* getSatelliteNameAPI(Simulation* simulation, int satInd);
LogFile*    getLogFilePointerAPI(Simulation* simulation);
const double* getPointerToParticleAttributeArrayAPI(Simulation* simulation, int partIndex, int attrIndex, bool originalData);

//Field Calculation
double getBFieldAtSAPI(Simulation* simulation, double s, double time);
double getEFieldAtSAPI(Simulation* simulation, double s, double time);

//Simulation Management Functions
Simulation* createSimulationAPI(double dt, double simMin, double simMax, const char* rootdir);
void initializeSimulationAPI(Simulation* simulation);
void __iterateSimCPUAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts);
void iterateSimulationAPI(Simulation* simulation, int numberOfIterations, int itersBtwCouts);
void freeGPUMemoryAPI(Simulation* simulation);
void terminateSimulationAPI(Simulation* simulation);
Simulation* loadCompletedSimDataAPI(const char* fileDir);
void setupExampleSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir);
void runExampleSimulationAPI(Simulation* sim, int iterations, int printEvery);
void runSingleElectronAPI(Simulation* sim, double vpara, double vperp, double s, double t_inc, int iterations, int printEvery);

//Fields Management
void setBFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString); //switch to comma delimited string of variables
void addEFieldModelAPI(Simulation* sim, const char* modelName, const char* doubleString);

//Particle Management
void createParticleTypeAPI(Simulation* simulation, const char* name, double mass, double charge, long numParts, const char* loadFileDir = "");

//Satellite Management
void    createSatelliteAPI(Simulation* simulation, int particleInd, double altitude, bool upwardFacing, const char* name);
int     getNumberOfSatellitesAPI(Simulation* simulation);
const double* getSatelliteDataPointersAPI(Simulation* simulation, int satelliteInd, int msmtInd, int attributeInd);

//CSV Writing
void writeCommonCSVAPI(Simulation* simulation);
```

*Note: this section refers in general to most of the functions in the above list.  Select functions that are worthy of further documentation will be listed immediately below this section.*

#### Input:

`simulation` - the pointer to the simulation instance created by `createSimulationAPI`

`others` - see documentation of [Simulation](./../Simulation/README.md) member functions for more specifics on other arguments to the various functions


#### Output:
See documentation of [Simulation](./../Simulation/README.md) member functions for more specifics on the output of the various functions.


#### Side-Effects:
See documentation of [Simulation](./../Simulation/README.md) member functions for more specifics on the side effects of the various functions.


---
```
void setupExampleSimulationAPI(Simulation* sim, int numParts, const char* loadFileDir);
void runExampleSimulationAPI(Simulation* sim, int iterations, int printEvery);
```
#### Inputs:
`sim` - described above

`numParts` - number of particles that Particle will contain

`loadFileDir` - where particle data will be loaded from.  *Note: name of data files must be preceded by "elec", as this is what the particle will be called.*

`iterations` - how many iterations to run the sim for

`printEvery` - how often to print a status report on the simulation.  Also checks whether all particles have escaped while printing.  A lower number here will lead to more printed status updates, and check whether all particles have escaped more often, incurring a slight performance penalty.  When in doubt, go with 500 or 1000.


#### Outputs:
None


#### Side-Effects:
Adds DipoleBLUT BField Model at ILAT 72.0 degrees, no EField Model, a Particle set named "elec" with the characteristics of the electron, and Satellites: 1 at sim top, 1 at sim bottom, 1 at 4e6 dist "s" along field line with detector facing up, 1 at same s with detector facing down, and 1 each of up/down detector at 3e6 s.


[Up a level](./../README.md)