Click the links below to read about the various components of this software:

[1. Simulation](Simulation/README.md) - Container class that contains equation of motion for particles, as well as a number of control functions, access functions, and other useful things.  Manages lifetime of all of the below classes.

[2. API](API/README.md) - Mostly extern C-style API for interfacing through Python or other language.

[3. BField](BField/README.md) - Abstract class for interfacing with various implementations of B Field models.

[4. EField](EField/README.md) - Abstract class for interfacing with various implementations of E Field models.  Can track numerous `EElem`s - E Field Elements.

[5. Particle](Particle/README.md) - Class that manages particles by tracking arrays of attributes, specified upon creation.  Also manages on GPU data, including cudaMalloc/cudaFree on initialization/destruction respectively.

[6. Satellite](Satellite/README.md) - Class that tracks whether or not a particle has passed a certain altitude from above or below.  Also manages on GPU data, including cudaMalloc/cudaFree on initialization/destruction respectively.

[7. FileIO](FileIO/README.md) - Namespaced functions that handle reading and writing various filetypes to/from disk.

[8. utils](utils/README.md) - Namespaced functions that accomplish a number of useful things, including generating a distribution of particles, generating a CSV from specified data, and loading/comparing particle distributions (including initial and final data, and satellite data).

[9. ErrorHandling](ErrorHandling/README.md) - Mostly macros to handle various error conditions, as well as make sure certain conditions are met.

[10. LogFile](LogFile/README.md) - A class with a number of member functions for writing a log file to disk along with elapsed time data.

[11. Examples](Examples/README.md) - Look here for examples of usage (Python and C++).

*Note: In this documentation, uppercase (and usually linked) names refer to classes, while lowercase names refer to non-class things.  For example: [Particle](Particle/README.md) refers to the class itself or an instance of the class which manages a large number of particles (lowercase).  particle(s) usually refers to a collection of attributes (ex: v_para, v_perp or mu, and s, as well as maybe time, index, etc) that represents a `real-world physical particle`.*