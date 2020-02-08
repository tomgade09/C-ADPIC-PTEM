# C-ADPIC-PTEM

Developer: Tom Gade

C-ADPIC-PTEM is a test particle simulation attempting to reproduce electron distributions gathered by FAST under the influence of a number of physical phenomena (mirror force, quasi-static potential structures, alfven waves, etc).  It is written with C++/CUDA and has convenient Python scripts for running.


## Compatibility
Right now, C-ADPIC-PTEM runs on Windows and is being tested on Linux.  As it is dependent on CUDA, it requires an NVIDIA graphics card with up-to-date(-ish, at least) drivers, as well as the CUDA libraries installed.  No other external dependencies for C++ exist.  Matplotlib and numpy are required for running in Python 3, if desired.


## Dependencies
CUDA; Python 3 x64 (Python 2 and 32 bit are not supported, but may work), matplotlib, numpy (for running through Python) (see above - Compatibility)

## Getting Started

#### 1. Download Repository

##### Platform Agnostic

  ```
  git clone https://github.com/tomgade09/C-ADPIC-PTEM
  cd C-ADPIC-PTEM
  ```
  
#### 2. Compile

##### Windows

Open the Visual Studio solution (`C-ADPIC-PTEM/vs/C-ADPIC-PTEM.sln`), ensure `Release` and `x64` is selected, and in the toolbar click `Build -> Build Solution`.  If post-processing is desired (Windows only), open `C-ADPIC-PTEM/utils/BinToCDF/vs/BinToCDF.sln`, ensure `Release` and `x64` is selected, and click `Build -> Build Solution`.

##### Linux
  
  ```
  ./configure
  make
  ```

Note: gcc compatible with the `-std=c++14` flag is required.


#### 3. Run an example simulation

##### 3a. Create a Particle Distribution
From a terminal in `%WHEREVER%/C-ADPIC-PTEM`

  ```
  cd python/utils
  python ./createParticleDistribution.py
  ```
  
If you open `createParticleDistribution.py` in a text editor or IDE, you'll see the pitch and energy ranges of the particle distribution.  You can modify these to your liking before running.  Five files should be created in a newly created folder `./out`.  Copy these to `C-ADPIC-PTEM/_in/data`


##### 3b. Run Simulation with Example Characteristics

  ```
  (from C-ADPIC-PTEM/python/utils)
  cd ..
  python ./simulation.py
  ```
  
The Python script will create the appropriate directories for you (`C-ADPIC-PTEM/_dataout/%DATE-TIME-GROUP%`) and save data after the fact.  See the documentation for the save output folder structure.  The characteristics of the simulation will be printed to the output (BField model, EField model, etc).


##### 3c. Run Post-Processing (Windows only, after compiling BinToCDF)

  ```
  cd ../lib
  ./BinToCDF.exe ".\..\_dataout\%DATE-TIME-GROUP%"
  ```

Replace `%DATE-TIME-GROUP%` with the name of the folder that you want to post-process.  This will output a CDF file in the same directory (`C-ADPIC-PTEM/lib`).  This CDF file can be further processed to the user's liking.  See documentation on the post-processing functions for more info about what data is output in the CDF file.


## Additional Documentation
[Read the documentation here](./docs/README.md)