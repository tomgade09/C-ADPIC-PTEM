# geoplasmasim

Developer: Tom Gade

geoplasmasim is a test particle simulation attempting to reproduce electron distributions gathered by FAST under the influence of a number of physical phenomena (mirror force, quasi-static potential structures, alfven waves, etc).  It is written with C++/CUDA and has convenient Python scripts for running.


## Compatibility
Right now, geoplasmasim runs on Windows and is being tested on Linux.  As it is dependent on CUDA, it requires an NVIDIA graphics card with up-to-date(-ish, at least) drivers, as well as the CUDA libraries installed.  No other external dependencies for C++ exist.  Matplotlib and numpy are required for running in Python 3, if desired.


## Dependencies
CUDA; Python 3 x64 (Python 2 and 32 bit are not supported, but may work), matplotlib, numpy (for running through Python) (see above - Compatibility)

## Getting Started

#### 1. Download Repository

##### Platform Agnostic

  ```
  git clone https://github.com/tomgade09/geoplasmasim
  cd geoplasmasim
  ```
  
#### 2. Compile

##### Windows

Open the Visual Studio solution, ensure Release and x64 is selected, and click Build -> Build Solution.

##### Linux
  
  ```
  ./configure
  make
  ```

Note: gcc compatible with the `-std=c++14` flag is required.

#### 3. Run an example simulation

From a terminal in `%WHEREVER%/geoplasmasim`

  ```
  cd python
  python ./simulation.py
  ```
  
The Python script will create the appropriate directories for you (`geoplasmasim/_dataout/%DATE-TIME-GROUP%`) and save data after the fact.  See the documentation for the save output folder structure.

## Additional Documentation
[Read the documentation here](./docs/README.md)