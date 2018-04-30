# geoplasmasim

Developer: Tom Gade

geoplasmasim is a test particle simulation attempting to reproduce electron distributions gathered by FAST under the influence of a number of physical phenomena (mirror force, quasi-static potential structures, alfven waves, etc).  It is written with C++/CUDA and has convenient Python scripts for running.


## Compatibility
Right now, geoplasmasim runs on Windows and is being tested on Linux.  As it is dependent on CUDA, requires an NVIDIA graphics card with up-to-date(-ish, at least) drivers.  No other external dependencies for C++ exist.  Matplotlib is required for Python.


## Dependencies
CUDA, matplotlib (for running through Python) (see above - Compatibility)

## Getting Started

#### Download Repository

  ```
  git clone https://github.com/tomgade09/geoplasmasim
  cd geoplasmasim
  ```

Open the Visual Studio solution, ensure Release and x64 is selected (works in x86 as well), and click Build -> Build Solution.


#### Run an example simulation
From a terminal in `%WHEREVER%/geoplasmasim`

  ```
  cd python
  python ./simulation.py
  ```
  
The Python script will create the appropriate directories for you (`geoplasmasim`/_dataout/%DATE-TIME-GROUP%`) and save data after the fact.  See the documentation for the save output folder structure.

## Additional Documentation
[Read the documentation here](./docs/README.md)