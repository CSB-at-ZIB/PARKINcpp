# PARKINcpp

## Quickstart

Basic instructions to get started with the code.

### Dependencies

Required:
* [cmake](http://www.cmake.org/cmake/resources/software.html)
* [SWIG](http://www.swig.org/)

MacOsX:  
Install system requirements with [homebrew](http://mxcl.github.com/homebrew/)
````
brew install gfortran
brew install swig
````

Windows:  
* [MinGW](http://www.mingw.org/)

### Build instructions:

Assuming that we want the binaries installed into a "build" directory at the root of the PARKINcpp tree :
````
cd PARKINcpp
mkdir build
cd build
````

Linux/OSX:
````
cmake -G "Unix Makefiles" ..
make
````

NOTE: On Windows you have to use the MinGW shell.

