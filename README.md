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

### Link to BioPARKIN

Presumingly, we want to make use of all the numerical code in BioPARKIN, the sister package of PARKINcpp, there is a last step to do (which, in priciple, can not handled automatically, due to the splitting of the two packages BioPARKIN and PARKINcpp) :
````
cp swig/parkin.py [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
cp swig/_parkin.so [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
````

Alternatively, it is possible to set symbolic links, as well :
````
cd [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
ln -s [Installation path of PARKINcpp]/PARKINcpp/build/swig/parkin.py .
ln -s [Installation path of PARKINcpp]/PARKINcpp/build/swig/_parkin.so .
````

Please mind both methods are possible, but mutually exclusive, so there is a choice: The first is more static, promissing more repeatable results, even if something in PARKINcpp would have changed, while the second is more dynamic, reflecting immediately and automatically(!) all changes after recompilation of PARKINcpp.  You have been informed!
