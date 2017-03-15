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

Presumingly, we want to make use of all the numerical code in BioPARKIN, the sister package of PARKINcpp, there is a last step to do (which, in principle, can not be handled automagically, due to the splitting of our two packages BioPARKIN and PARKINcpp into two repositories. But that's as it is...) :
````
cp swig/parkin.py [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
cp swig/_parkin.so [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
````

Alternatively, it is possible to set symbolic links instead :
````
cd [Installation path of BioPARKIN]/BioPARKIN/src/parkincpp
ln -s [Installation path of PARKINcpp]/PARKINcpp/build/swig/parkin.py .
ln -s [Installation path of PARKINcpp]/PARKINcpp/build/swig/_parkin.so .
````

Please note both methods are possible, but mutually exclusive. So there is a choice: The first is a bit more static, promising more repeatable results, even if something in PARKINcpp would have been changed, whereas the second can be considered to be a little bit more dynamic, reflecting immediately and automatically(!) all changes after any recompilation of PARKINcpp.  You have been informed!
