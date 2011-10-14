@echo off

SET SWIG=swig2.0
SET SWIGFLAGS=-c++ -python -threads
SET CXX=%1
SET CFLAGS=%2
SET INC=%3

:: remove "" from input
for /f "useback tokens=*" %%a in ('%CFLAGS%') do set CFLAGS=%%~a
for /f "useback tokens=*" %%a in ('%INC%') do set INC=%%~a

SET ifile=%4
SET cxxfile=%ifile:.i=.cxx%

echo."%SWIG% %SWIGFLAGS% -o %cxxfile% %4"
%SWIG% %SWIGFLAGS% -o %cxxfile% %4
echo."%CXX% %CFLAGS% %INC% -c %cxxfile% -o %5"
%CXX% %CFLAGS% %INC% -c %cxxfile% -o %5

