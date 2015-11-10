@echo off
rem MEXOPTS.BAT
rem
rem    Compile and link options used for building MEX-files
rem    using NVCC and the Microsoft Visual C++ compiler version 10.0.
rem
rem    Copyright 2012 The MathWorks, Inc.
rem
rem StorageVersion: 1.0
rem C++keyFileName: NVCCOPTS.BAT
rem C++keyName: nvcc
rem C++keyManufacturer: NVIDIA
rem C++keyVersion: 
rem C++keyLanguage: C++
rem C++keyLinkerName: Microsoft Visual C++ 2010
rem C++keyLinkerVersion: 10.0
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set VSINSTALLDIR=%VS100COMNTOOLS%\..\..
set VCINSTALLDIR=%VSINSTALLDIR%\VC
rem In this case, LINKERDIR is being used to specify the location of the SDK
set LINKERDIR='.registry_lookup("SOFTWARE\Microsoft\Microsoft SDKs\Windows\v7.0A" , "InstallationFolder").'
rem We assume that the CUDA toolkit is already on your path. If this is not the
rem case, you can set the environment variable MW_NVCC_PATH to the place where
rem nvcc is installed.
set PATH=%MW_NVCC_PATH%;%VCINSTALLDIR%\bin\amd64;%VCINSTALLDIR%\bin;%VCINSTALLDIR%\VCPackages;%VSINSTALLDIR%\Common7\IDE;%VSINSTALLDIR%\Common7\Tools;%LINKERDIR%\bin\x64;%LINKERDIR%\bin;%MATLAB_BIN%;%PATH%
rem Include path needs to point to a directory that includes gpu/mxGPUArray.h
set INCLUDE=%VCINSTALLDIR%\INCLUDE;%VCINSTALLDIR%\ATLMFC\INCLUDE;%LINKERDIR%\include;%INCLUDE%;%MATLAB%\toolbox\distcomp\gpu\extern\include
rem extern\lib\win64 points to gpu.lib: CUDA_LIB_PATH points to cudart.lib
set LIB=%VCINSTALLDIR%\LIB\amd64;%VCINSTALLDIR%\ATLMFC\LIB\amd64;%LINKERDIR%\lib\x64;%MATLAB%\extern\lib\win64;%LIB%;%CUDA_LIB_PATH%
set MW_TARGET_ARCH=win64

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=nvcc
set COMPFLAGS=-gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=\"sm_30,compute_30\" -c --compiler-options=/GR,/W3,/EHs,/D_CRT_SECURE_NO_DEPRECATE,/D_SCL_SECURE_NO_DEPRECATE,/D_SECURE_SCL=0,/DMATLAB_MEX_FILE,/nologo,/MD
set OPTIMFLAGS=--compiler-options=/O2,/Oy-,/DNDEBUG
set DEBUGFLAGS=--compiler-options=/Z7

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
rem Link with the standard mex libraries and gpu.lib.
set LIBLOC=%MATLAB%\extern\lib\win64\microsoft
set CUDALIBLOC=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v6.5\lib\x64
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib gpu.lib /LIBPATH:"%CUDALIBLOC%" cudart.lib cuda.lib /MACHINE:X64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /manifest /incremental:NO /implib:"%LIB_NAME%.x" /MAP:"%OUTDIR%%MEX_NAME%%MEX_EXT%.map"
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%LIB_NAME%.x" "%LIB_NAME%.exp"
set POSTLINK_CMDS1=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%;2" -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS2=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.map"
