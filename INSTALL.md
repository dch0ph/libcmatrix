
libcmatrix should in principle work with any Unix-like environment and an ANSI C++ compiler. It has only been thoroughly tested with GCC compilers, currently in a Ubuntu environment, and the configure script will almost certainly need tweaking for other compilers. 

**Windows**: It is likely that compilation will work for the Msys2 environment (win64 binaries). But given the availability of Windows Subsystem for Linux, there is litle incentive for supporting native Windows builds.

1. Adding support functionality

It is recommended to use a matrix support library if needing to deal with relatively large spin systems, . OpenBLAS seems to work out of the box easily (unlike ATLAS - the previously used external library). In a Ubuntu-like environment, this can be easily installed using `apt-get install openblas-dev`. You will need to know where the header and library files have been installed. 
Under current Ubuntu's, the relevant locations are `/usr/include/x86_64-linux-gnu` and `/usr/lib/x86_64-linux-gnu` respectively.

As noted below, Minuit support can be added, but is a little tricky, and is only relevant to a few use cases.

2. Configuration

`./configure` will run the autoconf script which creates `Makefile`, `include/config.h`, `Makefile_template` and `test/Makefile`
as a function of the system configuration.  Rather than tinker with the Makefiles produced, you should alter
compiler flags etc. by passing environment variables to `configure` e.g.
`CXX=g++ CXXFLAGS="-g +w2" ./configure`
(bash and related shells)
  
Recommended options (for optimisation): `-g -O3`  \[`-W -Wno-sign-compare` options are added automatically if gcc is being used\]

Optional components:

`--with-openblas` in the configure command line enables use of the OpenBLAS libraries for fast matrix multiplication etc. If the library has been compiled with an external library, then all programs  that use the library must be compiled these options.
The location of the header and library files must be set if these are not already in the compile/link path e.g.
`CPPFLAGS="-I<dir>/include" LDFLAGS="-L<dir>/lib" ./configure`

`--with-atlas` enables use of the ATLAS libraries for fast matrix multiplication etc. 
This quite challenging to compile and is no longer recommended.

`--with-acml` `--with-mkl` are alternatives to ATLAS/OpenBLAS for optimised for AMD- and Intel-based processors respectively.


`--with-threads` enables use of multi-threading.  A special version 
of library `libcmatrix_r.a` must be compiled, together with appropriate
options when compiling (see `Makefile_template`). This is not being actively developed/verified.
MPI (see below) is a simpler option in most cases. It is also not compatible with Windows.

`--with-MPI` enables use of the Message Passing Interface (MPI).
`mpi.h` and the `mpi` library must be must be accessible for this to work.

`--with-minuit` enables Minuit support for more complex optimisation. 
If `--with` or `--without` is not set explicitly, it is automatically included if the appropriate files are found. 
The configuration script checks for the current API (Minuit2). It no longer works with the older C++ Minuit API. Minuit is now quite hard to find. See below for some compilation tips. Recommended to try working without Minuit first.

`--with-sse` enables the use of SSE instructions, which would historically speed up complex arithmetric by a factor of ~2. This introduces more home-brewed code, and the SSE instructions are specific to Intel-type. Currently not recommended.

3. `make lib/libmatrix.a` will create the libcmatrix.a library. Other variants exist for debugging (-g) and profiling (-pg).  

If only compiling libcmatrix for linking with other programs such as pNMRsim, stop here.

There is no install script since libcmatrix is often only temporarily required for pNMRsim compilation.  If preferred, however, the contents of `lib` and `include` can be installed in say `/usr/local` and `/usr/local/include/libcmatrix` respectively.


4. Testing
`make all` in the `test` subdirectory will create a set of test programs.
These are mostly demonstrations of different aspects of the library rather than real stress tests.  Many programs require input; comments in the .cc files should give suitable values.  Don't worry if not all the programs compile - the libcmatrix API currently evolves for best interaction with pNMRsim and so the "standalone" interface is not so well defined.

**Optimisation**: If peak performance is important, the performance of key operations
such a matrix multiplication can be checked to ensure that appropriate "tuning"
options have been selected.  The `testops` programs determines the floating point "throughput" (in megaflops/s) of various
operations as a function of matrix size.  The peak performance should be of the same order as the CPU clock speed.
Note that for every operation that uses the external libtaries, there is a "cross-over" point,
with the external libtary only being used for larger matrices / vectors. `testops` can be used
to search for these cross-over points, which will be architecture dependent. The fact that 
external libraries are only effective in some cases is a motivation for avoiding them if not needed!



## Note on installing Minuit2

A standalone version of Minuit2 can be found [here](https://github.com/GooFit/Minuit2). The `cmake` build system is required (`apt-get install cmake` on Ubuntu). 
 
Unzip into `Minuit2-master`. Had to edit `CMakeLists.txt` to add 
`link_libraries("-lstdc++")` 
to get C++ programs to link properly .

```
mkdir build 
cd build
cmake ..  
```

The compiled libraries and headers can be installed, but it may be simpler just to add the relevant directories (`lib` and `inc` respectively) to the paths when compiling `libcmatrix`, `pNMRsim` etc.
