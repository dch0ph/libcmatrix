
libcmatrix should in principle work with any Unix-like environment and an ANSI C++ compiler. It has only been thoroughly tested with GCC compilers and the configure script will almost certainly need tweaking for other compilers.  


Windows: Compilation has most recently been tested for the Msys2 environment (win64 binaries). 


Configuration
-------------

`./configure` will run the autoconf script which creates `Makefile`, `config.h`, `Makefile_template` and `test/Makefile`
as a function of the system configuration.  Rather than tinker with the Makefiles produced, you should alter
compiler flags etc. by passing environment variables to configure e.g. `CXX=g++ CXXFLAGS="-g +w2" ./configure`
(bash and related shells).
  
Recommended options (for optimisation):
`g++` with `-O3`  [`-W` and `-Wno-sign-compare` options are added automatically if gcc is being used.]

Optional components:

`--with-openblas` in the configure command line enables use of the OpenBLAS
libraries for fast matrix multiplication etc. (recommended if, but only if,
dealing with relatively large spin systems). OpenBLAS seems to work out of the
box easily (unlike ATLAS - the previously used external library).

The location of the header and library files must be set if these are not already in the compile/link path e.g.
`CPPFLAGS="-I<dir>/include" LDFLAGS="-L<dir>/lib" ./configure`

If the library has been compiled with an external library, then all programs 
that use the library must be compiled these options.

`--with-atlas` enables use of the ATLAS libraries for fast matrix multiplication etc. 
This quite challenging to compile and is no longer recommended.

`--with-acml` or `--with-mkl` are alternatives to ATLAS for AMD- and Intel- based processors respectively.
Again `--with-openblas` is the recommended route.


`--with-threads` enables use of multi-threading.  A special version 
of library libcmatrix_r.a must be compiled, together with appropriate
options when compiling (see `Makefile_template`). Multi-threading
is not being actively developed/verified.  MPI (see below) is
a simpler option in most cases. This is not compatible with Windows
and is not recommended.

`--with-MPI` enables use of the Message Passing Interface (MPI).
`mpi.h` and `-lmpi` must be accessible for this to work. It is also generally necessary to use a wrapper script from compilation. 

`--with-minuit` enables Minuit support for more complex optimisation. 
If `--with` or `--without` is not set explicitly, it is automatically included if the appropriate files are found. 
The configuration script checks for the current API (Minuit2). Note that Minuit is now quite hard to find. See below for some compilation tips. Recommended to try working without Minuit first.

`--with-sse` enables the use of SSE instructions to speed up complex arithmetric by a factor of ~2. This introduces more
home-brewed code, but has been used now without problems for many years. The SSE instructions are specific to Intel-type
processors. 

Check the configuration options in `config.h`, but be sure you know what you are doing before changing anything!

Note that you can recreate the `Makefile`s etc. after edits of say `config.h` without re-running configure by running `config.status`.

Building
--------

`make lib/libmatrix.a` 
will create the libcmatrix library.

Other variants exist for debugging (-g) and profiling (-pg).  


If only compiling libcmatrix for linking with other programs such as `pNMRsim`, stop here.

There is no install script since libcmatrix is often only temporarily required for pNMRsim compilation.  If preferred, however, the contents of lib and include can be installed in say `/usr/local` and `/usr/local/include/libcmatrix` respectively.


Testing
-------

`make all` in the `test` subdirectory which will create a set of test programs.
These are mostly demonstrations of different aspects of the library rather than real stress tests.  Many programs require input; comments in the .cc
files should give suitable values.  Don't worry if not all the programs compile - the libcmatrix API currently evolves for best interaction with `pNMRsim` and so the "standalone" interface is not so well defined.

It is a good idea to check the performance of key operations
such a matrix multiplication to ensure that appropriate "tuning"
options have been selected.  The `testops` programs determine
the floating point "throughput" (in megaflops/s) of various
operations as a function of matrix size.  The peak performance
should be of the same order as the CPU clock speed.
Note that there is a "cross-over" point for every operation that uses the external libraries,
with the external libtary only being used for larger matrices / vectors. testops can be used
to search for these cross-over points, which will be architecture dependent. The fact that 
external libraries are only effective in some cases is a motivation for avoiding them if
not needed!


Read the documentation in `docs`!

<paul.hodgkinson@durham.ac.uk>


Note on installing Minuit2
--------------------------

A standalone version of Minuit2 can be obtained [here](https://github.com/GooFit/Minuit2). 

Download a zip file of the code and unzip into `Minuit2-master`, e.g. at the same directory level as `libcmatrix`.

`mkdir build` inside `Minuit2-master` to create the build directory. Inside this directory, run `cmake ..` 

Inside build: 

cmake ..  
to create Makefiles 

Had to edit CMakeLists.txt to add 

link_libraries("-lstdc++") 

to get C++ programs to link properly 

make install 

Install the header files and library. 



