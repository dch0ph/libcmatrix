# Generated by configure from Makefile.in
# Avoid editing this file directly
# Makefile for libcmatrix

CPPFLAGS= -msse3 -mfpmath=sse -I/home/devuser/libcmatrix/include -DLCM_SUPPRESS_VIEWS -DLCM_DISABLE_IMPLICIT_COMPLEX
CXX=g++
CXXFLAGS= -Wall -Wno-sign-compare -DLCM_USE_EXTERNAL
RANLIB=ranlib
# flags for optimisation
OPTFLAGS=-O2
# method to create archive
AR=ar ru

# You shouldn't need to alter anything below here

COMPILE=$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c
ALLOPTFLAGS=$(OPTFLAGS) -DNDEBUG

LOCALOBJS= local/CrystalSystem.o local/CrystalGeneric.o local/MoleculeStructure.o
COREOBJS= coredefs/complex.o coredefs/common.o coredefs/diagonal.o coredefs/diagonalise.o coredefs/power.o coredefs/invert.o coredefs/transforms.o coredefs/blocking.o coredefs/mixed.o coredefs/realtransforms.o coredefs/Chebyshev.o utils/ft.o utils/geometry.o utils/Lineshapes.o utils/constraints.o
NGPLOBJS=ngpl/Interpolate.o ngpl/simplex.o ngpl/levmarq.o 
NGPLUNSAFEOBJS=ngpl/statstuff.o
UNSAFEOBJS= coredefs/cmatrix.o coredefs/rmatrix.o NMR/MetaPropagation.o NMR/Sequence.o
NMROBJS= NMR/space_T.o NMR/wigner.o NMR/MAS.o NMR/InhomoMAS.o NMR/nonsecular.o NMR/proc.o NMR/superop.o NMR/HomoPropTD.o NMR/HomoPropFD.o NMR/InhomoProp.o NMR/StaticProp.o NMR/AsynchronousProp.o NMR/PhaseModulatedPropagator.o NMR/InversionGenerator.o
NONOPTOBJS = coredefs/writematrix-n.o coredefs/readmatrix-n.o utils/matlabio-n.o utils/readmatlabio-n.o utils/matlabcomposite-n.o utils/ttyio-n.o utils/timer-n.o NMR/spin_system-n.o  NMR/powder-n.o NMR/basespin_system-n.o NMR/NMR-n.o NMR/tensorop-n.o NMR/spinhalf_system-n.o utils/simpsonio-n.o
ALWAYSDEBUGOBJS = coredefs/common_g.o
#@MINUITOBJS@

GPLSAFEOBJS= $(COREOBJS) $(NMROBJS) $(LOCALOBJS) utils/Fork_controller.o $(LOCALOBJS)
ALLOBJS= $(UNSAFEOBJS) $(NGPLUNSAFEOBJS) $(GPLSAFEOBJS) $(NGPLOBJS)

ALLNORMOBJS=$(ALLOBJS) $(NONOPTOBJS) $(ALWAYSDEBUGOBJS)
ALLGPLOBJS= $(UNSAFEOBJS) $(GPLSAFEOBJS) $(NONOPTOBJS) $(ALWAYSDEBUGOBJS)
ALLPROFOBJS=$(ALLOBJS:.o=_p.o) $(NONOPTOBJS) $(ALWAYSDEBUGOBJS)
ALLGOBJS=$(ALLOBJS:.o=-g.o) $(NONOPTOBJS:-n.o=-g.o) $(ALWAYSDEBUGOBJS)
ALLTHROBJS=$(UNSAFEOBJS:.o=_r.o) $(ALLSAFEOBJS)

ROOT=..

first: lib/libcmatrix.a

#%.o: %.f
#	@F77@ @FFLAGS@ -c -o $@ $<

%_p.o: %.cc
	$(COMPILE) $(ALLOPTFLAGS) -pg -o $@ $<

%_g.o: %_g.cc
	$(COMPILE) -g -o $@ $<

%-g.o: %.cc
	$(COMPILE) -g -o $@ $<

%-n.o: %.cc
	$(COMPILE) -DNDEBUG -o $@ $<

%_r.o: %.cc
	$(COMPILE) $(ALLOPTFLAGS)  -DLCM_REENTRANT -o $@ $<

%.o: %.cc
	$(COMPILE) $(ALLOPTFLAGS) -o $@ $<

# clean up template repositories
tempclean:
#	-@CLEAN@
#	-cd coredefs; @CLEAN@
#	-cd utils; @CLEAN@
#	-cd NMR; @CLEAN@
#	-cd optim; @CLEAN@
#	-cd local; @CLEAN@

lib/libcmatrix.a: $(ALLNORMOBJS)
	$(AR) $@ $(ALLNORMOBJS)
	chmod a+rx $@
	$(RANLIB) $@

lib/libcmatrix_p.a: $(ALLPROFOBJS)
	$(AR)  $@ $(ALLPROFOBJS)
	chmod a+rx $@
	$(RANLIB) $@

lib/libcmatrix-g.a: $(ALLGOBJS)
	$(AR)  $@  $(ALLGOBJS)
	chmod a+rx $@
	$(RANLIB) $@

lib/libcmatrix_r.a: $(ALLTHROBJS)
	$(AR)  $@ $(ALLTHROBJS)
	chmod a+rx $@
	$(RANLIB) $@

lib/libcmatrix_gpl.a: $(ALLGPLOBJS)
	$(AR)  $@ $(ALLGPLOBJS)
	chmod a+rx $@
	$(RANLIB) $@

libcmatrix.tar.gz:
	cd $(ROOT) ; tar --exclude lib/* --exclude *.o --exclude=*.aux --exclude=*.log --exclude=*~ --exclude=Makefile --exclude=include/config.h --exclude=config.status -cvf libcmatrix.tar libcmatrix
	gzip -f $(ROOT)/libcmatrix.tar

libcmatrix_lite.tar.gz:
	cd $(ROOT) ; tar --exclude lib/* --exclude=*~ --exclude=*.aux --exclude=*.log --exclude=Makefile --exclude=include/config.h --exclude=config.status --exclude=local -cvf libcmatrix_lite.tar libcmatrix
	gzip -f $(ROOT)/libcmatrix_lite.tar

libcmatrix_gpl.tar.gz:
	cd $(ROOT) ; tar --exclude lib/* --exclude=*~ --exclude=*.aux --exclude=*.log --exclude=Makefile --exclude=include/config.h --exclude=config.status --exclude=ngpl -cvf libcmatrix_gpl.tar libcmatrix
	gzip -f $(ROOT)/libcmatrix_gpl.tar

archives: distclean libcmatrix.tar.gz libcmatrix_lite.tar.gz libcmatrix_gpl.tar.gz

normclean: tempclean
	-rm lib/libcmatrix.a $(ALLOBJS)

pclean: tempclean
	-rm lib/libcmatrix_p.a $(ALLPROFOBJS)

gclean: tempclean
	-rm lib/libcmatrix-g.a $(ALLGOBJS)

rclean: tempclean
	-rm lib/libcmatrix_r.a $(ALLTHROBJS)

testclean:
	-cd test; make clean

docclean:
	-cd docs; rm *.dvi *.log *.ilg

clean: normclean pclean gclean rclean testclean docclean
	-rm -r *.cache config.log

distclean: tempclean
	-rm $(ALLOBJS) $(ALLPROFOBJS) $(ALLGOBJS) $(ALLTHROBJS)
	-rm -r test/Makefile Makefile_template .make.* *.cache config.log

all: lib/libcmatrix.a lib/libcmatrix-g.a 