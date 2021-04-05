AC_INIT
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([include/config.h])

DEFAR="ar ru"
MAKEEXTRA=

m4_include(configure-compiler.ac)

NEEDG2C="no"

m4_include(configure-libraries.ac)

if test -z "$AR"; then
	AR=$DEFAR
fi

# removed configure test for constructor elision (moved to contrib/configure_old.ac) from here 19/1/16

EXTRAOBJS=""
if test -d "local"; then
	EXTRAOBJS="\$(LOCALOBJS)"
	AC_MSG_NOTICE([Compiling with local/periodic extensions])
else
	AC_MSG_WARN([Compiling without local extensions])
fi

if test x$with_mpi = xyes; then
	EXTRAOBJS="utils/cmatrix_MPI.o $EXTRAOBJS"
fi

AC_PROG_RANLIB
AC_SUBST(UNSAFELIBS)
AC_SUBST(OPTFLAGS)
AC_SUBST(AR)
#AC_CHECK_LIB(sunmath,sincos)
#AC_CHECK_HEADERS(sunmath.h)
#AC_CHECK_HEADER(sys/times.h)
#AC_CHECK_HEADER(sstream)

AC_CHECK_HEADER(sys/wait.h,
	AC_DEFINE(HAVE_FORK_CONTROLLER)
	EXTRAOBJS="utils/Fork_controller.o $EXTRAOBJS",
	AC_MSG_WARN("Warning: Fork_controller not enabled as usable fork/wait not found"))

AC_CHECK_HEADER(pthread.h,AC_DEFINE(HAVE_PTHREAD_H)
HAVEPTHREAD="yes",HAVEPTHREAD="no")

AC_ARG_WITH(threads,
	AS_HELP_STRING([--with-threads],[enable use of multi-threading (default=no)]),,with_threads=no)
if test x$with_threads = xyes; then
	if test $HAVEPTHREAD = no; then
		AC_MSG_FAILURE([pthread.h not found
Check CPPFLAGS (or CXX) is defined to include relevant directory])
	fi
	AC_CHECK_LIB(pthread,pthread_cond_wait,AC_DEFINE(HAVE_LIBPTHREAD),AC_MSG_FAILURE([pthread library (libpthread.a) not found
Check LDFLAGS is defined to include relevant directory]))
	MAKEEXTRA="lib/libcmatrix_r.a $MAKEEXTRA"
	EXTRAOBJS="coredefs/cmatrix_threads.o $EXTRAOBJS"
	AC_SUBST(MTFLAGS)
	AC_MSG_WARN("Warning: multi-threading support is not being maintained - no service guarantees!");
fi
AC_SUBST(MAKEEXTRA)
AC_SUBST(EXTRAOBJS)

#AC_CHECK_FUNCS(sysconf)
AC_CHECK_FUNCS(sincos isfinite)
AC_CHECK_SIZEOF(float,4)
if test "$ac_cv_sizeof_float" != 4; then
AC_MSG_WARN([float size <> 4.  Read/write of binary files won't work!])
fi 
AC_CHECK_SIZEOF(double,8)
if test "$ac_cv_sizeof_double" != 8; then
AC_MSG_WARN([double size <> 8.  Read/write of binary files won't work!])
fi
AC_EGREP_HEADER(u_int32_t,sys/types.h,AC_DEFINE(NEED_UNDERSCORE))
#AC_CHECK_PROG(CLEAN,CCadmin,CCadmin -clean,true)
CONFIGDIR=`pwd`
CONFIGLEAF=`basename $CONFIGDIR`
AC_SUBST(CONFIGDIR)
AC_SUBST(CONFIGLEAF)
if test x$NEEDG2C = xyes; then
	LIBS="$LIBS -lg2c"
fi
AC_CONFIG_FILES([Makefile Makefile_template test/Makefile])
AC_OUTPUT
