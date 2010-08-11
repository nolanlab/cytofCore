AC_INIT([cytofCore], 0.10)

# Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
    if test -z "${R_HOME}"; then
	echo "could not determine R_HOME"
	exit 1
    fi
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
AC_PROG_CXX	
AC_PROG_CXXCPP

AC_SYS_LARGEFILE

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
