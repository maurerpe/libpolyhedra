#############################################################################
# Copyright (C) 2021 Paul Maurer
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
#  
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#############################################################################

AC_PREREQ([2.69])
AC_INIT([libpolyhedra], [0.1], [junk4me46806@yahoo.com])
AM_INIT_AUTOMAKE([foreign subdir-objects])
LT_INIT
AC_CONFIG_SRCDIR([include/libpolyhedra.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
AC_SEARCH_LIBS([sqrtf], [m], [], [AC_MSG_ERROR([Missing required function sqrtf])])
AC_SEARCH_LIBS([pthread_mutex_lock], [pthread], [AC_DEFINE([HAVE_PTHREADS], [1], [Have pthreads support])])
AC_SEARCH_LIBS([BCryptGenRandom], [Bcrypt], [AC_DEFINE([HAVE_BCRYPTGENRANDOM], [1], [Have BCryptGenRandom support (windows)])])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h], [], [AC_MSG_ERROR([Missing required header])])
AC_CHECK_HEADERS([unistd.h pthread.h Windows.h Bcrypt.h Process.h locale.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t], [], [AC_MSG_ERROR([Missing required type ptrdiff_t])])

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset strcasecmp strdup strtoull], [], [AC_MSG_ERROR([Missing required function])])
AC_CHECK_FUNCS([getentropy CreateMutexA setlocale])

build_prog=true
AC_CHECK_FUNCS([getopt], [], [AC_MSG_WARN([Missing func getopt, command line utility will not be built.])
build_prog=false])
AM_CONDITIONAL([BUILD_PROG], [test "x$build_prog" = "xtrue"])

AC_CONFIG_FILES([Makefile
                 include/Makefile
                 lib/Makefile
                 src/Makefile])
AC_OUTPUT
