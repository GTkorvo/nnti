SET(TPL_ENABLE_JDK     ON      CACHE BOOL "")
SET(TPL_ENABLE_OpenSSL ON      CACHE BOOL "")
SET(TPL_ENABLE_Zlib    ON      CACHE BOOL "")
SET(TPL_ENABLE_TCL     ON      CACHE BOOL "")
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "")

# We don't have the MATLAB TPL for SEACAS
SET(TPL_ENABLE_MATLAB OFF CACHE BOOL "" FORCE)

SET(ANT_PATH         /usr/bin                                                     CACHE FILEPATH "")
SET(TCLSH_PATH       /usr/bin                                                     CACHE FILEPATH "")
SET(TCL_LIBRARY_NAMES "tcl8.5"                                                    CACHE STRING   "")
SET(JDK_INCLUDE_DIRS "/usr/lib/jvm/java/include;/usr/lib/jvm/java/include/linux"  CACHE FILEPATH "")
SET(JDK_LIBRARY_DIRS /usr/lib/jvm/java/jre/lib/amd64/server                       CACHE FILEPATH "")
SET(Netcdf_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/include CACHE FILEPATH "")
SET(Netcdf_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/netcdf-4.1.1/lib CACHE FILEPATH "")
SET(SILO_INCLUDE_DIRS   /opt/gcc-4.5.1/tpls/silo/include CACHE FILEPATH "")
SET(SILO_LIBRARY_DIRS   /opt/gcc-4.5.1/tpls/silo/lib     CACHE FILEPATH "")
SET(TBB_INCLUDE_DIRS    /opt/intel/Compiler/composerxe-2011.4.191/tbb/include                                      CACHE FILEPATH "")
SET(TBB_LIBRARY_DIRS    /opt/intel/Compiler/composerxe-2011.4.191/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21  CACHE FILEPATH "")
SET(LIBXML2_INCLUDE_DIRS  /usr/include/libxml2       CACHE FILEPATH "")
SET(LIBXML2_LIBRARY_DIRS  /usr/lib64                  CACHE FILEPATH "")
