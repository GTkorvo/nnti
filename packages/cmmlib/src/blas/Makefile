#
# Makefile
#

include ../../config.dir/package.mk

#
# BLAS_INCLUDE and BLAS_LIBS should be defined in cmmlib.mk, set up by configure script.
#

CFLAGS += -I. -I../../include $(BLAS_INCLUDE) $(FORTRAN_MANGLE)

 CFLAGS += -Wall

THE_LIB = ../../lib/libcmm.a

all: $(THE_LIB)
	cp $(INCS)    ../../include

SRCS = blas1.c blas2.c blas3.c myzgemm.c
INCS = cmmblas.h

OBJS=$(SRCS:.c=.o)

$(OBJS): $(INCS)

#
# This is a little complicated because we want to give error message
#
$(THE_LIB): $(OBJS)
	@(if [ ! -f ${THE_LIB} ] ; \
           then echo "Must recreate library. make ../base first" \
                exit 1; \
           else echo "$(AR) r $(THE_LIB) $(OBJS)";  \
                $(AR) r $(THE_LIB) $(OBJS);       \
         fi;)

doc:
	latex blas
	dvips -t letter -o blas.ps blas

clean:
	rm -rf core *.o *.a $(TESTS)
	rm -f *.log *.aux *.dvi *.ps

scrub:
	make clean

#
# Build test programs
#
TESTS = tcheck tperf
TEST_LIBS = $(BLAS_LIBS) $(FLIBS) $(LIBS)

test:
	@echo "Build cmmlib/blas test programs"
	make $(TESTS)

tcheck: tcheck.o $(THE_LIB)
	$(CC) $(CFLAGS) -o tcheck tcheck.o $(THE_LIB) $(TEST_LIBS)

tperf: tperf.o $(THE_LIB)
	$(CC) $(CFLAGS) -o tperf tperf.o $(THE_LIB) $(TEST_LIBS)


#
# Build examples
#
examples:
	@echo "Build cmmlib/base examples"

