#
# Makefile for 3d geometry routines
#

include ../../config.dir/package.mk

CFLAGS += -I. -I../../include

# CFLAGS += -Wall

THE_LIB = ../../lib/libcmm.a

all: $(THE_LIB)
	cp $(INCS)    ../../include

SRCS = geom.c
INCS = geom.h

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
TESTS = 
TEST_LIBS = $(BLAS_LIBS) $(FLIBS) $(LIBS)

test:
	@echo "Build cmmlib geom3 test programs"
	make $(TESTS)

#trand: trand.o $(THE_LIB)
#	$(CC) $(CFLAGS) -o trand trand.o $(THE_LIB) $(TEST_LIBS)

#
# Build examples
#
examples:
	@echo "Build cmmlib/base examples"

