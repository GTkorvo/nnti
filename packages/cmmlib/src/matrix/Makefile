#
# Makefile
#

include ../../config.dir/package.mk

CFLAGS += -I. -I../../include

CFLAGS += -Wall

THE_LIB = ../../lib/libcmm.a

all: $(THE_LIB)
	cp $(INCS)    ../../include

SRCS = c8chol.c c8tsolve.c                 \
  c8band.c  c8bandgersh.c  c8bandinvit.c   \
  c8bandoracle.c  c8bandqr.c c8bandeigen.c \
  c8house.c c8qr.c                         \
  c8condense.c c8eigen.c

INCS = cmmatrix.h cmmband.h

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
	latex matrix
	dvips -t letter -o matrix.ps matrix

clean:
	rm -rf core *.o *.a $(TESTS)
	rm -f *.log *.aux *.dvi *.ps

scrub:
	make clean

#
# Build test programs
#
TESTS = tchol tqr teigen
TEST_LIBS = $(BLAS_LIBS) $(FLIBS) $(LIBS)

test:
	@echo "Build cmmlib/matrix test programs"
	make $(TESTS)

tchol: tchol.o $(THE_LIB)
	$(CC) $(CFLAGS) -o tchol tchol.o $(THE_LIB) $(TEST_LIBS)

tqr: tqr.o $(THE_LIB)
	$(CC) $(CFLAGS) -o tqr tqr.o $(THE_LIB) $(TEST_LIBS)

teigen: teigen.o $(THE_LIB)
	$(CC) $(CFLAGS) -o teigen teigen.o $(THE_LIB) $(TEST_LIBS)


#
# Build examples
#
examples:
	@echo "Build cmmlib/matrix examples"

