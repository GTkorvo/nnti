#
# Makefile
#

include ../../config.dir/package.mk

CFLAGS += -I.

# CFLAGS += -Wall

THE_LIB = ../../lib/libcmm.a

all: $(THE_LIB)
	cp $(INCS)    ../../include

SRCS = misc.c clock.c mem.c str.c math.c
INCS = cmm.h phys.h

OBJS=$(SRCS:.c=.o)

$(OBJS): $(INCS)

$(THE_LIB): $(OBJS)
	@(if [ ! -f ${THE_LIB} ] ; \
           then echo "Must recreate library." \
                echo $(AR) cr $(THE_LIB) $(OBJS); \
                $(AR) cr $(THE_LIB) $(OBJS);      \
           else echo "$(AR) r $(THE_LIB) $(OBJS)";  \
                $(AR) r $(THE_LIB) $(OBJS);       \
         fi;)

#
# Build test programs
#
TESTS = test1 test2 test3
TEST_LIBS = $(BLAS_LIBS) $(FLIBS) $(LIBS)

test:
	@echo "Build cmmlib/base test programs"
	make $(TESTS)

#
# Build examples
#
examples:
	@echo "Build cmmlib/base examples"

clean:
	rm -rf core *.o *.a $(TESTS)

test1: test1.o $(THE_LIB)
	$(CC) $(CFLAGS) -o test1 test1.o $(THE_LIB) $(TEST_LIBS)

test2: test2.o $(THE_LIB)
	$(CC) $(CFLAGS) -o test2 test2.o $(THE_LIB) $(TEST_LIBS)

test3: test3.o $(THE_LIB)
	$(CC) $(CFLAGS) -o test3 test3.o $(THE_LIB) $(TEST_LIBS)

