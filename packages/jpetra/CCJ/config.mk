# When using javac or jikes, set these values correctly and use "make"
# or "make java". (Also make sure that your CLASSPATH is set correctly).     

default:	java
# default:	manta

# --------------- Configure must decide sensible values -----------------
JAVA_ROOT	= /System/Library/Frameworks/JavaVM.framework
JAVAC_COMMAND	= /usr/bin/javac
CCJ_ROOT	= /Users/jc/Desktop/research/sharedFileSystem/CCJ-0.1
MANTA_ROOT	= /usr/local/VU/manta
RMIC_LOWER	= 

# --------------- End of configurable part ------------------------------

JAVA_CFLAGS	= -classpath .:$(CCJ_ROOT):$(CCJ_ROOT)/jg-bench:$(JAVA_ROOT)/jre/lib/rt.jar
RMIC		= /usr/bin/rmic

CLASS_FILES	= $(SRC:%.java=%.class)
MAIN_CLASS	= $(MAIN:%.java=%.class)

MANTAC_COMMAND	= $(MANTA_ROOT)/mantac
MANTA_CFLAGS	=
MANTA_CFLAGS	+= $(JAVA_CFLAGS)
MANTA_CFLAGS	+= -keep
# MANTA_CFLAGS	+= -no_opt_packagers
# MANTA_CFLAGS	+= -g2

MANTA_LIB	= $(CCJ_ROOT)/libccj.a

JAVA_LDFLAGS	= -myrinet -no-shared-link -v
JAVA_LDLIBS	= $(MANTA_LIB)

MANTA_RMIC	= /usr/local/VU/jdk1.1/bin/rmic

MAIN_EXEC	= $(MAIN:%.java=%)

ifneq ($(MAIN),)
	MAIN_PACKAGER_LINE	:= $(shell grep -w package $(MAIN))
	ifeq ($(MAIN_PACKAGER_LINE),)
		MAIN_PACKAGE	:= ""
	else
		MAIN_PACKAGE_SC		:= $(filter-out package, $(MAIN_PACKAGER_LINE))
		MAIN_PACKAGE_DOT	:= $(subst ;,, $(MAIN_PACKAGE_SC))
		ifeq ($(MAIN_PACKAGE_DOT), ".")
			MAIN_PACKAGE	:= ""
		else
			MAIN_PACKAGE	:= $(MAIN_PACKAGE_DOT).
		endif
	endif
endif


O_FILES		= $(SRC:%.java=%.o)

CCJ_FILES= $(CCJ_ROOT)/CCJ/*.o
JGF_FILES= $(CCJ_ROOT)/jg-bench/jgfutil/*.o

MPIJAVA_FILES=/home/jason/CPE-01/mpiJava/mpiJava/src/Java/mpi/*.o /home/jason/CPE-01/mpiJava/mpiJava/lib/libmpijava.a /usr/local/VU/mpi/lib/i386/ch_panda4/optimized/libmpi.a

TRASH_FILES	=
TRASH_FILES	+= *core
TRASH_FILES	+= packagers_*.c packagers_*.h stub_*.h
TRASH_FILES	+= SUN_RMI_MASTER
TRASH_FILES	+= *.o *.jac *.s *.lasm
TRASH_FILES	+= *~
TRASH_FILES	+= $(OUT)
TRASH_FILES	+= *.class
TRASH_FILES	+= *_Skel.java
TRASH_FILES	+= *_Stub.java
TRASH_FILES	+= gmon.out
TRASH_FILES	+= a.out
TRASH_FILES	+= -r LASM_CACHE

$(CLASS_FILES) $(MAIN_CLASS):	$(SRC) $(MAIN)
	$(JAVAC_COMMAND) $(JAVA_CFLAGS) $(SRC) $(MAIN)

java:	$(CLASS_FILES)

#conversion rules

.SUFFIXES:
.SUFFIXES: .o .java .satin .class

.java.o:
	$(MANTAC_COMMAND) $(MANTA_CFLAGS) -c $*.java

# .java.class:
# 	$(JAVAC_COMMAND) $(JAVA_CFLAGS) $*.java

.java:
	$(MANTAC_COMMAND) $(MANTA_CFLAGS) $(JAVA_LDFLAGS) $(MAIN) $(O_FILES) $(JAVA_LDLIBS) -o $*

.PHONY: run
run:
	@echo $(JAVA_ROOT)/bin/java $(JAVA_CFLAGS) -DCCJ_HOSTS=\"\$$HOSTS\" $(MAIN:%.java=%) 'args'

.PHONY: prun
prun:
	@echo java: prun -subst -no-panda $(JAVA_ROOT)/bin/java "<nhosts>" $(JAVA_RFLAGS) $(JAVA_CFLAGS) -DCCJ_HOSTS=%h "$(MAIN_PACKAGE)$(MAIN:%.java=%)" 'args'
	@echo manta: prun -subst $(MAIN:%.java=%) "<nhosts>" -DCCJ_HOSTS=%h 'args'

$(MAIN_EXEC):	$(MAIN) $(O_FILES) $(MANTA_LIB)

.PHONY: manta
manta:	$(MAIN_EXEC)
