#
# Copyright 2001 Vrije Universiteit, The Netherlands.
# For full copyright and restrictions on use see the file COPYRIGHT in the
# top level of the CCJ distribution.
#

SRC =
SRC += CCJException.java
SRC += ColMember.java
SRC += ColGroup.java
SRC += ColGroupCentral.java
SRC += ColGroupCentralInterface.java
SRC += ColGroupMaster.java
SRC += ColMemberInterface.java
SRC += Partitionable.java
SRC += GroupAlreadyActiveException.java
SRC += Message.java
SRC += MyEntry.java
SRC += MyList.java
SRC += MyListInterface.java
SRC += NoSuchGroupException.java
SRC += NoSuchMemberException.java
SRC += RankAlreadyInUseException.java
SRC += Reducible.java
SRC += SendThread.java
SRC += ThreadPool.java
SRC += CcjHosts.java


MAIN =
OUT  =

STUB_SRC	=
STUB_SRC	+= ColGroupCentral
STUB_SRC	+= MyList

STUB_FILES	= $(STUB_SRC:%=%_Stub.java)
SKEL_FILES	= $(STUB_SRC:%=%_Skel.java)
STUB_OBJ	= $(STUB_SRC:%=%_Stub.o) $(STUB_SRC:%=%_Skel.o)

include ../config.mk

ifeq ($(RMIC_LOWER), 2)
%_Stub.java:	%.java
	(cd ..; $(RMIC) $(JAVA_CFLAGS) CCJ.$* -keepgenerated)
else
%_Stub.java:	%.java
	$(RMIC) $(JAVA_CFLAGS) CCJ.$* -keepgenerated
endif

.PHONY:		java
java: classes stubs

.PHONY:		classes
classes: $(CLASS_FILES)
	$(JAVAC_COMMAND) $(JAVA_CFLAGS) $(SRC)

.PHONY:		stubs
stubs: $(STUB_FILES)

.PHONY:		rmic_1_2_stubs
rmic_1_2_stubs:
	$(MAKE) $(MFLAGS) RMIC_LOWER=2 stubs

.PHONY:		manta_stubs
manta_stubs:
	$(MAKE) $(MFLAGS) RMIC_LOWER=1 RMIC=$(MANTA_RMIC) stubs
	$(MAKE) $(MFLAGS) $(STUB_OBJ)

.PHONY:		manta
manta: objs manta_stubs

.PHONY:		objs
objs: $(O_FILES) manta_stubs
	gcc -shared  -Wl,-soname,libccj.so -o ../libccj.so $(O_FILES) $(STUB_OBJ)
	ar urs $(MANTA_LIB) $(O_FILES) $(STUB_OBJ)

TRASH_FILES	+= $(STUB_FILES)
TRASH_FILES	+= $(SKEL_FILES)
TRASH_FILES	+= *.so

.PHONY: clean
clean:
	-rm -f $(TRASH_FILES)
