AM_CFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm
AM_CXXFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm 

bin_PROGRAMS = bwtformatdb reverse batmeth filter split
bwtformatdb_SOURCES= bwtformatdb.c BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h
batmeth_SOURCES=batmeth.cpp common.h map.cpp map.h commandline.cpp commandline.h batlib.cpp batlib.h BWT.c MiscUtilities.c MemManager.c TextConverter.c r250.c iniparser.c inistrlib.c dictionary.c DNACount.c HSP.c HSPstatistic.c karlin.c Timing.c Socket.c  
filter_SOURCES= filter.cpp 
split_SOURCES= split.cpp 
reverse_SOURCES = reverse.cpp 
copy:
	cp bwtformatdb ../bin
	cp reverse ../bin
	cp batmeth ../bin
	cp filter ../bin
	cp split ../bin
