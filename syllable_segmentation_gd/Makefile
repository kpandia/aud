SRCDIR = ./src
DSPLIB = ./front-end-dsp/src
DSP_DIR_FE = ./front-end-dsp/src/fe/
OBJDIR = ./obj
BINDIR = ./bin
NISTINC = ./nist/include/
NISTLIB = ./nist/lib/
CFLAGS = -g -c -Wall  
CC = gcc

all :   $(BINDIR)/SyllableBoundaries

$(BINDIR)/SyllableBoundaries : $(OBJDIR)/SyllableBoundaries.o \
	$(OBJDIR)/CommonFunctions.o \
	$(OBJDIR)/DspLibrary.o \
	$(OBJDIR)/SphereInterface.o \
	$(OBJDIR)/BatchProcessWaveform.o \
	$(OBJDIR)/InitAsdf.o \
        $(OBJDIR)/FeatureMappingToFunctions.o\
        $(OBJDIR)/QuickSort.o\
	$(OBJDIR)/GMM.o \
	$(OBJDIR)/VQ.o \
        $(OBJDIR)/PthLibrary.o \
        $(OBJDIR)/FmtsLibrary.o \
        $(OBJDIR)/HashTable.o 
	$(CC) -g -o $(BINDIR)/SyllableBoundaries  \
	-I$(DSPLIB)\
	-I$(NISTINC) \
	-L$(NISTLIB) \
	$(OBJDIR)/CommonFunctions.o \
	$(OBJDIR)/SyllableBoundaries.o \
	$(OBJDIR)/DspLibrary.o \
	$(OBJDIR)/SphereInterface.o \
	$(OBJDIR)/BatchProcessWaveform.o \
	$(OBJDIR)/InitAsdf.o \
        $(OBJDIR)/PthLibrary.o \
        $(OBJDIR)/FmtsLibrary.o \
        $(OBJDIR)/FeatureMappingToFunctions.o\
        $(OBJDIR)/HashTable.o \
        $(OBJDIR)/QuickSort.o \
	$(OBJDIR)/GMM.o \
	$(OBJDIR)/VQ.o \
	-lsp \
	-lutil \
	-lm


$(OBJDIR)/SyllableBoundaries.o: $(SRCDIR)/SyllableBoundaries.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/SyllableBoundaries.o $(SRCDIR)/SyllableBoundaries.c -I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/CommonFunctions.o : $(SRCDIR)/CommonFunctions.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/CommonFunctions.o $(SRCDIR)/CommonFunctions.c -I$(NISTINC) -I$(DSPLIB) -I$(DSP_DIR_FE)

$(OBJDIR)/DspLibrary.o : $(DSPLIB)/DspLibrary.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/DspLibrary.o $(DSPLIB)/DspLibrary.c \
	-I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/GMM.o : $(DSPLIB)/GMM.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/GMM.o $(DSPLIB)/GMM.c -I$(NISTINC) -I$(DSPLIB) -I$(DSP_DIR_FE)

$(OBJDIR)/VQ.o : $(DSPLIB)/VQ.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/VQ.o $(DSPLIB)/VQ.c -I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/BatchProcessWaveform.o : $(DSPLIB)/BatchProcessWaveform.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/BatchProcessWaveform.o $(DSPLIB)/BatchProcessWaveform.c -I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/PthLibrary.o : $(DSPLIB)/PthLibrary.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/PthLibrary.o $(DSPLIB)/PthLibrary.c -I$(NISTINC) -I$(DSPLIB)/fe

$(OBJDIR)/FmtsLibrary.o : $(DSPLIB)/FmtsLibrary.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/FmtsLibrary.o $(DSPLIB)/FmtsLibrary.c -I$(NISTINC)

$(OBJDIR)/SphereInterface.o : $(DSPLIB)/SphereInterface.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/SphereInterface.o $(DSPLIB)/SphereInterface.c -I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/InitAsdf.o : $(DSPLIB)/InitAsdf.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/InitAsdf.o $(DSPLIB)/InitAsdf.c -I$(NISTINC) -I$(DSPLIB)

$(OBJDIR)/HashTable.o : $(DSPLIB)/HashTable.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/HashTable.o $(DSPLIB)/HashTable.c \
	-I$(NISTINC) \
	-I$(DSPLIB)/fe 

$(OBJDIR)/QuickSort.o : $(DSPLIB)/QuickSort.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/QuickSort.o $(DSPLIB)/QuickSort.c \
	-I$(NISTINC) \
	-I$(DSPLIB)/fe 

$(OBJDIR)/FeatureMappingToFunctions.o : $(DSPLIB)/FeatureMappingToFunctions.c 
	$(CC) $(CFLAGS) -o $(OBJDIR)/FeatureMappingToFunctions.o $(DSPLIB)/FeatureMappingToFunctions.c \
	-I$(NISTINC) \
	-I$(DSPLIB)/fe 

clean: 
	rm -fr $(BINDIR)/SyllableBoundaries $(OBJDIR)/DspLibrary.o $(OBJDIR)/BatchProcessWaveform.o $(OBJDIR)/InitAsdf.o ./obj/HashTable.o $(OBJDIR)/FeatureMappingToFunctions $(OBJDIR)/SphereInterface.o $(OBJDIR)/SyllableBoundaries.o


