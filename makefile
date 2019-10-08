CC = mpic++ -g -gdwarf-2
CFLAGS = -std=c++11 -Wall
DFLAGS = -std=c++11 -Wall -ggdb
 
all: CodedTeraSort Splitter InputPlacement

clean:
	rm -f *.o
	rm -f CodedTeraSort Splitter Codegen BroadcastTest InputPlacement

cleanclean: clean
	rm -f /home/ubuntu/Input/*_*
	rm -f /home/ubuntu/Output/*_*
	rm -f /home/ubuntu/Tmp/*
	rm -f *.*~
	rm -f *~



CodedTeraSort: CodedMain.o CodedMaster.o CodedWorker.o Trie.o Utility.o PartitionSampling.o CodeGeneration.o
	$(CC) $(CFLAGS) -o CodedTeraSort CodedMain.o CodedMaster.o CodedWorker.o Trie.o Utility.o PartitionSampling.o CodeGeneration.o

Splitter: InputSplitter.o CodedConfiguration.h
	$(CC) $(CFLAGS) -o Splitter Splitter.cc InputSplitter.o

InputPlacement: InputPlacement.cc CodeGeneration.o CodedConfiguration.h
	$(CC) $(CFLAGS) -o InputPlacement InputPlacement.cc CodeGeneration.o
	

Trie.o: Trie.cc Trie.h
	$(CC) $(CFLAGS) -c Trie.cc

PartitionSampling.o: PartitionSampling.cc PartitionSampling.h CodedConfiguration.h
	$(CC) $(CFLAGS) -c PartitionSampling.cc

Utility.o: Utility.cc Utility.h
	$(CC) $(CFLAGS) -c Utility.cc

InputSplitter.o: InputSplitter.cc InputSplitter.h CodedConfiguration.h
	$(CC) $(CFLAGS) -c InputSplitter.cc

CodeGeneration.o: CodeGeneration.cc CodeGeneration.h
	$(CC) $(CFLAGS) -c CodeGeneration.cc



CodedMain.o: CodedMain.cc CodedConfiguration.h
	$(CC) $(CFLAGS) -c CodedMain.cc

CodedMaster.o: CodedMaster.cc CodedMaster.h CodedConfiguration.h
	$(CC) $(CFLAGS) -c CodedMaster.cc

CodedWorker.o: CodedWorker.cc CodedWorker.h CodedConfiguration.h
	$(CC) $(CFLAGS) -c CodedWorker.cc
