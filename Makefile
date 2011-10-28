BAMTOOLS_ROOT=/scratch/ngsvin2/RNA-seq/ruping/Tools/bamtools/
CXX=g++
BAMFLAGS=-lbamtools
CXXFLAGS=-lz -static -Wall -O3
PREFIX=./
SRC=./src
LIB=./lib
BIN=/breakpointer/
SOURCE_BP=breakpointer.cpp
SOURCE_BM=breakmis.cpp
BP=breakpointer
BM=breakmis

all: breakpointer breakmis breakvali pipeline

.PHONY: all

breakpointer:
	@mkdir $(PREFIX)/$(BIN)
	@echo "* compiling" $(SOURCE_BP)
	@$(CXX) $(SRC)/$(SOURCE_BP) -o $(PREFIX)/$(BIN)/$(BP) $(BAMFLAGS) $(CXXFLAGS) -I $(BAMTOOLS_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ 

breakmis:
	@echo "* compiling" $(SOURCE_BM)
	@$(CXX) $(SRC)/$(SOURCE_BM) -o $(PREFIX)/$(BIN)/$(BM) $(BAMFLAGS) $(CXXFLAGS) -I $(BAMTOOLS_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/

breakvali:
	@echo "* copy breakvali script" 
	@cp $(SRC)/breakvali.pl $(PREFIX)/$(BIN)/

pipeline:
	@echo "* copy pipeline script"
	@cp $(LIB) $(PREFIX)/$(BIN)/ -r
	@cp $(SRC)/breakpointer_run.pl $(PREFIX)/$(BIN)/
	@echo "* done."

clean:
	@echo "Cleaning up everthing."
	@rm -rf $(PREFIX)/$(BIN)/

.PHONY: clean