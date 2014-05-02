BAMTOOLS_ROOT=/ifs/home/c2b2/ac_lab/rs3412/tools/bamtools/
BOOST_ROOT=/ifs/home/c2b2/ac_lab/rs3412/tools/boost_1_54_0/
ZLIB_ROOT=/ifs/home/c2b2/ac_lab/rs3412/tools/zlib-1.2.8/
CXX=g++
BAMFLAGS=-lbamtools
CXXFLAGS=-lz -Wall
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
	@$(CXX) $(SRC)/$(SOURCE_BP) -o $(PREFIX)/$(BIN)/$(BP) $(BAMFLAGS) $(CXXFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/ -Wl,-rpath,$(BAMTOOLS_ROOT)/lib/:$(BOOST_ROOT)/lib/

breakmis:
	@echo "* compiling" $(SOURCE_BM)
	@$(CXX) $(SRC)/$(SOURCE_BM) -o $(PREFIX)/$(BIN)/$(BM) $(BAMFLAGS) $(CXXFLAGS) -I $(BAMTOOLS_ROOT)/include/ -I $(ZLIB_ROOT)/include/ -I $(BOOST_ROOT)/include/ -L $(BAMTOOLS_ROOT)/lib/ -L $(ZLIB_ROOT)/lib/ -L $(BOOST_ROOT)/lib/ -Wl,-rpath,$(BAMTOOLS_ROOT)/lib/:$(BOOST_ROOT)/lib/

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