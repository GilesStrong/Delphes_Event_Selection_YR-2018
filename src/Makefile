PREFIX=/usr/local
ROOTCFLAGS=$(shell $(PREFIX)/bin/root-config --cflags)
ROOTLIBS=$(shell $(PREFIX)/bin/root-config --libs)
DELPHES=$(HOME)/programs/Delphes-3.3.2
C++ = g++
C_FLAGS = -g
INCLUDES = -I$(PREFIX)/include -I$(DELPHES) -I$(DELPHES)/external -I./
LIBRARIES = -lGenVector -O2 -lboost_regex -lboost_filesystem -lboost_system -lDelphes -L$(DELPHES)
ROOT = $(ROOTCFLAGS) $(ROOTLIBS)
FLAGS =  $(INCLUDES) $(LIBRARIES) $(C_FLAGS) $(ROOT)
OBJS = main.o delphesReader.o  

delphes_event_selection: $(OBJS)
	$(C++) $(OBJS) -o delphes_event_selection $(FLAGS) 
	
main.o: main.cc main.hh finalstateCuts.hh delphesReader.h myMethods.hh
	$(C++) -c main.cc $(FLAGS)  
	
delphesReader.o: delphesReader.C delphesReader.h
	$(C++) -c delphesReader.C $(FLAGS)  

clean:                                              
	rm ./*.o ./*~ ./\#*