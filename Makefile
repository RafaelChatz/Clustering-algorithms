RM := rm -f
CC := g++ -std=c++11
CFLAGS := -c


LSH_SOURCE := cluster.cpp
SOURCES := $(wildcard *.cpp)
COMMON_SOURCES := $(filter-out $(LSH_SOURCE) $(CUBE_SOURCE),$(wildcard *.cpp))
OBJECTS := $(SOURCES:%.cpp=%.o)
COMMON_OBJECTS := $(COMMON_SOURCES:%.cpp=%.o)
FILES_TO_CLEAN := $(OBJECTS) 

all : cluster 

clean :
	$(RM) $(FILES_TO_CLEAN) lsh cube


$(OBJECTS) : %.o : %.cpp
	$(CC) $(CFLAGS) -o $@ $<

cluster : $(COMMON_OBJECTS) cluster.o
	$(CC) -o $@ $^ $(LDLIBS)

