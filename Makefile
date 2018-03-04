
CC=g++ 

CFLAGS=-Wall -std=c++11 
LFLAGS=-L/usr/lib -lm 

EXEC=MeshReconstructor
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

$(EXEC) : $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS)

%.o : %.cpp
	$(CC) -c $^ -o $@ $(CFLAGS)

clean :
	/bin/rm $(EXEC) *.o
