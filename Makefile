CC=gcc
CXX=g++

CXXFLAGS=-Wall -std=c++11 -fPIC -g
CFLAGS=-Wall -std=gnu99 -fPIC -g 
LIBS=-lstdc++

INCS=-I.
# CXXSRCS=$(wildcard ./*.cc)
SRCS=$(wildcard ./*.c)
OBJS=$(SRCS:.c=.o)
# CXXOBJS=$(SRCS:.cc=.o)

BIN=app


all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $< -o $@ $(CFLAGS)

%.o: %.cc
	$(CXX) -c $^ -o $@ $(CXXFLAGS)

%.o: %.c
	$(CC) -c $^ -o $@ $(INCS) $(CFLAGS)

clean:
	rm $(BIN) $(OBJS)
