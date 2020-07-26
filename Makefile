CC=gcc
CXX=g++

CXXFLAGS=-Wall -std=c++11 -fPIC
LIBS=-lstdc++

INCS=-I.
SRCS=$(wildcard ./*.cc)
OBJS=$(SRCS:.cc=.o)

BIN=app

all: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $< -o $@ $(CXXFLAGS)

%.o: %.cc
	$(CXX) -c $^ -o $@ $(CXXFLAGS)

clean:
	rm $(BIN) $(OBJS)
