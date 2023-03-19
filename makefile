TARGET = exe
CC = g++
LIBS = -lm -lboost_system -lboost_filesystem -lboost_iostreams -lpthread -lquadmath
HEAD = ./include
SRCS = ./src
INCDIRS = -I$(HEAD)
CFLAGS = -Wall -O3 -fPIC -std=c++20 -mlong-double-128 -fext-numeric-literals -ffast-math -funroll-loops $(INCDIRS)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp)) $(MOC_SOURCES:.cpp=.o)

#all: clean $(TARGET)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f $(TARGET)
	-rm -f $(SRCS)/*.o
	-rm -f $(HEAD)/*.o
