TARGET = exe
CC = g++
LIBS = -lm -lboost_system -lboost_filesystem -lboost_iostreams -lpthread
HEAD = ./include
SRCS = ./src
INCDIRS = -I$(HEAD)
CFLAGS = -Wall -O3 -fPIC -std=c++20 -fext-numeric-literals -ffast-math -funroll-loops $(INCDIRS)
.PHONY: clean

DEPS = $(wildcard $(HEAD)/*.hpp)
#MOCS = $(shell grep -l Q_OBJECT $(DEPS))
#MOC_SOURCES = $(MOCS:.hpp=.moc.cpp)
OBJS = $(patsubst %.cpp, %.o, $(wildcard $(SRCS)/*.cpp)) $(MOC_SOURCES:.cpp=.o)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

#%.moc.cpp: %.hpp
#	$(MOC) $(INCDIRS) $< -o $@

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

clean:
	-rm -f $(TARGET)
	-rm -f $(SRCS)/*.o
	-rm -f $(HEAD)/*.o
