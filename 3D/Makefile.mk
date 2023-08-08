CC			=g++-13
CFLAGS		=-c -Wall -O4 -std=c++17 -fopenmp -I/usr/local/include/eigen3 -I/usr/local/include/
LDFLAGS		=-fopenmp -std=c++17 -I../#-L/usr/local/opt/icu4c/lib #-Ofast
SOURCES		=./main.cpp
OBJECTS		=$(SOURCES:.cpp=.o)

DTYPE_FLAG  = -DUSE_real -DUSE_nHODLRdD # DUSE_nHODLRdD (nested) or -DUSE_snHODLRdD (semi-nested) or -DUSE_HODLRdD (non-nested)
#DTYPE_FLAG  = -DUSE_Helm # Helmholtz

EXECUTABLE	=./test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(DTYPE_FLAG) $< -o $@

clean:
	rm test *.o