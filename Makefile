CC = gcc
LD = gcc
CFLAGS = -O2 -fopenmp -fPIC 
LDFLAGS = -shared -fopenmp -llapacke -llapack -lblas -lm -lquadmath
OBJFILES = itersolve.o
SHARED = itersolve
SOURCE = itersolve.c

all: $(SHARED)

$(SHARED): $(OBJFILES)
	$(LD) -o $(SHARED) $(OBJFILES) $(LDFLAGS)

$(OBJFILES): $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -c -o $(OBJFILES) $(SOURCE) 

clean:
	rm -f $(OBJFILES) $(SHARED) *.jpg
