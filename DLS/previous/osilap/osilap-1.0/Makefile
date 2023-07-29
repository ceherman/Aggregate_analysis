CC      = gcc
FC      = gfortran
CFLAGS  = -O2 -Wall
LIBS    = -llapack -lblas -lm
LIBS    = libopenblas-mac.a -lpthread -static-libgfortran -static-libgcc

LDFLAGS = 
OBJS    = ilt.o nnls.o command.o output.o input.o BRDnewton_modify.o alpha_loop.o
PROGRAM = osilap

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:
	rm -f *.o *~ $(PROGRAM)

install: $(PROGRAM)
	install -s $(PROGRAM) $(DEST)
