CC      = gcc
LDLIBS = -lm

## Uncomment this to turn on compiler optimisation
OPTIM   = -O3

## Uncomment these to turn on fsanitize option
#LDLIBS = -lasan -lm
#OPTIM   = -Wextra -O3 -fsanitize=address 

## Uncomment the following to turn on debugging options
#DEBUG   = -g -std=gnu99

CFLAGS  = -Wall -Wextra $(OPTIM) $(DEBUG)

LIBRARY = psrgeom
LIBFILE = lib$(LIBRARY).a
HDRFILE = $(LIBRARY).h
OBJS = angle.o \
	   point.o \
	   lightcyl.o \
	   pulsar.o \
	   bfield.o \
	   cylinder.o

all: $(LIBFILE)

$(LIBFILE): $(OBJS)
	ar rcs $@ $^

install:
	cp $(LIBFILE) /usr/local/lib
	cp $(HDRFILE) /usr/local/include

clean:
	$(RM) *.o $(LIBFILE)
