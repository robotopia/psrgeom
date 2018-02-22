CC      = gcc
LDLIBS = -lm

## Uncomment this to turn on compiler optimisation
OPTIM   = -O2

## Uncomment these to turn on fsanitize option
#LDLIBS = -lasan -lm
#OPTIM   = -Wextra -O2 -fsanitize=address 

## Uncomment the following to turn on debugging options
#DEBUG   = -g -std=gnu99

CFLAGS  = -Wall -Wextra $(OPTIM) $(DEBUG)

LIBRARY = psrgeom
LIBFILE = lib$(LIBRARY).a
HDRFILE = $(LIBRARY).h
OBJS = psr_angle.o \
	   point.o \
	   lightcyl.o \
	   pulsar.o \
	   bfield.o \
	   cylinder.o \
	   dipole.o \
	   vfield.o

all: $(LIBFILE) man-pages

$(LIBFILE): $(OBJS)
	ar rcs $@ $^

man-pages:
	$(MAKE) -C man

install:
	cp $(LIBFILE) /usr/local/lib
	cp $(HDRFILE) /usr/local/include
	cp man/*.gz /usr/local/man/man3

clean-man-pages:
	$(RM) man/*.gz

clean:
	$(RM) *.o $(LIBFILE)
