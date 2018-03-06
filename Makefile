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

TARGETS = psr_fields

LIBRARY = psrgeom
LIBFILE = lib$(LIBRARY).a
HDRFILE = $(LIBRARY).h
OBJS = psr_angle.o \
	   point.o \
	   lightcyl.o \
	   pulsar.o \
	   fields.o \
	   cylinder.o \
	   dipole.o \
	   psrio.o

all: $(TARGETS) $(LIBFILE) man-pages

$(LIBFILE): $(OBJS)
	ar rcs $@ $^

$(TARGETS): $(OBJS)

man-pages:
	$(MAKE) -C man

install:
	cp $(LIBFILE) /usr/local/lib
	cp $(HDRFILE) /usr/local/include
	cp $(TARGETS) /usr/local/bin
	cp man/*.3.gz /usr/local/man/man3
	cp man/*.7.gz /usr/local/man/man7

clean-man-pages:
	$(RM) man/*.gz

clean:
	$(RM) *.o $(LIBFILE) $(TARGETS)
