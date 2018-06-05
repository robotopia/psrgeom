CC      = gcc
LDLIBS = -lGL -lGLU -lglut -lm #-lnewuoa -lnmead

## Uncomment this to turn on compiler optimisation
OPTIM   = #-O2

## Uncomment these to turn on fsanitize option
#LDLIBS = -lasan -lm
#OPTIM   = -Wextra -O2 -fsanitize=address 

## Uncomment the following to turn on debugging options
#DEBUG   = -g -std=gnu99

CFLAGS  = -Wall -Wextra $(OPTIM) $(DEBUG) -march=native

SMALL_TARGETS = psr_fields \
				psr_lines \
				psr_lineofsight \
				psr_lines_sparks \
				psr_emit \
				psr_cost_function \
				psr_polangle \
				psr_fitwidth \
				psr_trajectory \
				psr_visiblepoints \
				psr_mcpa \
				psr_findcaps \
				psr_jacksonbeam \
				psr_lofl \
				psr_profile \
				psr_jacksonbeamsingle \
				psr_beam

TARGETS = psrgeom

TESTS = calc_fields_test \
		transform_new_xz_test \
		bstep_error_test \
		spark_test \
		spark_random_test

LIBRARY = psrgeom
LIBFILE = lib$(LIBRARY).a
HDRFILE = $(LIBRARY).h
OBJS = psr_angle.o \
	   point.o \
	   pulsar.o \
	   fields.o \
	   cylinder.o \
	   dipole.o \
	   emitloc.o \
	   psrio.o \
	   numrec.o \
	   jackson.o \
	   photon.o \
	   gamma.o
	   #fitwidth.o

all: $(TARGETS) $(LIBFILE) man-pages #$(SMALL_TARGETS) $(TESTS)

$(LIBFILE): $(OBJS)
	ar rcs $@ $^

$(TARGETS) $(TESTS): $(OBJS)

man-pages:
	$(MAKE) -C man

install:
	cp $(LIBFILE) /usr/local/lib
	cp $(HDRFILE) /usr/local/include
	cp $(TARGETS) /usr/local/bin
	mkdir -p /usr/local/man/man3 /usr/local/man/man7
	cp man/*.3.gz /usr/local/man/man3
	cp man/*.7.gz /usr/local/man/man7

clean-man-pages:
	$(RM) man/*.gz

clean:
	$(RM) *.o $(LIBFILE) $(TARGETS)
