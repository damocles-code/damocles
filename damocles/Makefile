FC=gfortran
LD=gfortran
VERSION := $(shell git describe --always --tags --dirty || echo "2.0")
FFLAGS += -cpp -fPIC -fbounds-check -ffree-line-length-0  -fopenmp  -Jsource/ -DPREFIX=\"${PREFIX}\" -DVERSION=\"${VERSION}\" #-L/usr/lib #-L/usr/local/lib -L/Library/Frameworks/ -L/System/Library/Frameworks/ -L/usr/local/opt

# set prefix depending on OS
OS := $(shell uname)
ifeq ($(OS),Darwin)
  PREFIX=/usr/local
  FFLAGS += -L/usr/lib
else
  PREFIX=/usr
endif

ifeq ($(CO),debug)
  FFLAGS += -g -fbounds-check -Wall -Wuninitialized -ffpe-trap=zero,overflow,invalid,underflow,denormal
else
  FFLAGS += -O3
endif

.PHONY: all clean new install 

all:	damocles

new:    clean all

source/%.o: source/%.f90
	$(FC) $(FFLAGS) $^ -c -o $@

damocles: source/globals.o source/class_line.o source/class_geometry.o source/class_dust.o source/class_obs_data.o \
source/class_freq_grid.o source/rnglib.o source/class_grid.o source/electron_scattering.o source/input.o source/init_random_seed.o  source/random_routines.o source/vector_functions.o \
source/write_out.o source/class_packet.o source/BHmie.o  source/radiative_transfer.o \
source/model_comparison.o source/driver.o source/damocles_wrap.o source/damocles.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

clean:
	rm -f damocles source/*.o source/*.mod

install: damocles
	test -e ${DESTDIR}${PREFIX}/share/damocles || mkdir -p ${DESTDIR}${PREFIX}/share/damocles
	test -e ${DESTDIR}${PREFIX}/bin || mkdir -p ${DESTDIR}${PREFIX}/bin
	test -e ./output || mkdir -p ./output
	echo ${DESTDIR}${PREFIX}
	cp -R dustData ${DESTDIR}${PREFIX}/share/damocles/
	install damocles ${DESTDIR}${PREFIX}/bin

uninstall:
	rm -f ${DESTDIR}${PREFIX}/bin/damocles
	rm -rf ${DESTDIR}${PREFIX}/share/damocles
