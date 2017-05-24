FC=gfortran
LD=gfortran
FFLAGS=-fbounds-check -ffree-line-length-0  -Wunused-parameter 

.PHONY: all clean
all:	damocles

source/%.o: source/%.f90
	$(FC) $(FFLAGS) $^ -c -o $@

damocles: source/input.o source/electron_scattering.o source/random_routines.o source/vector_functions.o source/grain_sizes.o source/initialise.o source/init_packet.o source/BHmie.o source/init_random_seed.o \
source/propagate.o source/construct_grid.o source/model_comparison.o source/driver.o source/damocles_wrap.o source/damocles.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^


clean:
	rm damocles source/*.o	