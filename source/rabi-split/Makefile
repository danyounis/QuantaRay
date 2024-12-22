FC=gfortran# 9.4.0 / 11.2.0
H5=h5fc

PRECISION=DOUBLE

ifndef QRAY_MOD_PATH
QRAY_MOD_PATH=~/QuantaRay/source/modules
endif
MOD_OBJS=prec.o math.o quantum.o emfm.o

all: main tidy

main: $(MOD_OBJS) dataio.o main.o
	$(H5) $^ -o rabi_split -fopenmp

prec.o: $(QRAY_MOD_PATH)/prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: $(QRAY_MOD_PATH)/math.f08
	$(FC) -c $< -fopenmp

quantum.o: $(QRAY_MOD_PATH)/quantum.f08
	$(FC) -c $< -fopenmp

emfm.o: $(QRAY_MOD_PATH)/emfm.f08
	$(FC) -c $<

dataio.o: dataio.f90
	$(H5) -c $<

main.o: main.f08
	$(FC) -c $< -fopenmp

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod rabi_split
	@echo Done.
