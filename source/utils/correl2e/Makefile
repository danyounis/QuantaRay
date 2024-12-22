FC=gfortran# 9.4.0 / 11.2.0
H5=h5fc

PRECISION=DOUBLE

ifndef QRAY_MOD_PATH
QRAY_MOD_PATH=~/QuantaRay/source/modules
endif
MOD_OBJS=prec.o math.o

all: correl2e tidy

correl2e: $(MOD_OBJS) correl2e.o
	$(H5) $^ -o correl2e -fopenmp -llapack -lblas

prec.o: $(QRAY_MOD_PATH)/prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: $(QRAY_MOD_PATH)/math.f08
	$(FC) -c $< -fopenmp

correl2e.o: correl2e.f08
	$(H5) -c $< -fopenmp -llapack -lblas

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod correl2e
	@echo Done.
