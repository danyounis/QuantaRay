FC=gfortran# 9.4.0 / 11.2.0
H5=h5fc

PRECISION=DOUBLE

ifndef QRAY_MOD_PATH
QRAY_MOD_PATH=~/QuantaRay/source/modules
endif
MOD_OBJS=prec.o math.o quantum.o emfm.o vdm.o rochester.o

all: qcr ccr tidy

qcr: $(MOD_OBJS) dataio.o main_qcr.o
	$(H5) $^ -o qrayQ -fopenmp

ccr: $(MOD_OBJS) dataio.o main_ccr.o
	$(H5) $^ -o qrayC -fopenmp

prec.o: $(QRAY_MOD_PATH)/prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: $(QRAY_MOD_PATH)/math.f08
	$(FC) -c $< -fopenmp

quantum.o: $(QRAY_MOD_PATH)/quantum.f08
	$(FC) -c $< -fopenmp

emfm.o: $(QRAY_MOD_PATH)/emfm.f08
	$(FC) -c $<

vdm.o: $(QRAY_MOD_PATH)/vdm.f08
	$(FC) -c $<

rochester.o: $(QRAY_MOD_PATH)/rochester.f08
	$(FC) -c $<

dataio.o: dataio.f90
	$(H5) -c $<

main_qcr.o: main_qcr.f08
	$(FC) -c $< -fopenmp

main_ccr.o: main_ccr.f08
	$(FC) -c $< -fopenmp

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod qrayQ qrayC
	@echo Done.
