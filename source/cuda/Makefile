FC=nvfortran

PRECISION=SINGLE

ifndef QRAY_MOD_PATH
QRAY_MOD_PATH=~/QuantaRay/source/modules
endif
MOD_OBJS=prec.o math.o quantum.o emfm.o rochester.o cufft.o

CUDALIBS=-lcufft -Mcuda -L/opt/nvidia/hpc_sdk/Linux_x86_64/20.7/math_libs/11.0/lib64

all: 2d1e 1d2e tidy

2d1e: $(MOD_OBJS) qcr2d1e.o
	$(FC) $^ -o qray2d1eQ.cuda -fopenmp $(CUDALIBS)

1d2e: $(MOD_OBJS) qcr1d2e.o
	$(FC) $^ -o qray1d2eQ.cuda -fopenmp $(CUDALIBS)

prec.o: $(QRAY_MOD_PATH)/prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: $(QRAY_MOD_PATH)/math.f08
	$(FC) -c $<

quantum.o: $(QRAY_MOD_PATH)/quantum.f08
	$(FC) -c $<

emfm.o: $(QRAY_MOD_PATH)/emfm.f08
	$(FC) -c $<

rochester.o: $(QRAY_MOD_PATH)/rochester.f08
	$(FC) -c $<

cufft.o: $(QRAY_MOD_PATH)/cufft.cuf
	$(FC) -c $<

qcr2d1e.o: qcr2d1e.cuf
	$(FC) -c $< -fopenmp

qcr1d2e.o: qcr1d2e.cuf
	$(FC) -c $< -fopenmp

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod *.cuda
	@echo Done.
