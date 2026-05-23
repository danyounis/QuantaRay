# FC=gfortran # 9.4.0 / 11.2.0
FC=mpifort # 13.3.0

PRECISION=DOUBLE

# Do not print warnings
NOWARN=TRUE
# Do not perform some logical error checks
UNSAFE=FALSE

ifndef QRAY_MOD_PATH
QRAY_MOD_PATH=~/QuantaRay/source/modules
endif
MOD_OBJS=prec.o math.o quantum.o optimize.o qed_mcpm.o

CFLAGS = -D$(PRECISION)
ifeq ($(NOWARN),TRUE)
	CFLAGS += -DNOWARN
endif
ifeq ($(UNSAFE),TRUE)
	CFLAGS += -DUNSAFE
endif

make: clean wipe all

all: main tidy

main: $(MOD_OBJS) main.o
	$(FC) $^ -o prstrm -fopenmp

prec.o: $(QRAY_MOD_PATH)/prec.F08
	$(FC) -c -D$(PRECISION) $<

math.o: $(QRAY_MOD_PATH)/math.f08
	$(FC) -c -D$(PRECISION) $<

quantum.o: $(QRAY_MOD_PATH)/quantum.f08
	$(FC) -c $< -fopenmp

optimize.o: $(QRAY_MOD_PATH)/optimize.f08
	$(FC) -c -D$(PRECISION) $<

qed_mcpm.o: $(QRAY_MOD_PATH)/qed_mcpm.F08
	$(FC) -c $(CFLAGS) $<

main.o: main.f08
	$(FC) -c $< -fopenmp

tidy:
	rm -f *.o *.mod
	@echo Build successful.

clean:
	rm -f *.o *.mod prstrm
	@echo Done.
