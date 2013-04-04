FFLAGS = -O3 -fPIC
FC     = gfortran

OBJS  =  metan1.o linalg_metan1.o zibconst.o

libmetan1.a: $(OBJS)
	$(FC) -shared $(OBJS) -o libmetan1.so

clean:
	rm -f $(OBJS) libmetan1.so
