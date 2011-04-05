# Set your C compiler here...
GCC=gcc
# choose this line for x86 processors...
X86=-DX86
# choose this line for other processors...
#X86=
# end of user settings ###########################################

Ode: Ode.o formulc.o random.o sum.o dop853b.o
	$(GCC) Ode.o formulc.o  random.o sum.o dop853b.o -o Ode -lm

dop853b.o: dop853b.c dop853b.h
	$(GCC) -O3 -c dop853b.c
sum.o : sum.c 
	$(GCC) $(X86) -O3 -c sum.c
Ode-new.o : Ode-new.c 
	$(GCC) -O3 -c Ode-new.c
Ode.o : Ode.c 
	$(GCC) -O3 -c Ode.c
formulc.o : formulc.c formulc.h
	$(GCC) -O3 -c formulc.c
random.o: random.c
	$(GCC) -O3 -c random.c

html:
	update Ode.html-in > Ode.html

DIR=Ode/
tar: 
	(cd ..;  tar czvf Ode.tgz $(DIR)dop853b.c $(DIR)dop853b.h $(DIR)dopri5.h $(DIR)sum.c $(DIR)Ode.c $(DIR)formulc.h $(DIR)formulc.c $(DIR)makefile $(DIR)random.c $(DIR)sum.h $(DIR)ode.py $(DIR)examples $(DIR)Ode.html $(DIR)Ode.jpg $(DIR)updated.gif)

clean:  
	rm -f formulc *.o core
