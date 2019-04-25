CC= gcc
CFLAGS= -g -Wall -O2
INC=
LIB= -Lhtslib -lhts -lpthread -lz -lm
OBJS= bam.o deriv.o main.o map.o ped.o util.o mats.o

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(INC)

bsh-ref: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(INC) $(LIB)

clean:
	rm -f $(OBJS)

submodules:
	cd htslib && make libhts.a ; cd ..

# DO NOT DELETE THIS LINE -- make depend depends on it.

bam.o main.o: bam.h
bam.o deriv.o main.o ped.o: deriv.h
bam.o main.o map.o ped.o: map.h
bam.o deriv.o main.o ped.o: ped.h
bam.o main.o ped.o util.o: util.h
mats.o ped.o: mats.h
