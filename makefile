LINK=-lmkl_rt
DEBUG=-g
CC=g++

LFM_radar: util.o LFM_radar.o
	$(CC) -o LFM_radar util.o LFM_radar.o $(LINK) $(DEBUG)
	#icpc -debug -o0 LFM_radar util.o LFM_radar.o $(LINK)

LFM_radar.o: util.h LFM_radar.h LFM_radar.cpp
	$(CC) -c LFM_radar.cpp $(DEBUG)

MKLtest: util.o MKLtest.o
	$(CC) -o MKLtest util.o MKLtest.o $(LINK) $(DEBUG)

util.o: util.h util.cpp
	$(CC) -c util.cpp $(DEBUG)

MKLtest.o: util.h MKLtest.cpp
	$(CC) -c MKLtest.cpp $(DEBUG)

clean:
	rm  *.o
