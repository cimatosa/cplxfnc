PREFIX  =  /usr
CC      = /usr/bin/g++
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lm -larb -lflint


cplxfnc_check: cplxfnc_check.cpp cplxfnc.o
	$(CC) -o cplxfnc_check $(CFLAGS) cplxfnc_check.cpp cplxfnc.o $(LDFLAGS)


cplxfnc.o: cplxfnc.cpp
	$(CC) -c -o cplxfnc.o $(LDFLAGS) $(CFLAGS) -fPIC cplxfnc.cpp
	$(CC) -shared -o libcplxfnc.so cplxfnc.o

.PHONY: check
check:
	./cplxfnc_check

.PHONY: clean
clean:
	rm -v -rf *.o *.so

.PHONY: install
install:
	install -m 755 libcplxfnc.so $(PREFIX)/lib
	install -m 644 cplxfnc.hpp $(PREFIX)/include

