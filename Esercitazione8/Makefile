CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: Es08.1.exe Es08.2.exe

Es08.1.exe: Es08.1.o random.o
	$(CC) $^ -o $@

Es08.2.exe: Es08.2.o random.o
	$(CC) $^ -o $@

%.o: %.cpp %.h
	$(CC) -c $< -o $@ $(CFLAGS)

random.o: random.cpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o Es08.1.exe Es08.2.exe seed.out
