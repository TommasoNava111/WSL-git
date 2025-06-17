CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: Es02.1.exe Es02.2.exe

Es02.1.exe: Es02.1.o random.o
	$(CC) $^ -o $@

Es02.2.exe: Es02.2.o random.o
	$(CC) $^ -o $@

%.o: %.cpp %.h
	$(CC) -c $< -o $@ $(CFLAGS)

random.o: random.cpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o Es02.1.exe Es02.2.exe seed.out
