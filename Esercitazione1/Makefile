CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: Es01_1.exe Es01_2.exe Es01_3.exe

Es01_1.exe: Es01_1.o random.o
	$(CC) $^ -o $@

Es01_2.exe: Es01_2.o random.o
	$(CC) $^ -o $@

Es01_3.exe: Es01_3.o random.o
	$(CC) $^ -o $@

%.o: %.cpp %.h
	$(CC) -c $< -o $@ $(CFLAGS)

random.o: random.cpp
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o Es01_1.exe Es01_2.exe Es01_3.exe seed.out
