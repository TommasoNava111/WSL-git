CC = g++
CFLAGS = -Wall -O3 --std=c++11

Es03.1.exe : Es03.1.o random.o
	$(CC) random.o Es03.1.o -o Es03.1.exe
Es03.1.o : Es03.1.cpp
	$(CC) -c Es03.1.cpp -o Es03.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Es03.1.exe seed.out
