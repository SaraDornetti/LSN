CC = mpic++
CFLAGS = -O3 -std=c++17
AFLAGS = -larmadillo

all: main.exe

main.exe: main.o system.o individuo.o random.o
	$(CC) system.o individuo.o random.o main.o -o main.exe $(AFLAGS)

main.o: main.cpp system.h individuo.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

system.o: system.cpp system.h individuo.h random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

individuo.o: individuo.cpp individuo.h
	$(CC) -c individuo.cpp -o individuo.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o main.exe

