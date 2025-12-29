CC = g++
CFLAGS = -Wall -O3 -std=c++11

# Eseguibili
TARGETS = main.exe main2.exe

all: $(TARGETS)

main.exe: main.o random.o
	$(CC) $(CFLAGS) $^ -o $@

main2.exe: main2.o random.o
	$(CC) $(CFLAGS) $^ -o $@

# Regola generica: .cpp -> .o
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS) seed.out


