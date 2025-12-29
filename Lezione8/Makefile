CC = g++
CFLAGS = -Wall -O3 -std=c++11

OBJS = random.o
TARGETS = main.exe main2.exe main3.exe

all: $(TARGETS)

# ===== main =====
main.exe: main.o $(OBJS)
	$(CC) $^ -o $@

main.o: main.cpp random.h
	$(CC) $(CFLAGS) -c $< -o $@

# ===== main2 =====
main2.exe: main2.o $(OBJS)
	$(CC) $^ -o $@

main2.o: main2.cpp random.h
	$(CC) $(CFLAGS) -c $< -o $@

# ===== main3 =====
main3.exe: main3.o $(OBJS)
	$(CC) $^ -o $@

main3.o: main3.cpp random.h
	$(CC) $(CFLAGS) -c $< -o $@

# ===== random =====
random.o: random.cpp random.h
	$(CC) $(CFLAGS) -c $< -o $@

# ===== clean =====
clean:
	rm -f *.o *.exe seed.out
