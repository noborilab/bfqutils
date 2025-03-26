CC      ?=cc
CFLAGS  +=-std=gnu99
LDLIBS  +=-lm -lpthread -lz
ZDIR    ?=./libs/zlib
PREFIX  ?=/usr/local
BINDIR  ?=bin
MANDIR  ?=share/man/man1
TESTDIR  =./test

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG
debug: fqmerge

release: CFLAGS+=-DNDEBUG -O3
release: fqmerge

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@

objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

fqmerge: $(objects)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objects) -o $@ $(LDLIBS) 

clean:
	-rm -f src/*.o
