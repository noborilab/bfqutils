CC      ?=cc
LDLIBS  +=-lm
ZDIR    ?=./libs/zlib

ZOPS=
ZLIB=

ifeq ($(z_dyn),)
	ZLIB=$(ZDIR)/libz.a
else
	LDLIBS+=-lz
endif

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG
debug: fqmerge

release: CFLAGS+=-DNDEBUG -O3
release: fqmerge

libz/libz.a:
	(cd $(ZDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZDIR)

libz: libz/libz.a

clean/libz:
	$(MAKE) -C $(ZDIR) clean

clean/fqmerge:
	-rm -f src/*.o
	-rm -f ./fqmerge

clean: clean/libz clean/fqmerge

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@

objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

fqmerge: $(objects)
	$(CC) $(CFLAGS) $(LDFLAGS) $(objects) -o $@ $(ZLIB) $(LDLIBS) 

