#   bfqutils: Ben's FastQ Utilities
#   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

.PHONY: test

CC      ?=cc
CFLAGS  +=-std=gnu99
LDFLAGS +=-std=gnu99
LDLIBS  +=-lm
ZDIR    ?=./libs/zlib
ZLIB=
PREFIX  ?=/usr/local
BINDIR  ?=bin

ifeq ($(z_dyn),)
	ZLIB=$(ZDIR)/libz.a
else
	LDLIBS+=-lz
endif

release: CFLAGS+=-O3
release: LDFLAGS+=-O3
release: bfqutils

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer
debug: LDFLAGS+=-g3 -Og -fsanitize=address,undefined
debug: bfqutils

libz/libz.a:
	(cd $(ZDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZDIR)

libz: libz/libz.a

clean/libz:
	$(MAKE) -C $(ZDIR) clean

clean/bfq:
	-rm -f ./src/*.o
	-rm -f ./bfqutils

clean: clean/libz clean/bfq

src/bfqmerge.o: src/bfqmerge.c
	$(CC) $(CFLAGS) -c $^ -o $@

src/bfqstats.o: src/bfqstats.c
	$(CC) $(CFLAGS) -c $^ -o $@

src/bfqtrimpe.o: src/bfqtrimpe.c
	$(CC) $(CFLAGS) -c $^ -o $@

src/bfqtrimse.o: src/bfqtrimse.c
	$(CC) $(CFLAGS) -c $^ -o $@

src/bfqutils.o: src/bfqutils.c
	$(CC) $(CFLAGS) -c $^ -o $@

bfqutils: src/bfqmerge.o src/bfqstats.o src/bfqtrimpe.o src/bfqtrimse.o src/bfqutils.o
	$(CC) $(LDFLAGS) $^ -o $@ $(ZLIB) $(LDLIBS)

test: bfqutils
	(cd ./test && bash test.sh)

install: bfqutils
	install -p ./bfqutils $(PREFIX)/$(BINDIR)

uninstall:
	-rm -f $(PREFIX)/$(BINDIR)/bfqutils

