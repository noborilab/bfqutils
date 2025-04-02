#   fqmerge: FastQ merging tool for PE reads
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
LDLIBFS +=-lm
ZDIR    ?=./libs/zlib
ZLIB=

ifeq ($(z_dyn),)
	ZLIB=$(ZDIR)/libz.a
else
	LDLIBS+=-lz
endif

release: CFLAGS+=-O3
release: bfqmerge bfqstats

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer
debug: bfqmerge bfqstats

libz/libz.a:
	(cd $(ZDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZDIR)

libz: libz/libz.a

clean/libz:
	$(MAKE) -C $(ZDIR) clean

clean/bfq:
	-rm -f ./src/bfqmerge.o
	-rm -f ./src/bfqstats.o
	-rm -f ./bfqmerge
	-rm -f ./bfqstats

clean: clean/libz clean/bfq

src/bfqmerge.o: src/bfqmerge.c
	$(CC) $(CFLAGS) -c $^ -o $@

src/bfqstats.o: src/bfqstats.c
	$(CC) $(CFLAGS) -c $^ -o $@

bfqmerge: src/bfqmerge.o
	$(CC) $(CFLAGS) $^ -o $@ $(ZLIB)

bfqstats: src/bfqstats.o
	$(CC) $(CFLAGS) $^ -o $@ $(ZLIB)

test: bfqmerge bfqstats
	(cd ./test && bash test.sh)

