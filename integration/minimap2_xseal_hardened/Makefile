# minimap2_xseal: minimap2 with XSeal fused frontend (index.cpp replaces index.c).
CXX ?= clang++
CC ?= clang
CFLAGS ?= -O3 -g -Wall -Wc++-compat -DNDEBUG -DHAVE_KALLOC -march=native
CXXFLAGS ?= -O3 -g -std=c++20 -Wall -DNDEBUG -DHAVE_KALLOC -march=native -Iinclude -I.
LIBS ?= -lm -lz -lpthread

MM_C_SRCS := \
	align.c bseq.c esterr.c format.c hit.c jump.c kalloc.c kthread.c lchain.c \
	map.c misc.c options.c pe.c seed.c sketch.c splitidx.c sdust.c

MM_OBJS := $(MM_C_SRCS:.c=.o) \
	ksw2_ll_sse.o \
	ksw2_extz2_sse41.o ksw2_extz2_sse2.o \
	ksw2_extd2_sse41.o ksw2_extd2_sse2.o \
	ksw2_exts2_sse41.o ksw2_exts2_sse2.o \
	ksw2_dispatch.o \
	index.o main.o

PROG := minimap2_xseal

.PHONY: all clean

all: $(PROG)

$(PROG): $(MM_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(MM_OBJS) $(LIBS)

main.o: main.c
	$(CC) -c $(CFLAGS) $< -o $@

index.o: index.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

ksw2_ll_sse.o: ksw2_ll_sse.c
	$(CC) -c $(CFLAGS) -msse2 $< -o $@

ksw2_extz2_sse41.o: ksw2_extz2_sse.c
	$(CC) -c $(CFLAGS) -msse4.1 -DKSW_CPU_DISPATCH $< -o $@

ksw2_extz2_sse2.o: ksw2_extz2_sse.c
	$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $< -o $@

ksw2_extd2_sse41.o: ksw2_extd2_sse.c
	$(CC) -c $(CFLAGS) -msse4.1 -DKSW_CPU_DISPATCH $< -o $@

ksw2_extd2_sse2.o: ksw2_extd2_sse.c
	$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $< -o $@

ksw2_exts2_sse41.o: ksw2_exts2_sse.c
	$(CC) -c $(CFLAGS) -msse4.1 -DKSW_CPU_DISPATCH $< -o $@

ksw2_exts2_sse2.o: ksw2_exts2_sse.c
	$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $< -o $@

ksw2_dispatch.o: ksw2_dispatch.c
	$(CC) -c $(CFLAGS) -msse4.1 -DKSW_CPU_DISPATCH $< -o $@

clean:
	rm -f $(MM_OBJS) $(PROG)
