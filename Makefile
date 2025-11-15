# Makefile for TRTS Engine Reconstruction
#
# Builds the complete TRTS engine with strict rational arithmetic
# and deterministic propagation logic.

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -pedantic
LDFLAGS = -lgmp -lm

# Core library sources
CORE_SRCS = rational.c state.c config.c psi.c koppa.c engine.c simulate.c \
            config_loader.c analysis_utils.c

CORE_OBJS = $(CORE_SRCS:.c=.o)

# Main programs
PROGRAMS = trts_simulate trts_go_time

.PHONY: all clean

all: $(PROGRAMS)

# Core library
libtrts.a: $(CORE_OBJS)
	ar rcs $@ $^

# Simulation program (writes CSV files)
trts_simulate: libtrts.a
	$(CC) $(CFLAGS) -o $@ -c trts_simulate_main.c
	$(CC) $(CFLAGS) -o $@ trts_simulate_main.o libtrts.a $(LDFLAGS)

# Minimal CLI program
trts_go_time: libtrts.a
	$(CC) $(CFLAGS) -o $@ -c trts_go_time_main.c
	$(CC) $(CFLAGS) -o $@ trts_go_time_main.o libtrts.a $(LDFLAGS)

# Object file rules
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Dependencies (simplified - in production use makedepend or similar)
rational.o: rational.c rational.h
state.o: state.c state.h rational.h config.h
config.o: config.c config.h rational.h
psi.o: psi.c psi.h config.h state.h rational.h
koppa.o: koppa.c koppa.h config.h state.h rational.h
engine.o: engine.c engine.h config.h state.h rational.h
simulate.o: simulate.c simulate.h config.h state.h engine.h koppa.h psi.h rational.h
config_loader.o: config_loader.c config_loader.h config.h rational.h
analysis_utils.o: analysis_utils.c analysis_utils.h config.h state.h simulate.h rational.h

clean:
	rm -f $(CORE_OBJS) $(PROGRAMS) *.o libtrts.a
	rm -f events.csv values.csv

# Test targets
test: trts_simulate
	./trts_simulate --ticks 10

# Example: Run with golden ratio seeds
example_golden: trts_go_time
	./trts_go_time --ticks 50 --ups 3/2 --beta 5/3 --koppa 1/1

.PHONY: test example_golden
