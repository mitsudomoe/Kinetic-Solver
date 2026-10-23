CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Wshadow
LDLIBS  = -lm
PYTHON  ?= python3

BUILD   = build
SOURCES = src/grid.c src/quadrature.c src/kinetic.c src/rusanov.c src/scheme.c src/integrator.c src/simulation.c src/cases.c src/output.c
OBJECTS = $(SOURCES:src/%.c=$(BUILD)/%.o)
PROGRAM = $(BUILD)/kinetic-solver
HARNESS = $(BUILD)/flux-harness

.PHONY: all test docs doc-figures clean

all: $(PROGRAM)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c src/*.h | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGRAM): $(OBJECTS) $(BUILD)/main.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(HARNESS): tests/flux_harness.c $(OBJECTS)
	$(CC) $(CFLAGS) -Isrc $^ -o $@ $(LDLIBS)

test: $(PROGRAM) $(HARNESS)
	$(PYTHON) -m unittest discover -s tests -v

docs:
	mkdocs build --strict

doc-figures: $(PROGRAM)
	$(PYTHON) scripts/make_doc_figures.py

clean:
	rm -rf $(BUILD) output site
