PROG = main

CC = g++
CXXFLAGS = -std=c++17 -Wall -flto -Ofast
#DEBUGFLAGS = -g
#FASTFLAGS = -flto -Ofast


MAKEFLAGS = -s

INCPATH = inc/

HDRPATH = $(INCPATH)hdr/
LIBPATH = $(INCPATH)lib/
OBJPATH = $(INCPATH)obj/
SRCPATH = $(INCPATH)src/

EXT = $(notdir $(basename $(wildcard $(HDRPATH)*.hpp)))
HDR = $(addprefix $(HDRPATH),$(addsuffix .hpp, $(EXT)))
OBJ = $(addprefix $(OBJPATH),$(addsuffix .o, $(EXT)))
SRC = $(addprefix $(SRCPATH),$(addsuffix .cpp, $(EXT)))

PROG_OBJ = $(addprefix $(OBJPATH),$(addsuffix .o,$(PROG)))

.PHONY: all
all: $(PROG)

$(PROG): $(OBJ) $(PROG_OBJ)
	$(CC) $(CXXFLAGS) $^ -o $@

$(OBJ): $(SRC)
	mkdir -p $(OBJPATH)
	$(CC) $(CXXFLAGS) -c $(SRCPATH)$(addsuffix .cpp, $(notdir $(basename $@))) -I$(HDRPATH) -o $@

$(PROG_OBJ): $(PROG).cpp
	mkdir -p $(OBJPATH)
	$(CC) $(CXXFLAGS) -c $< -I$(HDRPATH) -o $@
	
.PHONY: clean
clean:
	rm -f $(OBJ) $(PROG_OBJ) $(PROG)

.PHONY: purge
purge:
	rm -f $(OBJPATH)*.o *.o

.PHONY: again
again:
	rm -f $(OBJ) $(PROG_OBJ) $(PROG)
	make

.PHONY: test
test:
	@echo "$(OBJ)"
