# Copyright (C) 2025 Simon Crase
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>
#
# Makefile snarfed from
# https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile


CPP_BASIC_FLAGS = -g -O3  -I/sw/include/root -Wall -std=c++23
CPPFLAGS  =  $(CPP_BASIC_FLAGS) 
LDFLAGS   = -g -O3
LDLIBS    = 
CC        = gcc
CXX       = g++
RM        = rm -f
MKDIR     = mkdir -p
SRCS      =  

TESTS     = 

OBJDIR = obj
		
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))
TEST_OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(TESTS))


MAIN      = rosalind.exe
TARGETS   = $(MAIN) 
TEST_MAIN = tests.exe

all : $(MAIN) $(TEST_OBJS)

run : all 
	${RM} configs/foo*.csv
	./$(MAIN)
	
tests : $(TEST_MAIN) all
	./$(TEST_MAIN)
	
clean :
	${RM} $(OBJDIR)/*.o *.stackdump

rebuild: clean all

depend: Makefile .depend

install: rebuild
	cp $(MAIN) /usr/local/bin
	mkdir -p ../logs
	mkdir -p ../config
	
.depend: $(SRCS) $(TESTS)  rosalind.cpp Makefile
	$(RM) ./.depend
	$(CXX) $(CPP_BASIC_FLAGS) -MM $^>>./.depend;
	sed -i -e 's/^.*:.*/$(OBJDIR)\/&/' .depend;

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) -c $< -o $@
		
$(MAIN): $(OBJS) $(OBJDIR)/rosalind.o 
	${CXX} $(LDFLAGS) -o $(MAIN) $(OBJDIR)/rosalind.o ${OBJS} ${LDLIBS}
	
$(TEST_MAIN): $(OBJS) $(OBJDIR)/tests.o $(TEST_OBJS)
	${CXX} $(LDFLAGS) -o $(TEST_MAIN) $(OBJDIR)/tests.o ${OBJS} $(TEST_OBJS) ${LDLIBS}
	
distclean: clean
	$(RM) *~ .depend

setup:
	-$(MKDIR) ../logs
	-$(MKDIR) ../config

include .depend
