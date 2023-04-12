# Copyright (C) 2023 Simon Crase. simon@greenweavez.nz
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
# Makefile snarfed from https://stackoverflow.com/questions/2481269/how-to-make-a-simple-c-makefile


CPPFLAGS    = -g -O3  -D_RUNTIME_CHECKS -pthread -I/sw/include/root  -std=gnu++11  -DVERSION="\"$(GIT_VERSION)\""
LDFLAGS     = -g -O3
LDLIBS      =
CC          = gcc
CXX         = g++
RM          = rm -f
MKDIR       = mkdir
SRCS        = newick.cpp


TESTS       = 	test-newick.cpp		
		
OBJS1      = $(subst .cpp,.o,$(SRCS)) 
OBJS       = $(subst .cc,.o,$(OBJS1)) 
TEST_OBJS  = $(subst .cpp,.o,$(TESTS))



MAIN      = qrtd.exe
TEST_MAIN = tests.exe 

TARGETS   = $(MAIN) 

all : $(TARGETS) $(TEST_OBJS)

run : all
	${RM}  *.stackdump
	$(MAIN)
	
clean :
	${RM} *.o *.stackdump

rebuild: clean all

depend: .depend

install: rebuild
	cp $(MAIN) /usr/local/bin
	
.depend: $(SRCS)  qrtd.cpp $(TESTS) 
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;
	sed -i -e 's/\/home\/Weka\/qrtd\///$g' .depend
	
$(MAIN): $(OBJS) qrtd.o
	${CXX} $(LDFLAGS) -o $(MAIN) qrtd.o ${OBJS} ${LDLIBS}

$(TEST_MAIN): $(OBJS) tests.o $(TEST_OBJS)
	${CXX} $(LDFLAGS) -o $(TEST_MAIN) tests.o ${OBJS} $(TEST_OBJS) ${LDLIBS}

tests : $(TEST_MAIN)
	${RM}  *.stackdump
	./$(TEST_MAIN)
	
distclean: clean
	$(RM) *~ .depend

setup:
	-$(MKDIR) configs
	-$(MKDIR) imgs
	-$(MKDIR) logs
	
include .depend