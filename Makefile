# The software kappa is a collection of programs to compute the homology of
# the moduli space of surfaces using the radial model.
# Copyright (C) 2015 - 2018  Felix Boes and Anna Hermann
# 
# This file is part of kappa.
# 
# kappa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# kappa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with kappa.  If not, see <http:#www.gnu.org/licenses/>.


LANG = en_US.UTF-8

CXXFLAGS      := -O3 -std=c++11 -D_GLIBCXX_USE_NANOSLEEP \
                 -Wextra -Wall -Wno-long-long -pedantic-errors
INCL          := -I.
BUILDDIR      := build
EXT           := cpp
SRCDIRS       := kappa libhomology
EXCLUDE       := test_ version
.DEFAULT_GOAL := compute_css

EXECUTE_BEFORE_BUILD_CHEAT := $(shell echo "const char* program_version_by_git = \"$(shell git rev-parse HEAD 2>/dev/null)\";" > kappa/version.cpp)

ifdef ENABLE_TCMALLOC
CXXFLAGS      += -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
LIBS          := -ltcmalloc -lpthread
else
LIBS          := -lpthread
endif

ifdef ENABLE_CUSTOM_LIBS
LIBS          += -L./lib/
INCL          += -I./include/
endif

LIBS          += -lboost_filesystem -lboost_system -lboost_iostreams \
                 -lboost_serialization -lboost_program_options \
                 -lboost_date_time -lgmpxx -lgmp

ifdef ENABLE_KAPPA_DEBUG
CXXFLAGS      += -g
endif

ifdef USE_LIBMAGICKXX
CXXFLAGS      += `Magick++-config --cppflags --cxxflags --ldflags --libs`
endif

ifndef CXX
CXX           := g++
endif

ifndef DOXYGEN
DOXYGEN       := doxygen
endif

ifeq ($(shell expr `$(CXX) -dumpversion` \<= 4.7.7), 1)
CXXFLAGS      := $(patsubst -std=c++11,-std=c++0x, $(CXXFLAGS))
endif
ifeq ($(shell expr `$(CXX) -dumpversion` \>= 4.9.0), 1)
CXXFLAGS      := $(CXXFLAGS) -fdiagnostics-color=auto  -fsanitize=undefined
endif
ifeq ($(shell expr `$(CXX) -dumpversion` \< 4.7), 1)
CXXFLAGS      := $(CXXFLAGS) -DBROKEN_VECTOR_IMPLEMENTATION -DBROKEN_UNIQUE_PTR_IMPLEMENTATION
endif
ifeq ($(shell expr `doxygen --version` \>= 1.8.7),1)
DOXYGENFLAGS  := $(DOXYGENFLAGS) -d Validate
endif

override BUILDDIR := $(strip $(BUILDDIR))
CXXSRC    := $(filter-out $(foreach d,$(SRCDIRS), $(foreach e,$(EXCLUDE), $(d)/$(e)%.$(EXT))), $(foreach d,$(SRCDIRS), $(wildcard $(d)/*.$(EXT))))
CXXOBJ    := $(patsubst %,build/%, $(CXXSRC:.$(EXT)=.o))
CXXDEP    := $(CXXOBJ:.o=.dep)
TAGS      := $(patsubst %,$(BUILDDIR)/%/.tag, $(SRCDIRS))

OBJ    := $(CXXOBJ)
STDOBJ := $(filter-out $(foreach d,$(SRCDIRS), build/$(d)/main_%.o), $(OBJ))

TARGETS := $(patsubst kappa/main_%.$(EXT),%, $(filter $(foreach d,$(SRCDIRS),$(d)/main_%.$(EXT)), $(foreach d,$(SRCDIRS), $(wildcard $(d)/*.$(EXT)))))

ifneq ($(MAKECMDGOALS),clean)
-include $(CXXDEP) $(CDEP)
endif

$(TARGETS): %: $(BUILDDIR)/kappa/main_%.o $(STDOBJ)
	$(CXX) $(CXXFLAGS) $(INCL) -o $@ kappa/version.cpp $(STDOBJ) $(BUILDDIR)/kappa/main_$@.o $(LIBS)

$(CXXOBJ): $(BUILDDIR)/%.o: %.$(EXT) $(BUILDDIR)/%.dep
	$(CXX) $(CXXFLAGS) $(INCL) -c $< -o $@

$(CXXDEP): $(BUILDDIR)/%.dep: %.$(EXT) $(TAGS)
	$(CXX) $(CXXFLAGS) $(INCL) -MM $< -MT $@ -MT $(<:.$(EXT)=.o) -o $@

overflowtest:
	$(CXX) -ftrapv overflowtest.cpp -o overflowtest_cxx
	clang++ -ftrapv overflowtest.cpp -o overflowtest_clang

%.tag:
	@mkdir -p $(dir $(@))
	@touch $@

.PHONY: doc
doc:
	$(DOXYGEN) ./Doxyfile $(DOXYGENFLAGS)
	#make --directory=./doc/latex/ pdf
	#@find ./doc/latex/ -regex ".*/refman\(_2on1\)?\.\(dvi\|ps\|pdf\)" -exec mv {} doc/ \;

.PHONY: clean
clean:
	rm -rf $(TARGETS) $(BUILDDIR)
	rm -rf overflowtest_cxx overflowtest_clang
	rm -rf version
	rm -rf doc/html/ doc/latex/ doc/*.pdf

.PHONY: clean_cache
clean_cache:
	rm -rf cache
