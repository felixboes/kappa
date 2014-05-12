LANG = en_US.UTF-8

CXXFLAGS      := -O3 -std=c++11 -D_GLIBCXX_USE_NANOSLEEP -g \
                 -Wextra -Wall -Wno-long-long -pedantic-errors
LIBS          := -lgmpxx -lgmp -lpthread \
                 -lboost_filesystem -lboost_system -lboost_iostreams \
                 -lboost_serialization -lboost_program_options
BUILDDIR      := build
EXT           := cpp
SRCDIRS       := kappa libhomology
EXCLUDE       := test_
.DEFAULT_GOAL := compute_homology

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
ifeq ($(shell expr `$(CXX) -dumpversion` \<= 4.7.0), 1)
CXXFLAGS      := $(CXXFLAGS) -DBROKEN_VECTOR_IMPLEMENTATION -DBROKEN_UNIQUE_PTR_IMPLEMENTATION
endif
ifeq ($(shell expr `doxygen --version` \>= 1.8.7),1)
DOXYGENFLAGS  := $(DOXYGENFLAGS) -d Validate
endif

INCL      := -I.
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

ifeq ($(MAKECMDGOALS),draw_differentials)
LIBS     := $(LIBS) Magick++-config --cppflags --cxxflags --ldflags --libs
CXXFLAGS := $(CXXFLAGS) -DCOMPILE_WITH_MAGICK
endif


$(TARGETS): %: $(BUILDDIR)/kappa/main_%.o $(STDOBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $(INCL) -o $@ $(STDOBJ) $(BUILDDIR)/kappa/main_$@.o

$(CXXOBJ): $(BUILDDIR)/%.o: %.$(EXT) $(BUILDDIR)/%.dep
	$(CXX) $(LIBS) $(CXXFLAGS) $(INCL) -c $< -o $@

$(CXXDEP): $(BUILDDIR)/%.dep: %.$(EXT) $(TAGS)
	$(CXX) $(CXXFLAGS) $(INCL) -MM $< -MT $@ -MT $(<:.$(EXT)=.o) -o $@

%.tag:
	@mkdir -p $(dir $(@))
	@touch $@

.PHONY: doc
doc:
	$(DOXYGEN) doxygen/Doxyfile $(DOXYGENFLAGS)
	make --directory=doxygen/latex/ pdf
	@find doxygen/latex/ -regex ".*/refman\(_2on1\)?\.\(dvi\|ps\|pdf\)" -exec mv {} doxygen/ \;

.PHONY: clean
clean:
	rm -rf $(TARGETS) $(BUILDDIR)
	rm -rf doxygen/html/ doxygen/latex/ doxygen/*.pdf
