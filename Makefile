LANG = en_US.UTF-8

CXXFLAGS      := -O3 -std=c++11 -D_GLIBCXX_USE_NANOSLEEP \
                 -Wextra -Wall -Wno-long-long -pedantic-errors
LIBS          := -lpthread \
                 -lboost_filesystem -lboost_system -lboost_iostreams \
                 -lboost_serialization -lboost_program_options \
                 -lboost_date_time
BUILDDIR      := build
EXT           := cpp
SRCDIRS       := kappa libhomology
EXCLUDE       := test_ version
.DEFAULT_GOAL := compute_css

EXECUTE_BEFORE_BUILD_CHEAT := $(shell echo "const char* program_version_by_git = \"$(shell git rev-parse HEAD 2>/dev/null)\";" > kappa/version.cpp)

ifdef ADV_OPTIMIZATION
LIBS          += -ltcmalloc lib/libgmpxx.a lib/libgmp.a
INCL          := -I. -I./include/
else
LIBS          += -lgmpxx -lgmp
INCL          := -I.
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
	rm -rf version
	rm -rf doc/html/ doc/latex/ doc/*.pdf

.PHONY: clean_cache
clean_cache:
	rm -rf cache
