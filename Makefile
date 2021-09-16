.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean))) #############

CFLAGS := -Wall -O3 -flto -fmax-errors=3 -Iinclude
# CFLAGS := -Wall -O0 -g -fmax-errors=3 -Iinclude
CXXFLAGS := -std=c++20 $(CFLAGS)

# generate .d files during compilation
DEPFLAGS = -MT $@ -MMD -MP -MF .build/$*.d

ROOT_CPPFLAGS := $(shell root-config --cflags \
  | sed 's/-std=[^ ]*//g;s/-I/-isystem /g')
ROOT_LIBDIR   := $(shell root-config --libdir)
ROOT_LDFLAGS  := $(shell root-config --ldflags) -Wl,-rpath,$(ROOT_LIBDIR)
ROOT_LDLIBS   := $(shell root-config --libs)

FIND_MAIN := \
  find src -type f -name '*.cc' \
  | xargs grep -l '^\s*int\s\+main\s*(' \
  | sed 's:^src/\(.*\)\.cc$$:bin/\1:'
EXE := $(shell $(FIND_MAIN))

all: $(EXE)

bin/hgamweb: $(patsubst %, .build/%.o, linalg wls json)
# C_hgamweb := -ffast-math
# LF_hgamweb := -static -static-libstdc++ -static-libgcc
L_hgamweb := -lgsl -lgslcblas

C_make_vars := $(ROOT_CPPFLAGS)
LF_make_vars := $(ROOT_LDFLAGS)
L_make_vars := -L$(ROOT_LIBDIR) -lCore -lRIO -lTree -lTreePlayer -lHist

# STRIP := hgamweb

#####################################################################

.PRECIOUS: .build/%.o lib/lib%.so

bin/%: .build/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(LF_$*) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)
	$(if $(filter $*,$(STRIP)), strip -s $@)

lib/lib%.so:
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(LF_$*) -shared $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

.build/%.o: src/%.cc
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(C_$*) -c $(filter %.cc,$^) -o $@

.build/%.o: src/%.c
	@mkdir -pv $(dir $@)
	$(CC) $(CFLAGS) $(DEPFLAGS) $(C_$*) -c $(filter %.c,$^) -o $@

-include $(shell [ -d '.build' ] && find .build -type f -name '*.d')

endif ###############################################################

clean:
	@rm -rfv bin lib .build

