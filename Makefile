.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean))) #############

CFLAGS := -Wall -O3 -flto -fmax-errors=3 -Iinclude
# CFLAGS := -Wall -O0 -g -fmax-errors=3 -Iinclude
CXXFLAGS := -std=c++20 $(CFLAGS)

# generate .d files during compilation
DEPFLAGS = -MT $@ -MMD -MP -MF .build/$*.d

FIND_MAIN := \
  find src -type f -name '*.cc' \
  | xargs grep -l '^\s*int\s\+main\s*(' \
  | sed 's:^src/\(.*\)\.cc$$:bin/\1:'
EXE := $(shell $(FIND_MAIN))

all: $(EXE)

bin/hgamweb: $(patsubst %, .build/%.o, linalg wls json)
# LF_hgamweb := -static-libstdc++ -static-libgcc
# L_hgamweb :=

#####################################################################

.PRECIOUS: .build/%.o lib/lib%.so

bin/%: .build/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(LF_$*) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

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

