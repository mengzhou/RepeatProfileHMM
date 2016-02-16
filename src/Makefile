ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(PROGRAM_ROOT)/smithlab_cpp
ifeq ("$(wildcard $(SMITHLAB_CPP))","")
$(error SMITHLAB_CPP not set and smithlab_cpp not found)
endif
all_subdirs += $(SMITHLAB_CPP)
lib_subdirs += $(SMITHLAB_CPP)
endif

PROGS = test

CXX = /home/rcf-40/mengzhou/panfs/tools/gcc-5.3.0/bin/g++
CXXFLAGS = -Wall -fmessage-length=50 -std=c++11 -static-libstdc++
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

INCLUDEDIRS = $(SMITHLAB_CPP)

INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lgsl -lgslcblas

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
INCLUDEARGS += -I/usr/local/include
LIBS += -L/usr/local/lib
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

all: $(PROGS)

install: $(PROGS)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, smithlab_utils.o smithlab_os.o OptionParser.o)\
  ProfileHMM.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
