
TARGET := 




############### DEFINE FLAGS ####################
DEF_FLAGS = 

#DEF_FLAGS += -D_MY_VERBOSE_LESS
DEF_FLAGS += -D_MY_VERBOSE
#DEF_FLAGS += -D_NLOGSCALE
#DEF_FLAGS += -D_MY_VERBOSE_MORE
#DEF_FLAGS += -D_MY_VERBOSE_TEDIOUS
#DEF_FLAGS += -D_OLD_IO


############### GIT VERSION? ####################
# Adds git version to DEF_FLAGS
IS_GIT_REPOSITORY := $(shell ((git status >/dev/null 2>&1) && echo yes) || (echo no))
ifeq ($(IS_GIT_REPOSITORY), yes)
    VERSION_STRING:=$(shell git rev-parse --short HEAD | tr -d '\n'; (git status --short | grep -v -q -e '^?') && echo '-mod')
    DEF_FLAGS += -DVERSION_STRING='"$(VERSION_STRING)"'
endif 

################# COMPILER ######################
CXX = g++

CXX_FLAGS =
CXX_FLAGS += -std=gnu++11
# CXX_FLAGS += -O2 
CXX_FLAGS += -g
# CXX_FLAGS += -Wall 
# CXX_FLAGS += -fmessage-length=0


############## BUILD DIRECTORY ##################
BUILD_DIR = build


################### LINKER ######################
LIB  = -lm
LIB += -fopenmp
# LIB += -lgsl -lgslcblas -L/usr/include/gsl

LIB_LINUX  = -lfftw3 -lfftw3_threads
LIB_LINUX += -lhdf5 -lhdf5_hl -lechelon -lechelon_hdf5 

LIB_WIN32  = ./lib/libfftw3-3.dll


#################### INCLUDES ###################
INC  = -I./include 
# INC += -I/usr/include/gsl
INC_LINUX =
INC_WIN32 = 


################### SOURCES #####################
SRC_DIR = src
SRC_EXT = cpp

SRC_COMMON   := $(shell find $(SRC_DIR)/common -type f -name \*.$(SRC_EXT))
SRC_KIELFLOW := $(shell find $(SRC_DIR)/kielflow -type f -name \*.$(SRC_EXT))
SRC_FRONTEND := $(shell find $(SRC_DIR)/frontend -type f -name \*.$(SRC_EXT))
SRC_BACKEND  := $(shell find $(SRC_DIR)/backend -type f -name \*.$(SRC_EXT))

SRC_ALL = $(SRC_COMMON) $(SRC_KIELFLOW) $(SRC_FRONTEND) $(SRC_BACKEND)

vpath = $(SRC_COMMON) $(SRC_KIELFLOW) $(SRC_FRONTEND) $(SRC_BACKEND) # %.cpp $(SRC_DIR)
#vpath %.cc $(SRC_ALL)
#vpath %.c $(SRC_ALL)

################### OBJECTS #####################
OBJ_COMMON   := $(patsubst $(SRC_DIR)/common/%,$(BUILD_DIR)/common/%,$(SRC_COMMON:.$(SRC_EXT)=.o))
OBJ_KIELFLOW := $(patsubst $(SRC_DIR)/kielflow/%,$(BUILD_DIR)/kielflow/%,$(SRC_KIELFLOW:.$(SRC_EXT)=.o))
OBJ_FRONTEND := $(patsubst $(SRC_DIR)/frontend/%,$(BUILD_DIR)/frontend/%,$(SRC_FRONTEND:.$(SRC_EXT)=.o))
OBJ_BACKEND  := $(patsubst $(SRC_DIR)/backend/%,$(BUILD_DIR)/backend/%,$(SRC_BACKEND:.$(SRC_EXT)=.o))


################## OS DETECT#####################
ifeq ($(OS),Windows_NT)
	LIB += $(LIB_WIN32)
	INC += $(INC_WIN32)
	EXC = run_WIN32
	RMV = del
endif

ifeq ($(shell uname), Linux)
	LIB += $(LIB_LINUX)
	INC += $(INC_LINUX)
	EXC = run_LINUX
	RMV = rm
endif


################### TARGETS #####################
all: $(TARGET)

status :
	@echo Version Git-repository : $(VERSION_STRING)
	@echo 
	@echo $(OBJ_COMMON)
	@echo $(SRC_COMMON)
	@echo $(SRC_KIELFLOW)


bin/kielflow : $(OBJ_COMMON) $(OBJ_KIELFLOW) 
	rm -f bin/kielflow
	$(CXX) $(CXX_FLAGS) $(DEF_FLAGS) $(INC) $(OBJ_COMMON) $(OBJ_KIELFLOW) -o ./bin/kielflow $(LIB)

bin/frontend : DEF_FLAGS += -D__FRONTEND__
bin/frontend : $(OBJ_COMMON) $(OBJ_FRONTEND)
	rm -f bin/frontend
	$(CXX) $(CXX_FLAGS) $(DEF_FLAGS) $(INC) $(OBJ_COMMON) $(OBJ_FRONTEND) -o ./bin/frontend $(LIB)

bin/backend : $(OBJ_COMMON) $(OBJ_BACKEND)
	rm -f bin/backend
	$(CXX) $(CXX_FLAGS) $(DEF_FLAGS) $(INC) $(OBJ_COMMON) $(OBJ_BACKEND) -o ./bin/backend $(LIB)


plot4 :
	./my_plots.sh -p splot_potential ./data/splot_data\ fields.dat
	./my_plots.sh -p splot_density ./data/splot_data\ fields.dat


# Objekt-Dateien erzeugen
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.$(SRC_EXT)
	@mkdir -p $(BUILD_DIR)/common
	@mkdir -p $(BUILD_DIR)/kielflow
	@mkdir -p $(BUILD_DIR)/backend
	@mkdir -p $(BUILD_DIR)/frontend
	$(CXX) $(CXX_FLAGS) $(DEF_FLAGS) $(INC) -c -o $@ $< -fopenmp -lm


windows :  $(OBJS)
	$(CXX) -o $(TARGET).exe $(OBJS) ./lib/libfftw3-3.dll


run/% : bin/%
	konsole --noclose -e ./$<

run/backend : bin/backend
	./$<
	$(MAKE) plot4
	@echo ""
	@echo "done"

run : $(EXC)
	@echo "target done"


run_WIN32 :
	start cmd.exe @cmd /k "call bin\main"


run_LINUX :
	konsole --noclose -e $(TARGET)


.PHONY : clean
clean :
	$(RMV) -r -f $(BUILD_DIR)/*

tools :
	g++ ./src/tools/MinMaxFinder.cpp -o ./bin/MinMaxFinder

-PHONY : konsole
konsole :
	konsole --noclose
