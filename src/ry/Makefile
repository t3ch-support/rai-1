BASE   = ../..
NAME   = $(shell basename `pwd`)
EXT    = $(shell python3-config --extension-suffix)
OUTPUT = lib$(NAME).so #$(EXT)

DEPEND = Core Algo Optim Gui Geo Kin KOMO PathAlgos

PYBIND = 1
PHYSX = 1

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/_make/generic.mk
