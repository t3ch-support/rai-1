BASE   = ../..
NAME   = $(shell basename `pwd`)
OUTPUT = lib$(NAME).so

OPTIM = fast_debug

DEPEND = Core Algo Kin Gui Geo KOMO Optim

GL = 1

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/_make/generic.mk
