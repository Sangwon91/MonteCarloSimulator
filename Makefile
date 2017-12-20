SRCS = $(wildcard *.cpp)

all: mc.x

mc.x:$(SRCS)
	gfortran -o $(@) -ffast-math -O2 -fimplicit-none -ffixed-form -ffree-form -fall-intrinsics  -Wtabs main.f90

clean:
	rm *.mod *.x
