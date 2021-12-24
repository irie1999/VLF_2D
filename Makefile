OBJS = allocate_memory.o  initialize_surface_impedance.o  update_E.o \
	current_source.o main.o  update_H.o output.o  update_H_PML.o \
	initialize_conductivity.o  update_D.o initialize_pml.o  \
	update_D_PML.o suffix.o input.o

main: $(OBJS)
	g++ -o $@  $(OBJS)

%.o: %.cpp fdtd2d.h
	g++ -c $< -Wall -O3 -I.
