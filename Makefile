CC = mpicc
FLAGS = -c
LINK_FLAGS =
OBJS = main.o column.o gauss_jordan.o mpi_wrappers.o
EXECUTABLE_NAME = gj

# Compile
all: $(OBJS)
	$(CC) $(OBJS) $(LINK_FLAGS) -o $(EXECUTABLE_NAME)
	rm -f $(OBJS)

main.o: main.c
	$(CC) $(FLAGS) main.c

column.o: column.c
	$(CC) $(FLAGS) column.c

gauss_jordan.o: gauss_jordan.c
	$(CC) $(FLAGS) gauss_jordan.c

mpi_wrappers.o: mpi_wrappers.c
	$(CC) $(FLAGS) mpi_wrappers.c

# Clean-up
clean:
	rm -f $(EXECUTABLE_NAME)
clean-all:
	rm -rf $(EXECUTABLE_NAME) $(OUTPUT)

run:
	mpiexec -n $(PN) ./$(EXECUTABLE_NAME) $(N) $(GROUP)

# Usage
help:
	@:
		$(info Run with: make run PN=NUMBER_OF_MPI_PROCESSES N=DIMENSION GROUP=FLAG(0 or 1))
