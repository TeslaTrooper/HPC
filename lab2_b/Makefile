#default build: gameoflife
#all: gameoflife
all: gameoflife

#gameoflife build depends on gameoflife.c
gameoflife: gameoflife.c
	gcc -O3 -std=gnu99 -fopenmp -o gameoflife gameoflife.c

#heat_equation build depends on heat_equation.c
heat_equation: heat_equation.c
	gcc -O3 -lpng -std=gnu99 -fopenmp -o heat_equation heat_equation.c
	rm -r heq
