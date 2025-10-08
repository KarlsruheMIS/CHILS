SHELL = /usr/bin/env bash

METIS_TO_JSON : main.c
	$(CC) -std=gnu17 -fopenmp -O3 -march=native -DNDEBUG -o $@ $^ 

.PHONY : clean
clean :
	rm -f METIS_TO_JSON