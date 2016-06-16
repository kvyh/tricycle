all: clean polyfit

polyfit:  polyfit.c
	gcc -g -Wall -O2 -I/usr/include -I/path/to/local/clibrary/include -L/path/to/local/clibrary/lib polyfit.c -lgslcblas -lgsl -lm -o polyfit

clean:
	rm -f polyfit
