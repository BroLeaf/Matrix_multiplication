CC=gcc
CFLAGS= -Wall -std=gnu99
TUNE= -O 
EXEC = parallel Strassen optimization

all: gen_test parallel Strassen Strassen_parallel

gen_test:
		$(CC) $(CFLAGS) -o gen_test gen_test.c

#sequence:
#		$(CC) $(TUNE) $(CFLAGS) -o matrix_seq matrix_seq.c

parallel:
		$(CC) $(TUNE) $(CFLAGS) -o matrix_mul matrix_mul.c -lpthread

Strassen:
		$(CC) $(TUNE) $(CFLAGS) -o matrix_Strassen matrix_Strassen.c

Strassen_parallel:
		$(CC) $(TUNE) $(CFLAGS) -fopenmp -o matrix_Strassen_parallel matrix_Strassen_parallel.c

optimization:
		$(CC) $(TUNE) $(CFLAGS) -o matrix_opt matrix_opt.c -lpthread

run: $(EXEC)
		@./matrix_mul < input.txt \
		&&./matrix_Strassen < input.txt &&./matrix_Strassen_parallel < input.txt \
		&& ./matrix_opt < input.txt

test: $(EXEC)
		perf stat --repeat 100 \
			-e cache-misses,cache-references,instructions,cycles \
			./matrix_opt < input.txt

clean:
		$(RM) gen_test matrix_seq matrix_mul matrix_Strassen matrix_Strassen_parallel matrix_opt
