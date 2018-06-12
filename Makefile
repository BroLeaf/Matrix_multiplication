CC=gcc
CFLAGS= -Wall -std=gnu99
TUNE= -O 

all: gen_test traditional Strassen

gen_test:
		@$(CC) $(CFLAGS) -o gen_test gen_test.c

traditional:
		@$(CC) $(TUNE) $(CFLAGS) -o matrix_mul matrix_mul.c -lpthread

Strassen:
		@$(CC) $(TUNE) $(CFLAGS) -o matrix_Strassen matrix_Strassen.c

clean:
		$(RM) gen_test matrix_mul matrix_Strassen
