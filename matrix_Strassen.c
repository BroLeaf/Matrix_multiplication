#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define thread_num 1

int r, c;
int **A, **B, **C;

static double diff_in_second (struct timespec t1, struct timespec t2)
{
    struct timespec diff;
    if (t2.tv_nsec-t1.tv_nsec < 0) {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec - 1;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec + 1000000000;
    } else {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec;
    }
    return (diff.tv_sec + diff.tv_nsec / 1000000000.0);
}

void get_matrix ()
{
    int i, j;
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            scanf("%d", &A[i][j]);
    scanf("%d%d", &c, &r);
    for (i = 0; i < c; i++)
        for (j = 0; j < r; j++)
            scanf("%d", &B[i][j]);
}

void display (int *A[])
{
    int i, j;
    for (i = 0; i < r; i++, printf("\n"))
        for (j = 0; j < c; j++)
            printf("%5d ", A[i][j]);
}


void *matrix_mul_thread (void *rank_arg)
{
    long rank = (long)rank_arg;
    for (int i = rank * r / thread_num; i < (rank+1) * r / thread_num; i++)
        for (int j = 0; j < c; j++)
            for (int k = 0; k < r; k++)
                C[i][j] += A[i][k] * B[k][j]; 
    pthread_exit(0);
    return NULL;
}

void mul (int *A[], int *B[], int *C[])
{
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                C[i][j] += A[i][k] * B[k][j];
}

void matrix_sub (int N, int *A[], int *B[], int *C[])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = A[i][j] - B[i][j];
}

void matrix_add (int N, int *A[], int *B[], int *C[])
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            C[i][j] = A[i][j] + B[i][j];
}

void Strassen (int N, int *A[], int *B[], int *C[])
{
    
    /* divide origin matrix into four pieces
            A11 | A12
            A21 | A22   
    */
    if (N == 2)
    {    
        mul(A, B, C);
    }
    else    
    {
        int **A11, **A12, **A21, **A22;
    int **B11, **B12, **B21, **B22;
    int **C11, **C12, **C21, **C22;
    int **P1, **P2, **P3, **P4, **P5, **P6, **P7;
    int **AA, **BB;
    A11 = (int **)malloc(N/2 * sizeof(int *));A12 = (int **)malloc(N/2 * sizeof(int *));A21 = (int **)malloc(N/2 * sizeof(int *));A22 = (int **)malloc(N/2 * sizeof(int *));
    B11 = (int **)malloc(N/2 * sizeof(int *));B12 = (int **)malloc(N/2 * sizeof(int *));B21 = (int **)malloc(N/2 * sizeof(int *));B22 = (int **)malloc(N/2 * sizeof(int *));
    C11 = (int **)malloc(N/2 * sizeof(int *));C12 = (int **)malloc(N/2 * sizeof(int *));C21 = (int **)malloc(N/2 * sizeof(int *));C22 = (int **)malloc(N/2 * sizeof(int *));
    P1 = (int **)malloc(N/2 * sizeof(int *));P2 = (int **)malloc(N/2 * sizeof(int *));P3 = (int **)malloc(N/2 * sizeof(int *));P4 = (int **)malloc(N/2 * sizeof(int *));
    P5 = (int **)malloc(N/2 * sizeof(int *));P6 = (int **)malloc(N/2 * sizeof(int *));P7 = (int **)malloc(N/2 * sizeof(int *));
    AA = (int **)malloc(N/2 * sizeof(int *));
    BB = (int **)malloc(N/2 * sizeof(int *));
    for (int i = 0; i < N/2; i++)
    {
        A11[i] = (int *)malloc(N/2 * sizeof(int));
        A12[i] = (int *)malloc(N/2 * sizeof(int));
        A21[i] = (int *)malloc(N/2 * sizeof(int));
        A22[i] = (int *)malloc(N/2 * sizeof(int));
        B11[i] = (int *)malloc(N/2 * sizeof(int));
        B12[i] = (int *)malloc(N/2 * sizeof(int));
        B21[i] = (int *)malloc(N/2 * sizeof(int));
        B22[i] = (int *)malloc(N/2 * sizeof(int));
        C11[i] = (int *)malloc(N/2 * sizeof(int));
        C12[i] = (int *)malloc(N/2 * sizeof(int));
        C21[i] = (int *)malloc(N/2 * sizeof(int));
        C22[i] = (int *)malloc(N/2 * sizeof(int));
        P1[i] = (int *)malloc(N/2 * sizeof(int));
        P2[i] = (int *)malloc(N/2 * sizeof(int));
        P3[i] = (int *)malloc(N/2 * sizeof(int));
        P4[i] = (int *)malloc(N/2 * sizeof(int));
        P5[i] = (int *)malloc(N/2 * sizeof(int));
        P6[i] = (int *)malloc(N/2 * sizeof(int));
        P7[i] = (int *)malloc(N/2 * sizeof(int));
        AA[i] = (int *)malloc(N/2 * sizeof(int));
        BB[i] = (int *)malloc(N/2 * sizeof(int));
    }
        for (int i = 0; i < N/2; i++)
            for (int j = 0; j < N/2; j++)
            {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j+N/2];
                A21[i][j] = A[i+N/2][j];
                A22[i][j] = A[i+N/2][j+N/2];
                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j+N/2];
                B21[i][j] = B[i+N/2][j];
                B22[i][j] = B[i+N/2][j+N/2];
            }
    
        
        // P1 = A11 x (B12 - B22)
        matrix_sub (N/2, B12, B22, BB);
        Strassen (N/2, A11, BB, P1);
        // P2 = (A11 + A12) x B22
        matrix_add (N/2, A11, A12, AA);
        Strassen (N/2, AA, B22, P2);
        // P3 = (A21 + A22) x B11
        matrix_add(N/2, A21, A22, AA);
        Strassen (N/2, AA, B11, P3);
        // P4 = A22 x (B21 - B11)
        matrix_sub (N/2, B21, B11, BB);
        Strassen (N/2, A22, BB, P4);
        // P5 = (A11 + A22) x (B11 + B22)
        matrix_add (N/2, A11, A22, AA);
        matrix_add (N/2, B11, B22, BB);
        Strassen (N/2, AA, BB, P5);
        // P6 = (A12 - A22) x (B21 + B22)
        matrix_sub (N/2, A12, A22, AA);
        matrix_add (N/2, B21, B22, BB);
        Strassen (N/2, AA, BB, P6);
        // P7 = (A11 - A21) x (B11 + B12)
        matrix_sub (N/2, A11, A21, AA);
        matrix_add (N/2, B11, B12, BB);
        Strassen (N/2, AA, BB, P7);
        // C11 = P6 + P5 + P4 - P2;
        matrix_add (N/2, P6, P5, AA);
        matrix_sub (N/2, P4, P2, BB);
        matrix_add (N/2, AA, BB, C11);
        // C12 = P1 + P2
        matrix_add (N/2, P1, P2, C12);
        // C21 = P3 + P4
        matrix_add (N/2, P3, P4, C21);
        // C22 = P1 + P5 - P3 - P7  
        matrix_add (N/2, P1, P5, AA);
        matrix_add (N/2, P3, P7, BB);
        matrix_sub (N/2, AA, BB, C22);
        // save the result
        for (int i = 0; i < N/2; i++)
            for (int j = 0; j < N/2; j++)
            {
                C[i][j] = C11[i][j];
                C[i][j + N/2] = C12[i][j];
                C[i + N/2][j] = C21[i][j];
                C[i + N/2][j + N/2] = C22[i][j];
            }   
        free(A11);free(A12);free(A21);free(A22);
        free(B11);free(B12);free(B21);free(B22);
        free(C11);free(C12);free(C21);free(C22);
        free(AA);free(BB);
        free(P1);free(P2);free(P3);free(P4);free(P5);free(P6);free(P7);
    }
    
}

int main()
{
    int i, j;
    struct timespec start, end;
    //pthread_t thread[thread_num];
    scanf("%d%d", &r, &c);
    /* declare 2D array of matrix */
    A = (int **)malloc(r * sizeof(int *));
    B = (int **)malloc(r * sizeof(int *));
    C = (int **)malloc(r * sizeof(int *));
    for (i = 0; i<r; i++)
    {
        A[i] = (int *)malloc(r*sizeof(int));
        B[i] = (int *)malloc(r*sizeof(int));
        C[i] = (int *)malloc(r*sizeof(int));
    }

    /* input matrix */
    get_matrix();
    /* initialize */
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            C[i][j] = 0;

    // Strassen mul
    clock_gettime(CLOCK_REALTIME, &start);
    Strassen (r, A, B, C);
    clock_gettime(CLOCK_REALTIME, &end);

    printf("execution time of Strassen : %lf sec\n", diff_in_second(start, end));
    //display (C);

    return 0;
}

