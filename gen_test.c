#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char* argv[])
{
    srand(time(NULL));
    int r = atoi(argv[1]), c = atoi(argv[2]), i = 0, j = 0;
    printf("%d %d\n", r, c);
    for(i = 0; i < r; i++)
    {
        for(j = 0; j < c; j++)
            printf("%d ", (rand()%999 + 1));
        printf("\n");
    }
    printf("%d %d\n", c, r);
    for(i = 0; i < c; i++)
    {
        for(j = 0; j < r; j++)
            printf("%d ", (rand()%999 + 1));
        printf("\n");
    }
    
    return 0;
}
