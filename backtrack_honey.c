#include <stdio.h>
#include <stdlib.h>
#include <time.h>

long count(int N)
{
    int i;
    int width, height;
    int *field;
    int pos, npos;
    int wlen;
    long count;
    int *moves;
    
    if (N == 1)
        return 3;
    if (N == 2)
        return 6;
    
    width = 2 * N + 1;
    height = ((N + 1) / 2) * 2 + 1;
    pos = N + width * ((N + 1) / 2);
    
    field = (int *)malloc(width * height * sizeof(int));
    for (i = 0; i < width * height; ++i)
        field[i] = -1;
    
    moves = (int *)malloc(6 * sizeof(int));
    moves[0] = moves[3] = 1;
    moves[1] = moves[4] = -1;
    moves[2] = width;
    moves[5] = -width;
    
    field[pos] = 0;
    ++pos;
    field[pos] = 0;
    ++pos;
    field[pos] = 0;
    wlen = 2;
    count = 0;
    
    while (1)
    {
        if (field[pos] % 4 != 3)
        {
            npos = pos + moves[field[pos] % 4 + 3 * (pos % 2)];
            if(field[npos]==-1)
            {
                if(wlen==N-1)
                {
                    ++count;
                    ++field[pos];
                }
                else
                {
                    ++wlen;
                    field[npos] = 4 * (field[pos] % 4);
                    pos = npos;
                }
            }
            else
            {
                ++field[pos];
            }
        }
        else
        {
            if(wlen==2)
                break;
            --wlen;
            npos = pos - moves[field[pos] / 4 + 3 * ((pos+1) % 2)];
            field[pos] = -1; 
            ++field[npos];
            pos = npos;
        }
    }
    
    free(field);
    free(moves);
    
    return count * 6;
}

int main(int argc, char *argv[])
{
    int n;
    long lcount;
    clock_t start, stop;
    float timed;
    
    for(n=1; n<=40; ++n)
    {
        start = clock();
        lcount = count(n);
        stop = clock();
        timed = (float)(stop - start) / CLOCKS_PER_SEC;
        printf("%d, %ld, %f,\n", n, lcount, (float)(stop - start) / CLOCKS_PER_SEC);
        
    }
}
