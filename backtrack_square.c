#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int ipow(int k, int n)
{
    int i;
    int a = 1;
    for(i=0; i<n; ++i)
    {
        a *= k;
    }
    return a;
}

int facc(int k)
{
    int i;
    int a = 1;
    for(i=1; i<=k; ++i)
        a *= i;
    return a;
}

long constrcount(int N, int dim)
{
    int i;
    int flen;
    int fsize;
    int *field;
    int spos;
    int pos, npos;
    int *moves;
    long count;
    int wlen;
    int *dimleaves;
    int lastleave;
    int dimity;
    
    if(N == 1)
        return 1;
    
    flen = 2*N + 1;
    fsize = ipow(flen, dim);
    
    spos = 0;
    for(i=0; i<dim; ++i)
        spos += N * ipow(flen, i);
    
    moves = (int *)malloc(2*dim * sizeof(int));
    for(i = 0; i<dim*2; ++i)
        moves[i] = ipow(flen, i/2) * (2*(i%2)-1);
    
    field = (int *)malloc(fsize * sizeof(int));
    for(i=0; i<fsize; ++i)
        field[i] = -1;
    
    
    dimleaves = (int *)malloc((dim) * sizeof(int));
    dimleaves[0] = 0;
    dimity = 0;
    lastleave = 0;
    
    
    pos = spos;
    field[pos] = 2*dim-1;
    count = 0;
    wlen = 0;
    
    
    while(1)
    {
        if(field[pos] % (2*dim+1) != 2*dim)
        {
            npos = pos + moves[field[pos] % (2*dim+1)];
            if(field[npos] == -1)
            {
                if(wlen==N-1)
                {
                    ++count;
                    ++field[pos];
                }
                else
                {
                    ++wlen;
                    if(dimity != dim)
                    {
                        if(field[pos] % (2*dim+1) == 2*(dim - dimity) - 1)
                        {
                            lastleave = wlen-1;
                            dimleaves[dimity] = lastleave;
                            ++dimity;
                        }
                        if(dimity != dim)
                            field[npos] = (2*dim+1) * (field[pos] % (2*dim+1)) + ((dim-dimity)*2 - 1);
                        else
                            field[npos] = (2*dim+1) * (field[pos] % (2*dim+1));
                    }
                    else
                    {
                        field[npos] = (2*dim+1) * (field[pos] % (2*dim+1));
                    }
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
            if(wlen==1)
                break;
            --wlen;
            if(wlen == lastleave)
            {
                --dimity;
                //printf("%d\n", dimity);
                lastleave = dimleaves[dimity-1];
            }
            npos = pos - moves[field[pos] / (2*dim+1)];
            field[pos] = -1;
            ++field[npos];
            pos = npos;
        }
    }
    
    
    free(field);
    free(moves);
    free(dimleaves);
    
    
    return count;
} 


long countwalk(int N, int dim)
{
    long *counts;
    int i;
    long count;
    
    counts = (long *)malloc(dim * sizeof(long));
    for(i=0; i<dim; ++i)
        counts[i] = constrcount(N, i+1);
    for(i=dim-1; i>0; --i)
        counts[i] -= counts[i-1];
    count = 0;
    for(i=0; i<dim; ++i)
        count += counts[i] * (long)((facc(dim) / facc(dim-i-1)) * ipow(2, i+1));
    return count;
}


int main(int argc, char *argv[])
{
    int n;
    long count;
    clock_t start, stop;
    float timed;
    
    for(n=1; n<=30; ++n)
    {
        start = clock();
        count = countwalk(n,2);
        stop = clock();
        timed = (float)(stop - start) / CLOCKS_PER_SEC;
        printf("%d, %ld, %f,\n", n, count, (float)(stop - start) / CLOCKS_PER_SEC);
        
    }
}
