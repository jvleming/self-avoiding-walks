// counts the number of Self Avoiding Walks (SAW) of length n on the honeycomb lattice
// implementation of Zbarsky's method: https://doi.org/10.1088/1751-8121/ab52b0
// the underlying transfer matrix method is: https://doi.org/10.48550/arXiv.1309.6709 
// author: Jannes Vleming
// 13 March 2023

// usage multithreading:
// compile like usual:
// gcc zbarsky.c -lm -O3
// for example you want to run 3 threads, run in command line:
// ./a.out 3 0
// ./a.out 3 1
// ./a.out 3 2
// this can be done simultaneously
// later add up the results

// parameters:
#define saw 1 // else SAP
#define zbar 1 // else jensen
#define honey 1 // else square

#define prune 1 // every how many steps, 0 no pruning, 
#define file 0 // else print
#define maxn 40 // when to stop
#define prime 189812533 // size of hashmap, must be prime

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>


#define ncross 8 // number of ints to store all crossings, each int has 32 bits to store 16 crossings
#define bl 16 // number of crossings that can be stored in datatype

int npol; // number of terms in polinomial: length of SAW
int intdiv; // 2^32 % prime
int numthreads; // number of threads
int thread; // thread index

// signature
struct sig
{
    struct sig *hashNext; // next in linked list in hashmap
    struct sig *iterNext; // next in linked list to cylce through
    int touch; // which sides has the SAW touched? 1 left, 2 up, 4 down, 8 right 
    unsigned int *cross; // crossings, 00=empty, 01=free, 10=upper, 11=lower
    unsigned long *poly; // generating polynomial
    int len; // number of (empty/nonempty) crossings
    int hash;
};

struct sig **hashMap;

// allocates memory for signature
struct sig *getSig()
{
    struct sig *sig;
    if (!(sig = malloc(sizeof(struct sig) + ncross * sizeof(unsigned int) + npol * sizeof(unsigned long))))
        printf("memory allocation failed\n");
    sig->cross = (unsigned int *)(sig + 1);
    sig->poly = (unsigned long *)(sig->cross + ncross);
    
    sig->hashNext = 0;
    sig->iterNext = 0;
    return sig;
}

// creates empty signature
struct sig *zeroSig(int len)
{
    int i;
    struct sig *sig = getSig();
    for (i = 0; i < ncross; ++i)
        sig->cross[i] = 0;
    for (i = 1; i < npol; ++i)
        sig->poly[i] = 0;
    sig->poly[0] = 1;
    sig->touch = 0;
    sig->len = len;
    return sig;
}

// copies signature, shifts if len increased, shifts polynomial
struct sig *copySig(struct sig *a, int len, int shift)
{
    int i;
    struct sig *out = getSig();
    if (len > a->len) // shift crossings
    {
        for (i = 0; i < ncross; ++i)
        {
            out->cross[i] = a->cross[i] << 2;
            if (i)
                out->cross[i] |= a->cross[i-1] >> (bl - 1) * 2;
        }
    }
    else // just copy crossings
    {
        for (i = 0; i < ncross; ++i)
            out->cross[i] = a->cross[i];
    }
    for (i = 0; i < shift; ++i) // shifts polyonial by shift
        out->poly[i] = 0;
    for (i = shift; i < npol; ++i)
        out->poly[i] = a->poly[i - shift];
    out->touch = a->touch;
    out->len = len;
    return out;
}

// calculates hash of signature
void hash(struct sig *a)
{
    int i;
    long l = a->touch;
    for (i = 0; i < ncross; ++i) // crossings are seen as one enormous integer and then divided
    {
        l = (l * intdiv + (a->cross[i] % prime)) % prime; 
    }
    a->hash = (int)l;
}

// stores signature in hashmap and linked list, adds together corresponding sigs
void store(struct sig *a)
{
    int i, check;
    
    hash(a);
    struct sig *cur = hashMap[a->hash];
    if (cur) // the hashmap is not empty at the hash of a
    {
        while(1)
        {
            // check if the signatures are the same
            check = a->touch == cur->touch;
            for (i = 0; i < ncross; ++i)
                check = check && a->cross[i] == cur->cross[i];
            if (check)
            {
                // add polynomial and discard a
                for (i = 0; i < npol; ++i)
                    cur->poly[i] += a->poly[i];
                free(a);
                return;
            }
            
            // look further in hashmap linked list
            if (!cur->hashNext)
            {
                cur->hashNext = a;
                break;
            }
            cur = cur->hashNext;
        }
    }
    else
    {
        hashMap[a->hash] = a; // add a to hashmap
    }
    a->iterNext = hashMap[-1]; // add a to cycle linked list
    hashMap[-1] = a;
}

// generates signatures continuing from sig a, edge pos-1,pos are affected 
unsigned long moveBd(struct sig *src, int pos, int freeDir)
{
    unsigned long count = 0;
    
    int i;
    int check, from;
    int dir, depth;
    struct sig *out;
    
    // what boundaries would a signature touch here?
    int ntouch = 1;
    if (!pos)
        ntouch |= 0b10;
    else if (pos == src->len - 1)
        ntouch |= 0b100;
    if (freeDir & 0b100) // right edge of rectangle, can close
        ntouch |= 0b1000;
    
    
    
    int nlen = src->len;
    int npos = pos;
    if (ntouch & 0b10) // first row increase length of signature
    {
        ++nlen;
        ++npos;
    }
    else if (ntouch & 0b100) // last row, decrease length of signature
        --nlen;
    int lnpos = npos - 1;
    
    // a is incoming edge from above, b incoming from left
    int a = 0;
    if (pos)
        a = src->cross[(pos - 1) / bl] >> 2 * ((pos - 1) % bl) & 0b11;
    int b = src->cross[pos / bl] >> 2 * (pos % bl) & 0b11;
    
    int c;
    
    // if not started yet
    if (!src->touch)
    {
        // case: empty
        if (!(ntouch & 0b100))
        {
            out = zeroSig(nlen);
            store(out);
        }
        
        if (saw)
        {
            // case: free to right
            if (freeDir & 0b01)
            {
                out = zeroSig(nlen);
                out->cross[lnpos / bl] |= 1u << 2 * (lnpos % bl);
                out->touch = ntouch;
                store(out);
            }
            
            // case: free to below
            if (freeDir & 0b10)
            {
                out = zeroSig(nlen);
                out->cross[npos / bl] |= 1u << 2 * (npos % bl); 
                out->touch = ntouch;
                store(out);
            }
            
            // case: free in both directions
            if ((freeDir & 0b11) == 0b11)
            {
                out = zeroSig(nlen);
                out->cross[lnpos / bl] |= 1u << 2 * (lnpos % bl);
                out->cross[npos / bl] |= 1u << 2 * (npos % bl); 
                out->touch = ntouch;
                store(out);
            }
        }
        else
        {
            // in case of SAP: open loop
            if ((freeDir & 0b11) == 0b11)
            {
                out = zeroSig(nlen);
                out->cross[lnpos / bl] |= 2u << 2 * (lnpos % bl);
                out->cross[npos / bl] |= 3u << 2 * (npos % bl); 
                out->touch = ntouch;
                store(out);
            }
        }
    }
    else // we already started
    {
        if (!a && !b)
        {
            // case: continue empty
            out = copySig(src, nlen, 0);
            store(out);
            
            // case: catch free edge or loop
            
            // we first walk up in the signature, looking for possibilities
            dir = -1;
            i = pos - 2;
            depth = 0;
            while (1)
            {
                // if we reach top or get enclosed in loop, start again at i and go down
                if (dir == -1 && (i < 0 || depth < 0))
                {
                    dir = 1;
                    i = pos + 1;
                    depth = 0;
                }
                // if we reach bottom or get enclosed again, stop
                if (depth < 0 || i >= src->len)
                    break;
                
                
                from = i + npos - pos; // i in new (shifted) coordinates
                c = src->cross[i / bl] >> 2 * (i % bl) & 0b11; // edge at i
                if (!depth)
                {
                    // catch free
                    if (c == 1)
                    {
                        // catch right
                        if (freeDir & 0b01)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2u * (from % bl));
                            out->cross[lnpos / bl] |= (dir == -1 ? 3u : 2u) << 2 * (lnpos % bl);
                            out->cross[from / bl] |= (dir == -1 ? 2u : 3u) << 2 * (from % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                        
                        // catch below
                        if (freeDir & 0b10)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                            out->cross[npos / bl] |= (dir == -1 ? 3u : 2u) << 2 * (npos % bl);
                            out->cross[from / bl] |= (dir == -1 ? 2u : 3u) << 2 * (from % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                        
                        // catch and release
                        if((freeDir & 0b11) == 0b11)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                            out->cross[npos / bl] |= 1u << 2 * (npos % bl);
                            out->cross[lnpos / bl] |= (dir == -1 ? 3u : 2u) << 2 * (lnpos % bl);
                            out->cross[from / bl] |= (dir == -1 ? 2u : 3u) << 2 * (from % bl);
                            out->touch |= ntouch;
                            store(out);
                            
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                            out->cross[npos / bl] |= (dir == -1 ? 3u : 2u) << 2 * (npos % bl);
                            out->cross[lnpos / bl] |= 1u << 2 * (lnpos % bl);
                            out->cross[from / bl] |= (dir == -1 ? 2u : 3u) << 2 * (from % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                    }
                    
                    if((freeDir & 0b11) == 0b11)
                    {
                        // catch loop middle
                        if (dir == -1 && c == 2)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[lnpos / bl] |= 3u << 2 * (lnpos % bl);
                            out->cross[npos / bl] |= 2u << 2 * (npos % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                        
                        // catch loop from below
                        if (dir == -1 && c == 3)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                            out->cross[from / bl] |= 2u << 2 * (from % bl);
                            out->cross[lnpos / bl] |= 3u << 2 * (lnpos % bl);
                            out->cross[npos / bl] |= 3u << 2 * (npos % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                        
                        // catch loop from above
                        if (dir == 1 && c == 2)
                        {
                            out = copySig(src, nlen, 0);
                            out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                            out->cross[from / bl] |= 3u << 2 * (from % bl);
                            out->cross[lnpos / bl] |= 2u << 2 * (lnpos % bl);
                            out->cross[npos / bl] |= 2u << 2 * (npos % bl);
                            out->touch |= ntouch;
                            store(out);
                        }
                    }
                }
                
                if (c == 2)
                    depth += dir;
                if (c == 3)
                    depth -= dir;
                i += dir;
            }
            
        }
        else if (a && !b || !a && b)
        {
            // case: stop free edge
            if ((a | b) == 1)
            {
                from = a ? pos - 1 : pos;
                check = 1;
                for (i = 0; i < src->len; ++i) // check if otherwise empty
                {
                    if (i != from && src->cross[i / bl] & 3u << 2 * (i % bl))
                    {
                        check = 0;
                        break;
                    }
                }
                if (check)
                {
                    if ((src->touch | ntouch) == 0b1111) // signature can be completed
                    {
                        count += src->poly[npol - 1];
                    }
                }
                else // end free edge, signature not empty
                {
                    out = copySig(src, nlen, 1);
                    from += npos - pos;
                    out->cross[from / bl] &= ~(3u << 2 * (from % bl));
                    out->touch |= ntouch;
                    store(out);
                }
            }
            
            // case: continue edge left
            if (freeDir & 0b01)
            {
                out = copySig(src, nlen, 1);
                out->cross[lnpos / bl] &= ~(3u << 2 * (lnpos % bl));
                out->cross[npos / bl] &= ~(3u << 2 * (npos % bl));
                out->cross[lnpos / bl] |= (a | b) << 2 * (lnpos % bl);
                out->touch |= ntouch;
                store(out);
            }
            
            // case: continue edge below
            if (freeDir & 0b10)
            {
                out = copySig(src, nlen, 1);
                out->cross[lnpos / bl] &= ~(3u << 2 * (lnpos % bl));
                out->cross[npos / bl] &= ~(3u << 2 * (npos % bl));
                out->cross[npos / bl] |= (a | b) << 2 * (npos % bl);
                out->touch |= ntouch;
                store(out);
            }
        }
        else if (a == 2 && b == 3)
        {
            // case: close loop
            check = 1; // check otherwise empty
            for (i = 0; i < src->len; ++i)
            {
                if (i != pos && i != pos - 1 && src->cross[i / bl] & 3u << 2 * (i % bl)) 
                {
                    check = 0;
                    break;
                }
            }
            if (check) // we close loop and signature
            {
                if ((src->touch | ntouch) == 0b1111)
                {
                    count += src->poly[npol - 2];
                    
                }
            }
            else // we close the loop, but signature is nonempty
            {
                out = copySig(src, nlen, 2);
                out->cross[lnpos / bl] &= ~(3u << 2 * (lnpos % bl));
                out->cross[npos / bl] &= ~(3u << 2 * (npos % bl));
                out->touch |= ntouch;
                store(out);
            }
        }
    }
    return count;
}

// returns 1 if signature cannot contribute to count, otherwise 0
int allPrune (struct sig *a, int pat, int col, int ncol, int w, int row, int site)
{
    int i, j, d, minn, cr;
    int delx, dely, complete, touch, totcomplete, mintouch;
    int d1, d2, d3;
    
    // function is called after boundary has been moved, thus update site, row
    ++site;
    if (site == ncol - col)
    {
        site = 0;
        ++row;
    }
    
    // not started yet, can go
    if (!a->touch)
        return 0;
    
    // minn is first nonempty element of polynomial
    for (minn = 0; minn < npol && !a->poly[minn]; ++minn);
    if (minn == npol)
        return 1; // empty polynomial, no contribution
    
    // get some memory for 3 arrays
    static int arlen = 0;
    static int *depth = 0; // how many loops/ free edges around this point
    static int *x = 0; // x coordinate of point to the right/below of crossing
    static int *y = 0; // y coordinate ...
    if (arlen < a->len)
    {
        arlen = a->len;
        free(depth);
        if (!(depth = malloc(3 * arlen * sizeof(int))))
            printf("memory allocation failed\n");
        x = depth + arlen;
        y = x + arlen;
    }
    
    // assign depth, x, y at every point in boundary
    d = 0;
    for (i = 0; i < a->len; ++i)
    {
        cr = a->cross[i / bl] >> 2 * (i % bl) & 0b11;
        if (!cr)
            depth[i] = d;
        else if (cr == 2)
            depth[i] = ++d;
        else if (cr == 3)
            depth[i] = --d;
        else
            depth[i] = d + 1;
        
        if (row)
        {
            if (i < row)
            {
                x[i] = ncol + 1;
                y[i] = i;
            }
            else if (i < row + ncol - col - site)
            {
                x[i] = ncol - (i - row);
                y[i] = row;
            }
            else if (i == row + ncol - col - site)
            {
                x[i] = col + site + 1;
                y[i] = row;
            }
            else if ( i < row + ncol - col + 1)
            {
                x[i] = ncol - (i - row - 1);
                y[i] = row + 1;
            }
            else
            {
                x[i] = col + 1;
                y[i] = i - (ncol - col);
            }
        }
        else
        {
            if (i < 1)
            {
                x[i] = col + site + 1;
                y[i] = 0;
            }
            else if (i < site + 1)
            {
                x[i] = col + site + 1 - i;
                y[i] = 1;
            }
            else
            {
                x[i] = col + 1;
                y[i] = i - site;
            }
        }
    }
    
    totcomplete = 0; // number of edges needed to complete loops, +1 per free edge
    mintouch = a->touch & 0b1000 ? 0 : npol; // number of aditional edges needed to touch right border
    for (i = 0; i < a->len; ++i)
    {
        cr = a->cross[i / bl] >> 2 * (i % bl) & 0b11;
        
        
        if (cr == 1)
        {
            ++totcomplete;
            if (mintouch && depth[i] == 1)
            {
                touch = w - 1 - x[i];
                if (honey && touch && !(pat / 2))
                    touch += touch - 1; // for the zig-zag
                if (touch < mintouch)
                    mintouch = touch;
            }
        }
        
        if (cr == 2)
        {
            d1 = depth[i];
            d2 = depth[i];
            d3 = depth[i];
            for (j = i + 1; depth[j] >= depth[i]; ++j)
            {
                if (x[j] == ncol + 1)
                {
                    if (depth[j] > d1)
                        d1 = depth[j];
                }
                else if (x[j] > col + 1)
                {
                    if (depth[j] > d2)
                        d2 = depth[j];
                }
                else
                {
                    if (depth[j] > d3)
                        d3 = depth[j];
                }
            }
            d1 -= depth[i];
            d2 -= depth[i];
            d3 -= depth[i];
            if (honey && !(pat / 2))
            {
                d2 *= 2;
            }
            else if (honey)
            {
                d1 *= 2;
                d3 *= 2;
            }
            
            if(ncol + d1 >= w || col + d3 >= w)
                return 1;
            
            delx = d1 + d3 + abs((x[i] + d1) - (x[j] + d3)); // number of horizontal edges at least
            dely = d2 ? d2 + row - y[i] + abs(y[j] - (row + d2)): y[j] - y[i]; // vertical edges
            
            // zig zag on honey
            if (honey && !(pat / 2))
                dely = dely > delx - 1 ? dely : delx - 1;
            else if (honey)
                delx = delx > dely - 1 ? delx : dely - 1;
            
            complete = delx + dely + 2;
            totcomplete += complete;
            
            // calculate extra edges to touch border
            if (depth[i] == 1)
            {
                delx = w - 1 - x[i] + w - 1 - x[j];
                dely = d2 ? d2 + row - y[i] + abs(y[j] - (row + d2)): y[j] - y[i];
                
                if (honey && !(pat / 2))
                    dely = dely > delx - 1 ? dely : delx - 1;
                else if (honey)
                    delx = delx > dely - 1 ? delx : dely - 1;
                touch = delx + dely + 2 - complete;
                if (mintouch > touch)
                    mintouch = touch;
            }
        }
    }
    if (minn + totcomplete + mintouch > npol)
        return 1;
    return 0;
}

int maxsig; // maximal number of simultaneous signatures, for every n


// counts number of SAW of length n in rectangle,
// w width, h height, k as in zbarsky paper
// patern: 0, 1 vertical, 2, 3, horizontal,
// skip: bits are columns to skip
unsigned long countRectSkip(int n, int w, int h, int pat, int k, int skip)
{
    int i, pos;
    int crossCount;
    int freeDir;
    int col, ncol, row, site;
    struct sig *start, *cur, *icur;
    int q = n / k;
    npol = n;
    int doprune;
    int sigcount;
    int iprune = 0;
    
    unsigned long count = 0;
    
    // other thread will take care of this rect
    if (numthreads != 1 && rand() % numthreads != thread)
        return 0;
    
    start = zeroSig(h); // begin with nothing
    col = -1; // current added column
    ncol = 0; // column to add, add columns from left to right
    while(1)
    {
        for (row = 0; row < h; ++row) // add row by row from up to down
        {
            for (site = 0; site < ncol - col; ++site) // at every row, add the ncol-col sites one by one
            {
                freeDir = 0b11; // 1: can go right, 2: can go down, 4: can close
                
                if (row == h - 1) // cannot go down if at bottom
                    freeDir &= 0b01;
                if (honey && (col + 1 + row + site + pat) % 2) // restrictions of the honeycomb lattice
                    freeDir &= pat / 2 ? 0b01 : 0b10;
                if (col + site == w - 2) // all the way to the right: cannot go further right, can close
                {
                    freeDir &= 0b10;
                    freeDir |= 0b100;
                }
                
                // for every signature, generate new signatures with moved boundary
                for (cur = start; cur; cur = cur->iterNext)
                {
                    count += moveBd(cur, row ? row + ncol - col - site : 0, freeDir);
                    sigcount++;
                }
                
                // cleanup
                sigcount = 0;
                for (cur = start; cur; cur = icur)
                {
                    icur = cur->iterNext;
                    free(cur); // delete all old signatures
                    ++sigcount;
                }
                start = hashMap[-1];
                hashMap[-1] = 0;
                for (cur = start; cur; cur = cur->iterNext)
                {
                     hashMap[cur->hash] = 0; // clean up hashmap
                     ++sigcount;
                }
                if (sigcount > maxsig)
                    maxsig = sigcount;
                
                // pruning
                if (prune && ++iprune == prune)
                {
                    iprune = 0;
                    icur = 0;
                    for (cur = start; cur;)
                    {
                        doprune = allPrune(cur, pat, col, ncol, w, row, site);
                        if (doprune) // delete signature
                        {
                            if (icur)
                            {
                                icur->iterNext = cur->iterNext;
                                free(cur);
                                cur = icur->iterNext;
                            }
                            else
                            {
                                start = cur->iterNext;
                                free(cur);
                                cur = start;
                            }
                        }
                        else
                        {
                            icur = cur;
                            cur = icur->iterNext;
                        }
                    }
                }
                
            }
            
            // throw signatures away if too many crossings
            if(k != 1 && (ncol || skip & 1))
            {
                icur = 0;
                for (cur = start; cur;)
                {
                    crossCount = 0;
                    for (i = 0; i < row + 1; ++i)
                    {
                        if (cur->cross[i / bl] >> 2 * (i % bl) & 0b11)
                        {
                            ++crossCount;
                        }
                    }
                    if (crossCount > q) // delete signature
                    {
                        if (icur)
                        {
                            icur->iterNext = cur->iterNext;
                            free(cur);
                            cur = icur->iterNext;
                        }
                        else
                        {
                            start = cur->iterNext;
                            free(cur);
                            cur = start;
                        }
                    }
                    else
                    {
                        icur = cur;
                        cur = icur->iterNext;
                    }
                }
            }
            
        }
        
        // select new column to add
        col = ncol;
        if (++ncol == w)
            break;
        while (ncol != w - 1 && !(skip >> ncol % k & 1))
            ++ncol;

    }
    return count;
}

// counts the number of SAW of length in w*h rect, pattern pat
unsigned long countRect(int n, int w, int h, int pat)
{
    int i;
    int skip, pm;
    int k = 1;
    if (zbar)
        k = (int)(sqrt(n * log(n)) + 0.5);
    
    unsigned long count = 0;
    
    // cycle through all nonempty subsets of columns
    for (skip = 1; skip < 1 << k; ++skip)
    {
        pm = -1;
        for (i = 0; i < k; ++i)
        {
            pm *= skip >> i & 1 ? -1 : 1; // calculate parity of number of not skipped columns
        }
        count += countRectSkip(n, w, h, pat, k, skip) * pm;
    }
    return count;
}

// counts the number of SAW of length n
unsigned long count(int n)
{
    int iw, ih, pat;
    unsigned long count = 0;
    
    if (honey)
    {
        for (iw = 1; iw <= (n + 1) / 2; ++iw)
        {
            for (ih = 1; ih <= n - iw; ++ih)
            {
                if (iw % 2 && !(ih % 2))
                {  
                    if (iw >= ih)
                    {
                        count += countRect(n, iw + 1, ih + 1, 0);
                        count += countRect(n, iw + 1, ih + 1, 1);
                    }
                    else
                    {
                        count += countRect(n, ih + 1, iw + 1, 2);
                        count += countRect(n, ih + 1, iw + 1, 3);
                    }
                }
                else
                {
                    if (iw >= ih)
                        count += 2 * countRect(n, iw + 1, ih + 1, 0);
                    else
                        count += 2 * countRect(n, ih + 1, iw + 1, 2);
                }
            }
        }
    }
    else // square grid
    {
        for (iw = 1; iw <= n; ++iw)
        {
            for (ih = 1; iw + ih <= n && ih <= iw; ++ih)
            {
                count += countRect(n, iw + 1, ih + 1, -1) * (iw == ih ? 1 : 2);
            }
        }
    }
    
    if (saw)
    {
        if (numthreads == 1 || thread == 0)
            count += 2; // this is for rectangles of area 0
        if (!honey)
            count *= 2; // symmetry
    }
    return count;
}

int main (int argc, char *argv[])
{
    int n;
    clock_t start, stop;
    unsigned long ct;
    int i;
    FILE *fp;
    char filename[100];
    char s[1000];
    
    numthreads = 1;
    if(argc > 2)
    {
        numthreads = atoi(argv[1]);
        thread = atoi(argv[2]);
    }
    if (file)
    {
        mkdir("out", S_IRWXU);
        sprintf(filename, "out/%d of %d.txt", thread, numthreads);
        fp = fopen(filename, "w");
    }
    
    if (!(hashMap = malloc((prime + 1) * sizeof(struct sig *))))
        printf("memory allocation failed\n");
    for (i = 0; i < prime; ++i)
        hashMap[i] = 0;
    ++hashMap;
    intdiv = 1;
    for (i = 0; i < bl * 2; ++i)
        intdiv = (intdiv << 1) % prime;
    
    sprintf(s, "%s on %s, %s, thread %d of %d\n", saw ? "SAW" : "SAP",
                honey ? "honeycomb" : "square", zbar ? "Zbarsky" : "Jensen", thread, numthreads);
    if (file)
    {
        fprintf(fp, s);
        fflush(fp);
    }
    else
        printf(s);
    
    for (n = 2; n <= maxn; n += saw ? 1 : 2)
    {
        maxsig = 0;
        start = clock();
        ct = count(n);
        stop = clock();
        sprintf(s, "%d, %ld, %f, %d,\n", n, ct, (float)(stop - start) / CLOCKS_PER_SEC, maxsig);
        if (file)
        {
            fprintf(fp, s);
            fflush(fp);
        }
        else
            printf(s);
    }
    
    if (file)
        fclose(fp);
    free(hashMap - 1);
}
