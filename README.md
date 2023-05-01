# self-avoiding-walks
An algorithm to count the number of self-avoiding walks in sub-exponential time.

In zbarsky.c we exactly enumerate the number of self-avoiding walks on the honeycomb lattice, using Zbarsky's method, in sub-exponential time.
The code also works for counting self-avoiding polygons, the square lattice, and Jensen's more traditional transfer matrix method, for this change the parameters in the code.

The report Self_avoiding_walks.pdf explains the method in detail.

Also added are some results obtained by running the program on a supercomputer.
Furthermore backtrack_honey.c and backtrack_square are backtracking algorithms for the honeycombgrid and square grid respectively.
