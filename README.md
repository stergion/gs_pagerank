# gs_pagerank

A C implementation of the PageRank algorithm with Gauss-Seidel iterations. There are 2 versions:

- The serial version *pagerank_gs.c*
- The parallel version using OpenMP *pagerank_gs_omp.c*

### Compiling and Running the programs

To compile and run the programs:

1. Download the two `.c` files and the **data** directory and save them in the same directory.
2. Extract the **datasets.rar** file
3. Compile the programs:
  * Serial: `gcc -O3 pagerank_gs.c -o pagerank_gs -lm`
  * Parallel: `gcc -O3 pagerank_gs_omp.c -o pagerank_gs_omp -lm -fopenmp`
4. Run the programs:
  * For the serial version provide as an argument the directory containing the *adj_list* and *nodes* files. E.g.
  `$ ./pagerank_gs ./data/datasets/__computational_complexity`
  * For the parallel version provide as the first argument the directory containing the *adj_list* and *nodes* files and as the second argument the number of CPU threads. E.g.
  `$ ./pagerank_gs_omp ./data/datasets/__computational_complexity 4`
