/* 	$Id: util.h,v 1.2 2017/07/07 07:59:16 leonid Exp $	 */
#ifndef UTIL_H
#define UTIL_H

/************************************************************************ 
 * routines to calculate distance of a quantum LDPC code 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>
 * currently implemented: a version of the random window algorithm in 
 * A. Dumer, A. A. Kovalev, and L. P. Pryadko "Distance verification..."
 * in IEEE Trans. Inf. Th., vol. 63, p. 4675 (2017). 
 * doi: 10.1109/TIT.2017.2690381
 ************************************************************************/


#define SWAPINT(a,b) do{ int t=a; a=b; b=t; } while(0)

#define ERROR(fmt,...)                                                 \
  do{                                                                  \
    printf("#:[31;1m *** ERROR: " fmt " ***[0m \n",##__VA_ARGS__); \
    exit(-1);                                                          \
  }                                                                    \
  while(0)

/**
 * macros from nauty.h
 * SETWD(pos) gives the setword in which pos is located
 * SETBT(pos) gives the location of bit pos in a setword
 */
#define SETWD(pos) ((pos)>>6)
#define SETBT(pos) ((pos)&0x3F)
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */

#define FIRSTBIT(x) __builtin_ctzll(x) // number of trailing zeros 

#ifdef __POPCNT__ 

static inline int m4ri_bitcount(word w){
  return __builtin_popcountll(w);  
}

#else /* no __POPCNT__ */

#define MASK(c)    (((uint64_t)(-1)) / (__M4RI_TWOPOW(__M4RI_TWOPOW(c)) + 1))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (__M4RI_TWOPOW(c))) & MASK(c))

static inline int m4ri_bitcount(word w)  {
   uint64_t n = __M4RI_CONVERT_TO_UINT64_T(w);
   n = COUNT(n, 0);
   n = COUNT(n, 1);
   n = COUNT(n, 2);
   n = COUNT(n, 3);
   n = COUNT(n, 4);
   n = COUNT(n, 5);
   return (int)n;
}

#endif /* __POPCNT__ */


typedef struct{
  int css; /* 1: css, 0: non-css -- currently not supported */
  int steps; /* how many random decoding steps */
  int debug; /* debug information */ 
  char *finP, *finG; /* generators */
  int wmax; /* max cluster size to try */
  int wmin; /* min distance below which we are not interested at all */
  int seed;/* rng seed, set=0 for automatic */
  int dist; /* target distance of the code */
  int dist_max; /* distance actually checked */
  int dist_min; /* distance actually checked */
  int n0;  /* code length, =n for css, (n/2) for non-css */
  int n; /* actual n = matrix size*/
} params_t; 

/**
 * sparse binary matrix in compressed-row form (CSR, nz=-1) or 
 * List-Of-Pairs (nz pairs).
 * use mzp_compress() to convert from LOP to CSR. 
 */
typedef struct{    /*  */
  int rows ;	    /* number of rows */
  int cols ;	    /* number of columns */
  int nz ;	    /* # of entries in triplet matrix */
  int nzmax ;	    /* # allocated size */
  int *p ;	    /* row pointers (size rows+1) OR row indices */
  int *i ;	    /* col indices, size nzmax */
} csr_t ;

extern params_t prm;



#if defined(__cplusplus) && !defined (_MSC_VER)
extern "C" {
#endif
  
/**
 * total number of set bits in the mzd array A.  
 * Uses built-in popcount if available, otherwise fast macros from 
 * nauty.h (Nauty library by Brendan McKay)
 */ 
size_t mzd_weight(const mzd_t *A);

/**
* nextelement(set1,m,pos) = the position of the first element in set set1   
* which occupies a position greater than pos.  If no such element exists,   
* the value is -1.  pos can have any value less than n, including negative  
* values.                                                                   
*  
* near verbatim copy from naututil.c (Nauty library by Brendan McKay)
*/

int nextelement(word *set1, int m, int pos);

/**
 * Copy of mzd_gauss_delayed from mzd.c (m4ri package) except additionally 
 * returns the list of pivot columns 
 */
rci_t mzd_gauss_naive(mzd_t *M, mzp_t *q, int full);

/**
 * Convert CSR sparse binary matrix to MZD
 * allocate dst if needed (must be correct size or NULL)
 */
mzd_t *mzd_from_csr(mzd_t *dst, const csr_t *p);

/**
 * sparse-S by dense B multiplication
 * C=C+S*B; allocate C if needed.
 * if clear=1 set C=0 first 
 */
  mzd_t * csr_mzd_mul(mzd_t *C, const csr_t *S, const mzd_t *B, int clear);

/**
 * return uniformly distributed random number in the range [0,...,max-1] 
 * uses RAND internally 
 * \todo Replace by a better generator 
 */
int rand_uniform(const int max);

/**
 * replace pivot q with a random pivot of same length, 
 * *** note: LAPACK style pivot permutations! ***
 * return pointer to q.
 * input: perm -- existing permutation
 */ 
mzp_t * mzp_rand(mzp_t *q);

/**
 * print out the permutation (only needed under windows)
 */
void mzp_out(mzp_t const *p);

/**
 * apply pivot p to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p(mzp_t *q, const mzp_t *p,rci_t start);

/**
 * apply pivot p (transposed) to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p_trans(mzp_t *q, const mzp_t *p,const rci_t start);

/**
 * kill a CSR matrix 
 */
csr_t *csr_free(csr_t *p);

/**
 * initialize a CSR matrix 
 * check existing size and (re)allocate if  needded 
 */
csr_t *csr_init(csr_t *mat, int rows, int cols, int nzmax);

/**
 *  compress a sparse binary matrix from List-of-Pairs to CSR
 */ 
void csr_compress(csr_t *mat);

/**
 *  output a CSR matrix (List-of-Pairs or CSR) 
 */ 
void csr_out(const csr_t *mat);

/**
 * read sparse matrix into a (binary) CSR (all entries default to 1)
 * from sparse integer Matrix Market matrix file (.mmx)
 * (re)allocate mat if needed
 * use transpose=1 to transpose.
 */
csr_t *csr_mm_read(char *fin, csr_t *mat, int transpose);

/** 
 * Permute columns of a CSR matrix with permutation perm.
 */
csr_t *csr_apply_perm(csr_t *dst, csr_t *src, mzp_t *perm);



#if defined(__cplusplus) && !defined (_MSC_VER)
}
#endif


#endif /* UTIL_H */
