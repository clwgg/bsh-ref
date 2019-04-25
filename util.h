#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>

/*! @typedef
  @abstract Struct containing all output file handles.
  @field  pfout     Output PED file (should always exist)
  @field  adfout    Output adstat file (only with option -t)
  @field  safout    Output sastat file (only with option -s)

  @discussion:
          adstat: Contains output of shared derived allele test
          sastat: Contains sampling statistics from bam file sampling

    Current format of sastat file:
    name | name/nrep | num pos | cov | base sampled | rel read pos | ctot | PED bases
 */
typedef struct outf_t {
  FILE *pfout;
  FILE *adfout;
  FILE *safout;
} outf_t;

/*! @function
  @abstract Number of rows in file.
  @param  fn    File name.
  @return       Number of rows.
 */
int nrows(char *fn);

/*! @function
  @abstract Number of columns in file.
  @param  fn        File name.
  @param  offset    Offset to substract from column count.
  @return           Number of columns - offset.
*/
int ncols(char *fn, int offset);

/*! @function
  @abstract Close file handles in the outf_t struct.
  @param  outf    outf_t struct to close file handles in.
 */
void close_outf(outf_t *outf);

/*! @function
  @abstract Number of positions in an int.
  @param  n    Integer input.
*/
int numPlaces(int n);

#endif
