#ifndef DERIV_H
#define DERIV_H

#include "ped.h"

/*
  Several items related to the shared derived allele test are defined within
  ped.h, where the data is read into the pedl_t **pedk struct. The flags one can
  set in the file defining the individuals participating in the test are
  repeated here:

    flags for shared derived test: 1 to define individual as ancestral
                                   2 to define individual as a test sample
                                   3 to define individual as a query sample
 */

/*! @function
  @abstract Count informative and shared derived sites.
  @param  peda     Struct of array with saved PED lines to test (against).
  @param  i        ID number of site being investigated.
  @param  s        Base sampled, which is being tested.
  @param  infor    Array of informative site counts per test combination.
  @param  share    Array of shared derived site counts per test combination.
  @return          Status indicator.
 */
int count_deriv(peda_t *peda, int i, char s, int *infor, int *share);

/*! @function
  @abstract Print the header of the shared derived output file.
  @param  cdc      Number of combinations to be tested.
  @param  peda     Struct of array with saved PED lines.
  @param  adout    File handle of the adstat output file.
  @return          Status indicator.
 */
int print_deriv_header(int cdc, peda_t *peda, FILE *adout);

/*! @function
  @abstract Print the final shared and informative counts.
  @param  cdc      Number of combinations which were tested.
  @param  name     Name of the sample which was tested.
  @param  infor    Array of informative site counts per test combination.
  @param  share    Array of shared derived site counts per test combination.
  @param  adout    File handle of the adstat output file.
  @return          Status indicator.
 */
int print_deriv(int cdc, char *name, int *infor, int *share, FILE *adout);

/*! @function
  @abstract Print header and test PED lines, if applicable.
  @param  peda     Struct of array with saved PED lines.
  @param  l        Number of positions (from map file).
  @param  adout    File handle of the adstat output file.
  @return          Status indicator.
 */
int ped_deriv(peda_t *peda, int l, FILE *adout);

#endif
