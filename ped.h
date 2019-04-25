#ifndef PED_H
#define PED_H

#include <stdio.h>
#include <stdint.h>

#include "map.h"
#include "mats.h"

/*! @typedef
  @abstract One line of PED information.
  @field  fid     Family ID name from PED.
  @field  iid     Individual ID name from PED.
  @field  pos     Pointer to pos_t positional data (def: map.h).
  @field  a       2xc or 1xc matrix of base information.
  @field  c       Length of pos array.
  @field  flag    Flag from shared derived file.

  @discussion PED lines are parsed in read_ped and output by ped_line(). After
  printing, they are either discarded, or kept in the peda_t struct (see below).
*/
typedef struct pedl_t {

  char *fid;
  char *iid;
  pos_t *pos;
  mats_m *a;
  int c;
  int flag;

} pedl_t;

/*! @typedef
  @abstract Struct to store an array of pedl_t PED lines.
  @field  pedk    Pointer to array of pedl_t PED lines.
  @field  c       Length of pedk array.
  @field  flag    Flag to distinguish different uses of this data.

  @discussion Gets filled with PED information in read_ped() depending on which
  lines are selected for later analysis.
*/
typedef struct peda_t {
  pedl_t **pedk;
  int c;
  int flag;
} peda_t;

/*! @typedef
  @abstract Connecting name with flag for shared derived test.
  @field  flag    Flag set in the file defining the shared derived test.
  @field  fid     FID of the PED individual carrying that flag.
  @field  iid     IID of the PED individual carrying that flag.

  @discussion Information for the shared derived test. The flag in f is passed
  on to the pedl_t.flag if fid and iid correspond to those in pedl_t.
  Once the pedl_t.flag is set, the PED line is not discarded but instead kept in
  pedl_t **pedk, which is an array of pointers to pedl_t PED lines. With that,
  all information relevant for the shared derived test is contained within that array.

          flags for shared derived test: 1 to define individual as ancestral
                                         2 to define individual as a test sample
                                         3 to define individual as a query sample
 */
typedef struct ad_t {
  int  flag;
  char *fid;
  char *iid;
} ad_t;

/*! @function
  @abstract Print a PED line from a pedl_t struct.
  @param  line    pedl_t PED line struct to print.
  @param  a       Number of alleles to sample and print (1/2).
  @param  out     File handle of the PED output file.
  @return         Status indicator.

  @discussion The sampling of one allele if a == 1 is also done here. In this
  case, the genotype matrix will be reduced from two rows to one.
 */
int ped_line(pedl_t *line, int a, FILE *out);

/*! @function
  @abstract Free pedl_t PED line struct.
  @param  line    pedl_t struct to free.
 */
void free_pedl(pedl_t *line);

/*! @function
  @abstract Free ad_t struct array.
  @param  adl    Pointer to the array to free.
  @param  adc    Length of array in adl.
 */
void free_adl(ad_t *adl, int adc);

/*! @function
  @abstract Read PED file and related information.
  @param  fn      File name of input PED file.
  @param  l       Number of positions, length of array in map.
  @param  a       Number of alleles to sample (1/2).
  @param  out     Output PED file handle.
  @param  map     Array with position information (from map file).
  @param  adfn    File name for shared derived test (-t option).
  @param  peda    Pointer to struct where select PED lines get stored.
  @return         Status indicator -or- number of PED lines stored in pedk.

  @discussion Reads not only the PED file, but also the shared derived test
  file. Depending on that file, PED lines get stored in peda if applicable. Also
  stores allele information from PED panel across individuals in map for later
  use (see also pos_t documentation in map.h).
 */
int read_ped(char *fn, int l, int a, FILE *out, pos_t *map, char *adfn, peda_t *peda);

/*! @function
  @abstract Free the peda_t array struct.
  @param  peda    Struct with pedl_t array to be freed.
 */
void free_peda(peda_t *peda);

#endif
