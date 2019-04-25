#ifndef BAM_H
#define BAM_H

#include <stdint.h>

#include "map.h"
#include "ped.h"
#include "util.h"
#include "deriv.h"

#include "htslib/htslib/sam.h"

/*! @typedef
  @abstract Base from bam at variable site.
  @field  ctot    Did the read have terminal CtoT at 5p (1), 3p (2) or both (3).
  @field  base    Base sequence information (translate using t_seq()).
  @field  pos     Position of base in read.
  @field  qlen    Length of the read.

  @discussion struct filled in read_pos(). sample_pos() samples a base_t at a
   given snp_t.
 */
typedef struct base_t {
  uint8_t ctot;
  uint8_t base;
  uint16_t pos;
  uint16_t qlen;
} base_t;

/*! @typedef
  @abstract Position of reference variability from bam.
  @field  bases    Array of base_t bases at this site.
  @field  c        Length of 'bases' array.
  @field  ped      Reference information from PED file.

  @discussion struct filled mostly in read_pos(). ped field filled in
  read_bases() with information gathered in read_ped().
 */
typedef struct snp_t {
  base_t *bases;
  int c;
  char ped;
} snp_t;

/*! @typedef
  @abstract Information related to bam files.
  @field  in           bam file handle.
  @field  h            bam header file handle.
  @field  idx          bam index file handle.
  @field  itr          Iterator targeting variable positions (from map file).
  @field  min_mapQ     Minimum mapping quality option.
  @field  name         Sample name parsed from blist input.
  @field  snpbases     Array of snp_t with length equal to map file line count.
  @field  ctot_flag    Option flag whether to filter by terminal CtoTs.

  @discussion
          ctot_flag: 0 for no filtering
                     1 for one end filtering
                     2 for both ends filtering
 */
typedef struct aux_t {

  samFile *in;
  bam_hdr_t *h;
  hts_idx_t *idx;
  hts_itr_t *itr;
  int min_mapQ;
  char *name;

  snp_t *snpbases;

  int ctot_flag;

} aux_t;

/*! @function
  @abstract Fill aux_t struct array from input bam list.
  @param  blist        Filename of list with bam files and names.
  @param  mapQ         Minimum mapping quality filter.
  @param  data         aux_t struct array to be filled.
  @param  l            Number of positions in map file (for snp_t array init.).
  @param  ctot_flag    Option flag for CtoT filtering (-c option).
  @return              Status indicator.
 */
int read_aux(char *blist, int mapQ, aux_t *data, int l, int ctot_flag);

/*! @function
  @abstract Free the aux_t struct array.
  @param  data     Pointer to aux_t array.
  @param  b        Length of aux_t array.
  @param  l        Number of positions (from map file).
 */
void free_aux(aux_t *data, int b, int l);

/*! @function
  @abstract Read base information from bam files into snp_t array.
  @param  data       Single aux_t struct to be read from.
  @param  l          length of map_arr pos_t array.
  @param  map_arr    pos_t array from map file.
  @return            Status indicator.
 */
int read_bases(aux_t *data, int l, pos_t *map_arr);

/*! @function
  @abstract Sampling of bases and their output.
  @param  data    Single aux_t struct to be sampled from.
  @param  l       Number of positions (length aux_t.snpbases).
  @param  a       Sample one or two alleles.
  @param  nrep    Sample each file multiple (nrep) times.
  @param  out     outf_t struct with outfile handles.
  @param  peda    Struct of PED line array for shared derived test.
  @return         Status indicator.

  @discussion This function also allocates the variables for and, if applicable,
  calls the shared derived test for each bam file. That's why it also requires
  the data stored in peda_t *peda.
 */
int bam_line(aux_t *data, int l, int a, int nrep, outf_t *out, peda_t *peda);

#endif
