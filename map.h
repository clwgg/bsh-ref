#ifndef MAP_H
#define MAP_H

/*! @typedef
  @abstract One informative reference position.
  @field  id    Reference sequence name.
  @field  rp    Position of site in reference sequence.
  @field  s     All alleles present in panel at this position.

  @discussion In load_map(), the id and rp fields are filled for one pos_t array
  of all positions. Additionally, the pos_t array created in load_map() receives
  information in the s field, about what alleles are present at each position in
  the entirety of the PED file. This is later used to copy this information into
  the ped field of the snp_t struct defined in bam.h.
 */
typedef struct pos_t {

  char *id;
  int rp;
  char s;

} pos_t;

/*! @function
  @abstract Load position information from map file.
  @param  fn     File name of map file.
  @param  map    Array of pos_t, with length = number of positions.
  @return        Status indicator -or- number of lines/entries read.
 */
int load_map(char *fn, pos_t *map);

/*! @function
  @abstract Free pos_t array.
  @param  map    Array of pos_t structs.
  @param  l      Length of array in map.
 */
void free_map(pos_t *map, int l);

#endif
