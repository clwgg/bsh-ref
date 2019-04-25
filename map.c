#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

#include "map.h"

#include "htslib/htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int load_map(char *fn, pos_t *map)
{

  gzFile file;
  kstream_t *ks;
  kstring_t str = {0,0,NULL};

  file = gzopen(fn, "r");
  if (!file) return 0;

  ks = ks_init(file);

  int i = 0;

  while (ks_getuntil(ks, '\n', &str, 0) >= 0) { // go line by line

    char *id, *pos, *token;                     // keep first(id) and last(pos) element
    id = strtok(str.s, " \t");
    while((token = strtok(NULL, " \t"))) {
      pos = token;
    }
    int rp = atoi(pos);                         // the following assumes 1-based coordinates
    map[i].rp = rp;
    map[i].id = malloc(sizeof *id * str.m);
    memset(map[i].id, 0, sizeof *id * str.m);
    memcpy(map[i].id, id, sizeof *id * str.m);

    ++i;

  }

  gzclose(file);
  ks_destroy(ks);
  free(str.s);

  return i;
}

void free_map(pos_t *map, int l)
{

  if (!map) return;

  int i;
  for(i = 0; i < l; ++i) {
    if(map[i].id) free(map[i].id);
  }

  free(map);
}
