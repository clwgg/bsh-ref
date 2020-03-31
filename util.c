#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <limits.h>

#include "util.h"

#include "htslib/htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int nrows(char *fn)
{

  gzFile file;
  kstream_t *ks;
  kstring_t str = {0,0,NULL};

  file = gzopen(fn, "r");
  if (!file) return -1;

  ks = ks_init(file);

  int i = 0;

  while (ks_getuntil(ks, '\n', &str, 0) >= 0) {
    ++i;
  }

  gzclose(file);
  ks_destroy(ks);
  free(str.s);

  return i;
}

int ncols(char *fn, int offset)
{

  gzFile file;
  kstream_t *ks;
  kstring_t str = {0,0,NULL};

  file = gzopen(fn, "r");
  if (!file) return -1;

  ks = ks_init(file);

  ks_getuntil(ks, '\n', &str, 0);

  char *token;
  token = strtok(str.s, " \t");
  int i = 1;
  while((token = strtok(NULL, " \t"))) {
    ++i;
  }

  gzclose(file);
  ks_destroy(ks);
  free(str.s);

  return i - offset;
}


void close_outf(outf_t *outf)
{
  if (outf->pfout) fclose(outf->pfout);
  if (outf->adfout) fclose(outf->adfout);
  if (outf->safout) fclose(outf->safout);
}

int numPlaces(int n)
{
  if (n < 0) return numPlaces ((n == INT_MIN) ? INT_MAX: -n);
  if (n < 10) return 1;
  return 1 + numPlaces (n / 10);
}
