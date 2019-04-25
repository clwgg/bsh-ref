#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

#include "map.h"
#include "deriv.h"
#include "util.h"
#include "ped.h"
#include "mats.h"

#include "htslib/htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int ped_line(pedl_t *line, int a, FILE *out)
{

  fprintf(out, "%s %s 0 0 0 -9", line->fid, line->iid);

  mats_m *g = NULL;
  if (a == 1) {
    g = init_m(line->c, 1);
  }

  int i;
  for(i = 0; i < line->c; ++i) {
    if (a == 2) {
      fprintf(out, " %c %c", get_m(line->a, 0, i), get_m(line->a, 1, i));
    }
    else if (a == 1) {
      int r = rand() % 2;
      mats_s c = col_s(line->a, i);
      char s = get_s(c, r);
      fprintf(out, " %c %c", s, s);

      set_m(g, 0, i, s);
    }
  }
  fprintf(out, "\n");

  if (a == 1) {
    free_m(line->a);
    line->a = g;
  }

  return 0;
}

void free_pedl(pedl_t *line)
{
  free_m(line->a);
  free(line->fid);
  free(line->iid);
  free(line);
}

void free_adl(ad_t *adl, int adc)
{
  int j;
  for(j = 0; j < adc; ++j) {
    if(adl[j].fid) free(adl[j].fid);
    if(adl[j].iid) free(adl[j].iid);
  }
  free(adl);
}

int read_ped(char *fn, int l, int a, FILE *out, pos_t *map, char *adfn, peda_t *peda)
{

  gzFile file;
  kstream_t *ks;
  kstring_t str = {0,0,NULL};

  int adc = 0;
  ad_t *adl = NULL;
  if (adfn) {
    adc = nrows(adfn);
    if (adc) {
      adl = malloc(sizeof *adl * adc);
      memset(adl, 0, sizeof *adl * adc);
    }
    if (adl) {
      FILE *adf = NULL;
      adf = fopen(adfn, "r");
      if (adf) {
        int j, ret;
        for (j = 0; j < adc; ++j) {
          ret = fscanf(adf, "%d %ms %ms", &adl[j].flag, &adl[j].fid, &adl[j].iid);
          if (ret != 3) {
            free_adl(adl, adc);
            fclose(adf);
            return -3;
          }
        }
        fclose(adf);
      }
    }
  }

  int fc = 0, ret;

  file = gzopen(fn, "r");
  if (!file) return -1;

  ks = ks_init(file);

  while (ks_getuntil(ks, '\n', &str, 0) >= 0) {

    pedl_t *line = NULL;
    line = malloc(sizeof *line);
    memset(line, 0, sizeof *line);

    line->pos = map;
    line->a = init_m(l, 2);

    char *token;
    token = strtok(str.s, " ");
    line->fid = malloc(strlen(token) + 2);
    strcpy(line->fid, token);

    token = strtok(0, " ");
    line->iid = malloc(strlen(token) + 2);
    strcpy(line->iid, token);

    int i;

    for(i = 0; i < adc; ++i) {
      if(strcmp(line->fid, adl[i].fid) == 0 && strcmp(line->iid, adl[i].iid) == 0) {
        line->flag = adl[i].flag;
        break;
      }
    }

    for(i = 0; i<4; ++i) {
      token = strtok(0, " ");
    }

    i = 0;
    while((token = strtok(0, " "))) {

      set_m(line->a, 0, i, *token);

      token = strtok(0, " ");
      set_m(line->a, 1, i, *token);

      if (strchr("ACGT0", get_m(line->a, 0, i)) == NULL ||
          strchr("ACGT0", get_m(line->a, 1, i)) == NULL ) {
        line->c = l;  // Set count to l to allow freeing of entire line->pos
        free_pedl(line);
        ret = -4;
        goto err;
      }

      // keep information on alleles in reference panel in map[i].s
      mats_s buf = col_s(line->a, i);
      char ref[4] = {'A','C','G','T'};
      int j,k;
      for (k = 0; k < 4; ++k) {
        for (j = 0; j < 2; ++j) {
          if (get_s(buf, j) == ref[k]) map[i].s |= 1 << k;
        }
      }

      ++i;
    }
    line->c = i;

    if (i != l) {
      line->c = l;  // Set count to l to allow freeing of entire line->pos
      free_pedl(line);
      ret = -2;
      goto err;
    }

    ped_line(line, a, out);

    if (line->flag) {
      peda->pedk[fc] = line;
      ++fc;
    } else {
      free_pedl(line);
    }

  }

  ret = fc;

 err:
  free_adl(adl, adc);
  gzclose(file);
  ks_destroy(ks);
  free(str.s);

  return ret;
}

void free_peda(peda_t *peda)
{
  if (!peda->pedk) return;

  int j;
  for(j = 0; j < peda->c; ++j) {
    if(peda->pedk[j]) free_pedl(peda->pedk[j]);
  }
  free(peda->pedk);
}
