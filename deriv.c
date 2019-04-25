#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "ped.h"
#include "deriv.h"

int count_deriv(peda_t *peda, int i, char s, int *infor, int *share)
{
  int j, k;

  char ref[4] = {'A','C','G','T'};
  int ac = 0;                          // Keep track of test sample
  unsigned long long pbase[4] = {0,0,0,0};
  char anc = '0';
  for(j = 0; j < peda->c; ++j) {
    if (peda->pedk[j]->flag == 1 && anc == '0') {
      /* Set initial ancestral base */
      anc = get_m(peda->pedk[j]->a, 0, i);
      continue;
    } else if (peda->pedk[j]->flag == 1 && anc != '0' && anc != get_m(peda->pedk[j]->a, 0, i)) {
      /* Ancestral base inconsistent */
      return 1;
    } else if (peda->pedk[j]->flag == 2) {
      for(k = 0; k < 4; ++k) {
        if (get_m(peda->pedk[j]->a, 0, i) == ref[k]) {
          pbase[k] |= 1 << ac;         // Set bit for test sample ac+1 in pbase
        }
      }
      ++ac;
    }
  }
  if (anc == '0') return 1;            // Need ancestral base established

  int sbase[4] = {0,0,0,0};
  for(k = 0; k < 4; ++k) {
    /* Keep track of sampled base and ancetral base */
    if (s == ref[k]) sbase[k] |= 1;
    if (anc == ref[k]) sbase[k] |= 2;
  }

  /*
    The main trick below is, that the decimal representation of the bit
    array in pbase is at the same time the index in the infor and share
    arrays. Example: derived base A is shared in test samples 1 and 2
                     => bit array pbase[0] = .... 00000011
                     this is 3 in decimal
                     infor[3] and share[3] correspond to combination 1-2
  */
  for(k = 0; k < 4; ++k) {
    if (pbase[k] && (sbase[k] >> 1) ^ 1) {
      /* Base is set, and different from ancestral (XOR with 1) */
      ++infor[pbase[k]];
      if (sbase[k] & 1) {
        /* Sampled base is the same as derived base in pbase[k] combination. */
        ++share[pbase[k]];
      }
    }
  }
  return 0;
}

int print_deriv_header(int cdc, peda_t *peda, FILE *adout)
/* Print all combinations of FIDs from test samples in the right order. */
{
  int i,j;

  fprintf(adout, "fid\tiid");
  for(i = 1; i < cdc; ++i) {
    /* Loop over combinations */
    int ac = 0;               // Keep track of test sample
    int c = 0;                // Keep track of samples in combination
    fprintf(adout, "\t");
    for(j = 0; j < peda->c; ++j) {
      if (peda->pedk[j]->flag != 2) continue;
      if ((i >> ac) & 1) {    // Is sample ac part of combination i?
        if (c) fprintf(adout, "-");
        fprintf(adout, "%s_%s", peda->pedk[j]->fid, peda->pedk[j]->iid);
        ++c;
      }
      ++ac;
    }
  }
  fprintf(adout, "\n");

  return 0;
}

int print_deriv(int cdc, char *name, int *infor, int *share, FILE *adout)
{
  int i;

  fprintf(adout, "%s", name);
  for (i = 1; i < cdc; ++i) {
    fprintf(adout, "\t%d/%d", share[i], infor[i]);
  }
  fprintf(adout, "\n");

  return 0;
}

int ped_deriv(peda_t *peda, int l, FILE *adout)
{
  int i, j;

  int dc = 0;
  int sc = 0;
  for (i = 0; i < peda->c; ++i) {
    if (peda->pedk[i]->flag == 2) ++dc;
    if (peda->pedk[i]->flag == 3) ++sc;
  }

  if (dc == 0) return 1;

  int cdc = pow(2, dc);

  if (dc > CHAR_BIT * sizeof(unsigned long long)) {
    /* Too many test populations to keep in ULL bit array */
    return -1;
  }

  print_deriv_header(cdc, peda, adout);

  if (sc == 0) return 2;

  int **infor = NULL;
  int **share = NULL;
  infor = malloc(sizeof *infor * sc);
  share = malloc(sizeof *share * sc);
  for (i = 0; i < sc; ++i) {
    infor[i] = calloc(cdc, sizeof *infor[i]);
    share[i] = calloc(cdc, sizeof *share[i]);
  }

  for (i = 0; i < l; ++i) {
    int csc = 0;
    for (j = 0; j < peda->c; ++j) {
      if (peda->pedk[j]->flag == 3) {
        if (get_m(peda->pedk[j]->a, 0, i) == '0') continue;
        count_deriv(peda, i, get_m(peda->pedk[j]->a, 0, i), infor[csc], share[csc]);
        ++csc;
      }
    }
  }

  int csc = 0;
  for (j = 0; j < peda->c; ++j) {
    if (peda->pedk[j]->flag == 3) {
      char *name = NULL;
      name = calloc(strlen(peda->pedk[j]->fid) + strlen(peda->pedk[j]->iid) + 4, sizeof *name);
      memcpy(name, peda->pedk[j]->fid, strlen(peda->pedk[j]->fid));
      strcat(strcat(name, "\t"), peda->pedk[j]->iid);
      print_deriv(cdc, name, infor[csc], share[csc], adout);
      ++csc;
      free(name);
    }
  }

  for (i = 0; i < sc; ++i) {
    free(infor[i]); free(share[i]);
  }
  free(infor); free(share);
  return 0;
}
