#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <math.h>
#include <limits.h>

#include "map.h"
#include "ped.h"
#include "util.h"
#include "deriv.h"
#include "bam.h"

#include "htslib/htslib/sam.h"

#include "htslib/htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int read_aux(char *blist, int mapQ, aux_t *data, int l, int ctot_flag)
{

  gzFile file;
  kstream_t *ks;
  kstring_t str = {0,0,NULL};

  file = gzopen(blist, "r");
  if (!file) return 0;

  ks = ks_init(file);

  int i = 0;

  while (ks_getuntil(ks, '\n', &str, 0) >= 0) { // go line by line

    if (*str.s == '#') {
      ++i; continue;
    }

    char *name, *bfile;
    name = strtok(str.s, "\t");
    bfile = strtok(0, "\t");

    data[i].name = malloc(sizeof *name * str.m);
    memset(data[i].name, 0, sizeof *name * str.m);
    memcpy(data[i].name, name, sizeof *name * str.m);

    data[i].min_mapQ = mapQ;

    data[i].in = sam_open(bfile, "r");
    if (!data[i].in) {
      gzclose(file);
      ks_destroy(ks);
      free(str.s);

      return -2;
    }

    data[i].h = sam_hdr_read(data[i].in);

    data[i].idx = sam_index_load(data[i].in, bfile);
    if (!data[i].idx) {
      gzclose(file);
      ks_destroy(ks);
      free(str.s);

      return -3;
    }

    data[i].itr = NULL;

    // initialize snpbases
    data[i].snpbases = malloc(sizeof *data[i].snpbases * l);
    memset(data[i].snpbases, 0, sizeof *data[i].snpbases * l);

    // set flags
    data[i].ctot_flag = ctot_flag;

    ++i;
  }

  gzclose(file);
  ks_destroy(ks);
  free(str.s);

  return i;
}

void free_aux(aux_t *data, int b, int l)
{

  if (!data) return;

  int i, j;
  for(i = 0; i < b; ++i) {
    if(data[i].h) bam_hdr_destroy(data[i].h);
    if(data[i].idx) hts_idx_destroy(data[i].idx);
    if(data[i].itr) hts_itr_destroy(data[i].itr);
    if(data[i].in) sam_close(data[i].in);
    if(data[i].name) free(data[i].name);

    if(data[i].snpbases) {
      for(j = 0; j < l; ++j) {
        if(data[i].snpbases[j].bases) free(data[i].snpbases[j].bases);
      }
      free(data[i].snpbases);
    }
  }

  free(data);
}

int read_bam(void *data, bam1_t *b)
{

  aux_t *aux = (aux_t*)data;
  int ret;

  while (1) {
    ret = aux->itr ? sam_itr_next(aux->in, aux->itr, b) : sam_read1(aux->in, aux->h, b);
    if (ret < 0) break;
    if ( (int)b->core.qual < aux->min_mapQ ) continue;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FDUP) ) continue;
    break;
  }
  return ret;

}

char t_cig(int i)
{
  char c = 0;
  switch (i) {
  case 0: c = 'M'; break;
  case 1: c = 'I'; break;
  case 2: c = 'D'; break;
  case 3: c = 'N'; break;
  case 4: c = 'S'; break;
  case 5: c = 'H'; break;
  case 6: c = 'P'; break;
  case 7: c = '='; break;
  case 8: c = 'X'; break;
  }
  return c;
}

char t_seq(int i)
{
  char c = 0;
  switch (i) {
  case 1:  c = 'A'; break;
  case 2:  c = 'C'; break;
  case 4:  c = 'G'; break;
  case 8:  c = 'T'; break;
  case 15: c = 'N'; break;
  }
  return c;
}

char t_comp(char c)
{
  char out = 0;
  switch (c) {
  case 'A': out = 'T'; break;
  case 'C': out = 'G'; break;
  case 'G': out = 'C'; break;
  case 'T': out = 'A'; break;
  case 'a': out = 't'; break;
  case 'c': out = 'g'; break;
  case 'g': out = 'c'; break;
  case 't': out = 'a'; break;
  }
  return out;
}

uint8_t term_ctot(bam1_t *b)
{

  uint8_t out = 0;

  uint32_t *cig = bam_get_cigar(b);
  uint8_t *pMD = bam_aux_get(b, "MD");
  uint8_t *seq = bam_get_seq(b);
  char *md = bam_aux2Z(pMD);

  char op;
  char base, ref;
  // Check beginning of read (may be original 5p or 3p)
  op = t_cig(bam_cigar_op(cig[0]));
  if (strchr("M=X", op)) {
    /* Only look if CIGAR shows (mis-)match of first base */
    if (md[0] == '0') {
      /* Only look, if MD shows mismatch at first base */
      base = t_seq(bam_seqi(seq, 0));
      ref = md[1];
      if (bam_is_rev(b)) {
        /* If the read is reversed, take complement
           of bases and set counter for 3p end
           if the reversed bases show C->T */
        base = t_comp(base);
        ref = t_comp(ref);
        if (strchr("Cc", ref) && strchr("Tt", base)) {
          out |= 2;
        }
      }
      else if (strchr("Cc", ref) && strchr("Tt", base)) {
        /* If the read isn't reversed and shows C->T, set
           counter for 5p end */
        out |= 1;
      }
    }
  }
  // Check end of read (may be original 3p or 5p)
  op = t_cig(bam_cigar_op(cig[b->core.n_cigar - 1]));
  if (strchr("M=X", op)) {
    /* Conditionals as above, only with last base */
    if (md[strlen(md)-1] == '0' && isalpha(md[strlen(md)-2])) {
      base = t_seq(bam_seqi(seq, b->core.l_qseq - 1));
      ref = md[strlen(md)-2];
      if (bam_is_rev(b)) {
        /* If the read is reversed, take complement
           of bases and set counter for 5p end
           if the reversed bases show C->T */
        base = t_comp(base);
        ref = t_comp(ref);
        if (strchr("Cc", ref) && strchr("Tt", base)) {
          out |= 1;
        }
      }
      else if (strchr("Cc", ref) && strchr("Tt", base)) {
        /* If the read isn't reversed and shows C->T, set
           counter for 3p end */
        out |= 2;
      }
    }
  }

  return out;
}

int read_pos(int n, const bam_pileup1_t *plp, snp_t *snp)
{
  int i,j = 0;

  snp->bases = malloc(sizeof *snp->bases * n);
  memset(snp->bases, 0, sizeof *snp->bases * n);

  for (i = 0; i < n; ++i) {

    const bam_pileup1_t *p = plp + i;
    uint8_t *seq = bam_get_seq(p->b);
    uint8_t nuc = bam_seqi(seq, p->qpos);

    if (p->is_del || nuc == 15) continue; // skip deletions and Ns

    snp->bases[j].base = nuc;

    snp->bases[j].ctot = term_ctot(p->b);

    if (p->b->core.l_qseq <= UINT16_MAX) {
      snp->bases[j].pos = p->qpos;
      snp->bases[j].qlen = p->b->core.l_qseq;
    } else {
      fprintf(stderr, "Reads longer than %d are not supported for positional data.\n", UINT16_MAX);
    }

    ++j;
  }

  snp->c = j;

  return j;

}

int read_bases(aux_t *data, int l, pos_t *map_arr)
{
  bam_plp_t iter;
  const bam_pileup1_t *plp;
  int tid, p, n;

  int i;
  for(i = 0; i < l; ++i) {

    data->snpbases[i].ped = map_arr[i].s; // add ped panel information gathered in read_ped() to snp_t

    if (data->itr) {
      hts_itr_destroy(data->itr);
      data->itr = NULL;
    }
    data->itr = sam_itr_queryi(data->idx, bam_name2id(data->h, map_arr[i].id), map_arr[i].rp - 1, map_arr[i].rp);

    iter = bam_plp_init(read_bam, (void*)data);  // pileup iterator (not to be confused with data.itr)

    int c = 0;
    while ((plp = bam_plp_auto(iter, &tid, &p, &n)) != 0) {
      // This loop will return pileup for all reads overlapping data.itr
      ++c;                            // Keep track if any reads overlap data.itr
      if (p+1 != map_arr[i].rp) continue;   // Keep only pileup on top of SNP, take care with 0 vs 1-base

      read_pos(n, plp, &data->snpbases[i]);

    }
    bam_plp_destroy(iter);
  }

  return 0;
}

char sample_pos(aux_t *data, int i, FILE *saout, int nrep)
{
  snp_t *snp = &data->snpbases[i];

  base_t final = {0,0,0,0}; // The sampled base_t is kept here

  if (data->ctot_flag == 0) {
    /* Without CtoT filtering, sample a random base from snp->bases */
    int r = rand() % snp->c;
    final = snp->bases[r];
  }
  else if (data->ctot_flag > 0) {
    /* With CtoT filtering, create a new base_t array *arr to
       keep 'eligible' bases in, depending on the filtering (1/2) */
    base_t *arr = calloc(snp->c, sizeof *arr);
    int n = 0, i;
    for (i = 0; i < snp->c; ++i) {
      if (data->ctot_flag == 1) {
        if (snp->bases[i].ctot > 0) {
          arr[n] = snp->bases[i];
          ++n;
        }
      }
      else if (data->ctot_flag == 2) {
        if (snp->bases[i].ctot == 3) {
          arr[n] = snp->bases[i];
          ++n;
        }
      }
    }

    if (n == 0) {
      /* No eligible bases post filtering */
      final.base = 0;
    } else {
      /* Sample base from arr of filtered bases */
      int r = rand() % n;
      final = arr[r];
    }

    free(arr);
  }

  char c = '0';

  if (final.base == 0)
    return c;
  else
    c = t_seq(final.base);

  if (saout) {
    if (nrep)
      fprintf(saout, "%s\t%d", data->name, nrep);
    else
      fprintf(saout, "%s\t%s", data->name, data->name);

    fprintf(saout, "\t%d\t%d\t%c", i+1, snp->c, c);
    if (final.pos / (double)final.qlen <= 0.5) {
      fprintf(saout, "\t%d", final.pos + 1);
    } else {
      fprintf(saout, "\t%d", final.pos + 1 - final.qlen - 1);
    }

    fprintf(saout, "\t%d", final.ctot);

    char ref[4] = {'A','C','G','T'};
    int k;
    for (k = 0; k < 4; ++k) {
      if (snp->ped & 1 << k) {
        fprintf(saout, "\t%c", ref[k]);
      }
    }
    fprintf(saout, "\n");
  }

  return c;
}

int bam_line(aux_t *data, int l, int a, int nrep, outf_t *out, peda_t *peda)
{

  FILE *pout = out->pfout;
  int i;

  int dc = 0;
  for (i = 0; i < peda->c; ++i) {
    // Count test samples
    if (peda->pedk[i]->flag == 2) ++dc;
  }
  int cdc = pow(2, dc); // Number of possible combinations of test samples

  if (dc > CHAR_BIT * sizeof(unsigned long long)) {
    /* Too many test populations to keep in ULL bit array */
    return -1;
  }

  int *infor = NULL;
  int *share = NULL;
  infor = calloc(cdc, sizeof *infor);
  share = calloc(cdc, sizeof *share);

  if (nrep)
    fprintf(pout, "%s %d 0 0 0 -9", data->name, nrep);
  else
    fprintf(pout, "%s %s 0 0 0 -9", data->name, data->name);

  for(i = 0; i < l; ++i) {

    char s, s1, s2;
    if (data->snpbases[i].c == 0) {
      fprintf(pout, " 0 0");
      continue;
    } else if (a == 1) {
      s = sample_pos(data, i, out->safout, nrep);
      fprintf(pout, " %c %c", s, s);
    } else if (a == 2) {
      s1 = sample_pos(data, i, out->safout, nrep);
      s2 = sample_pos(data, i, out->safout, nrep);
      fprintf(pout, " %c %c", s1, s2);
    }

    if (dc && a == 1 && s != '0') {
      count_deriv(peda, i, s, infor, share);
    }

  }
  fprintf(pout, "\n");

  if (dc) {
    char *name = NULL;
    name = calloc(strlen(data->name) + numPlaces(nrep) + 4, sizeof *name);
    memcpy(name, data->name, strlen(data->name));
    strcat(name, "\t");
    sprintf(name + strlen(name), "%d", nrep);
    print_deriv(cdc, name, infor, share, out->adfout);
    free(name);
  }

  free(share); free(infor);

  return 0;
}
