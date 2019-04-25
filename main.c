#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "bam.h"
#include "map.h"
#include "ped.h"
#include "util.h"
#include "deriv.h"

static int usage(char **argv)
{
  printf("\nUsage: %s [options] -m in.map [-p in.ped | -b in.lst]\n\n", argv[0]);
  printf("Options:\n");
  printf("\t-m\tmap file containing positions (required)\n");
  printf("\t-p\tped file containing known individuals (at least one of -p and -b)\n");
  printf("\t-b\tbam list file containing new individuals (at least one of -p and -b)\n");
  printf("\t-o\toutfile basename (default: out)\n\n");

  printf("\t-a\tsample 1 allele (pseudo-haploid) or 2 (default: 1)\n");
  printf("\t-r\tsample each bam file -r times (default: 1)\n");
  printf("\t-q\tmin map quality (default: 0)\n");
  printf("\t-s\tcreate stat file for bam sampling\n");
  printf("\t-c\trestrict base sampling to reads with terminal C-to-T in one (1) or both (2) ends (default: 0)\n");
  printf("\t  \t(Warning: this is so far only implemented for single stranded libraries, which also have C-to-T at 3p)\n\n");

  printf("\t-t\tancestral vs. derived file (needs -p)\n\n");

  return 1;
}

void init_random(void)
{
  char *env = getenv("RANDOM_SEED");
  if (env == NULL || strcmp(env, "") == 0) {
    srand(time(NULL));
  } else {
    srand(atoi(env));
  }
}

int main(int argc, char **argv)
{

  init_random();

  int elem;
  int mapQ = 0;
  char *mapf = NULL;
  char *pedf = NULL;
  char *out = "out";
  int a = 1;
  int repb = 0;

  char *adfile = NULL;
  char *blist = NULL;

  int sa_flag = 0;
  int ctot_flag = 0;


  while (( elem = getopt(argc, argv, "b:q:m:p:a:r:o:t:sc:") ) >= 0) {
    switch(elem) {
      case 'q': mapQ = atoi(optarg); break;
      case 'm': mapf = optarg; break;
      case 'b': blist = optarg; break;
      case 'p': pedf = optarg; break;
      case 'o': out = optarg; break;
      case 't': adfile = optarg; break;
      case 'a': a = atoi(optarg); break;
      case 'r': repb = atoi(optarg); break;
      case 's': sa_flag = 1; break;
      case 'c': ctot_flag = atoi(optarg); break;
    }
  }

  if (!mapf) {
    fprintf(stderr, "Map file is required.\n\n");
    return usage(argv);
  }
  if (!(blist || pedf)) {
    fprintf(stderr, "Either ped file or bam list file is required.\n\n");
    return usage(argv);
  }
  if (a != 1 && a != 2) {
    fprintf(stderr, "-a has to be either 1 or 2.\n\n");
    return usage(argv);
  }
  if (a != 1 && adfile) {
    fprintf(stderr, "ancestral vs. derived so far only works with -a = 1\n\n");
    return usage(argv);
  }
  if (adfile && !pedf) {
    fprintf(stderr, "ped file input is needed for ancestral vs. derived\n\n");
    return usage(argv);
  }
  if (ctot_flag > 2 || ctot_flag < 0) {
    fprintf(stderr, "C-to-T flag (-c) should be 0,1 or 2\n\n");
    return usage(argv);
  }


  int ret = 0;
  int adc = 0;
  int b = 0;
  int l = 0;

  peda_t peda = {NULL, 0, 0};
  pos_t *map_arr = NULL;
  aux_t *data = NULL;


  outf_t outf = {NULL, NULL, NULL};

  char *outp = malloc( sizeof(char) * (strlen(out) + 10) );
  sprintf(outp, "%s.ped", out);
  outf.pfout = fopen(outp, "w+");

  if (!outf.pfout) {
    fprintf(stderr, "Could not open output ped file\n\n");
    goto err;
  }


  l = nrows(mapf);
  if (l == -1) {
    fprintf(stderr, "Could not open map file.\n\n");
    goto err;
  }
  if (ncols(mapf, 0) != 4) {
    fprintf(stderr, "map file expects 4 columns - first: sequence name; last: base coordinate.\n\n");
    goto err;
  }

  map_arr = malloc(sizeof *map_arr * l);
  memset(map_arr, 0, sizeof *map_arr * l);

  ret = load_map(mapf, map_arr);
  if (ret == 0) {
    fprintf(stderr, "Could not open/read map file.\n\n");
    goto err;
  }

  if (adfile) {
    adc = nrows(adfile);
    if (adc == -1) {
      fprintf(stderr, "Could not open file from -t option.\n\n");
      goto err;
    }
    if (adc) {
      peda.pedk = calloc(adc, sizeof *peda.pedk);
      peda.c = adc;
      peda.flag = 1;

      sprintf(outp, "%s.adstat", out);
      outf.adfout = fopen(outp, "w+");

      if (!outf.adfout) {
        fprintf(stderr, "Could not open output adstats file\n\n");
        goto err;
      }
    }
  }

  if (pedf) {
    int s = ncols(pedf, 6);
    if (s == -1) {
      fprintf(stderr, "Could not open ped file.\n\n");
      goto err;
    }
    if (s/l != 2) {
      fprintf(stderr, "ped file not properly biallelic.\n\n");
      goto err;
    }

    ret = read_ped(pedf, l, a, outf.pfout, map_arr, adfile, &peda);
    if (ret >= 0 && ret != adc) {
      fprintf(stderr, "Inconsistent parsing of file from -t option.\n\n");
      goto err;
    }
    if (ret == -3) {
      fprintf(stderr, "File from -t option misformed.\n\n");
      goto err;
    }
    if (ret == -1) {
      fprintf(stderr, "Could not open ped file.\n\n");
      goto err;
    }
    if (ret == -2) {
      fprintf(stderr, "Column number inconsistencies in ped file\n\n");
      goto err;
    }
    if (ret == -4) {
      fprintf(stderr, "ped file not of 'ACGT0' format.\n\n");
      goto err;
    }


    if (adc) {
      ret = ped_deriv(&peda, l, outf.adfout);
      if (ret == -1) {
        fprintf(stderr, "Too many test populations in shared derived test.\n\n");
        goto err;
      }
    }

  }

  if (blist) {
    b = nrows(blist);
    if (b == -1) {
      fprintf(stderr, "Could not open list file.\n\n");
      goto err;
    }

    data = malloc(sizeof *data * b);
    memset(data, 0, sizeof *data * b);

    ret = read_aux(blist, mapQ, data, l, ctot_flag);
    if (ret == 0) {
      fprintf(stderr, "Could not open/read list file.\n\n");
      goto err;
    }
    else if (ret == -2) {
      fprintf(stderr, "Could not open bam file.\n\n");
      goto err;
    }
    else if (ret == -3) {
      fprintf(stderr, "Bam file needs to be indexed.\n\n");
      goto err;
    }

    if (sa_flag) {
      sprintf(outp, "%s.sastat", out);
      outf.safout = fopen(outp, "w+");

      if (!outf.safout) {
        fprintf(stderr, "Could not open output sample stats file\n\n");
        goto err;
      }
    }

    int i;
    for(i = 0; i < b; ++i) {
      if (!data[i].in) continue;

      ret = read_bases(&data[i], l, map_arr);

      int j = 0;
      if (repb) j = 1;
      while (j <= repb) {
        ret = bam_line(&data[i], l, a, j, &outf, &peda);
        j++;
        if (ret == -1) {
          fprintf(stderr, "Too many test populations in shared derived test.\n\n");
          goto err;
        }
      }
    }
  }

  free(outp);
  close_outf(&outf);
  free_aux(data, b, l);
  free_map(map_arr, l);
  free_peda(&peda);

  return 0;

 err:
  if (outp) free(outp);
  close_outf(&outf);
  free_aux(data, b, l);
  free_map(map_arr, l);
  free_peda(&peda);
  return usage(argv);

}
