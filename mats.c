#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mats.h"

mats_m* init_m(int ncol, int nrow)
{
  mats_m *m = NULL;

  m = malloc(sizeof *m);

  m->m = calloc(nrow * ncol, sizeof m->m[0]);

  m->c = calloc(nrow, sizeof m->c[0]);
  m->r = calloc(ncol, sizeof m->r[0]);

  m->nrow = nrow;
  m->ncol = ncol;

  return m;
}

void free_m(mats_m *m)
{
  free(m->m);
  free(m->r);
  free(m->c);
  free(m);
}

void set_m(mats_m *m, int r, int c, int e)
{
  m->m[r * m->ncol + c] = e;
}

int get_m(mats_m *m, int r, int c)
{
  return m->m[r * m->ncol + c];
}

////////////////////

mats_v* init_v(int n)
{
  mats_v *v = NULL;

  v = malloc(sizeof *v);

  v->v = calloc(n, sizeof v->v[0]);

  v->n = n;

  return v;
}

void free_v(mats_v *v)
{
  free(v->v);
  free(v);
}

void set_v(mats_v *v, int n, int e)
{
  v->v[n] = e;
}

int get_v(mats_v *v, int n)
{
  return v->v[n];
}

////////////////////

mats_s col_s(mats_m *m, int c)
{
  int i;
  for (i = 0; i < m->nrow; ++i) {
    m->c[i] = &m->m[i * m->ncol + c];
  }

  mats_s s;
  s.v = m->c;
  s.n = m->nrow;

  return s;
}

mats_s row_s(mats_m *m, int r)
{
  int i;
  for (i = 0; i < m->ncol; ++i) {
    m->r[i] = &m->m[r * m->ncol + i];
  }

  mats_s s;
  s.v = m->r;
  s.n = m->ncol;

  return s;
}

void set_s(mats_s s, int n, int e)
{
  *s.v[n] = e;
}

int get_s(mats_s s, int n)
{
  return *s.v[n];
}

mats_v* s_write_v(mats_s s)
{
  mats_v *v = NULL;
  int n = s.n;

  v = malloc(sizeof *v);

  v->v = calloc(n, sizeof v->v[0]);

  int i;
  for (i = 0; i < n; ++i) {
    v->v[i] = *s.v[i];
  }

  v->n = n;

  return v;
}

////////////////////

void zeroc_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      set_m(m, i, j, '0');
    }
  }
}

void zeroi_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      set_m(m, i, j, 0);
    }
  }
}

void zeroc_v(mats_v *v)
{
  int i;
  for (i = 0; i < v->n; ++i) {
    set_v(v, i, '0');
  }
}

void zeroi_v(mats_v *v)
{
  int i;
  for (i = 0; i < v->n; ++i) {
    set_v(v, i, 0);
  }
}

void randomc_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      set_m(m, i, j, rand() % 26 + 'A');
    }
  }
}

void randomi_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      set_m(m, i, j, rand() % 100 + 100);
    }
  }
}

void printc_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      printf("  %c", get_m(m, i, j));
    }
    printf("\n");
  }
  printf("\n");
}

void printc_v(mats_v *v)
{
  int i;
  for (i = 0; i < v->n; ++i) {
    printf("  %c", get_v(v, i));
  }
  printf("\n\n");
}

void printi_m(mats_m *m)
{
  int i, j;
  for (i = 0; i < m->nrow; ++i) {
    for (j = 0; j < m->ncol; ++j) {
      printf("  %d", get_m(m, i, j));
    }
    printf("\n");
  }
  printf("\n");
}

void printi_v(mats_v *v)
{
  int i;
  for (i = 0; i < v->n; ++i) {
    printf("  %d", get_v(v, i));
  }
  printf("\n\n");
}

////////////////////

int rbind_m(mats_m *sink, mats_m *src)
{

  if (sink->ncol != src->ncol) {
    return -1;
  }

  sink->m = realloc(sink->m, sizeof sink->m[0] * sink->nrow * sink->ncol +
                             sizeof  src->m[0] *  src->nrow *  src->ncol);

  memcpy(sink->m + (sink->nrow * sink->ncol),
         src->m, sizeof  src->m[0] *  src->nrow *  src->ncol);

  sink->c = realloc(sink->c, sizeof sink->c[0] * sink->nrow +
                             sizeof  src->c[0] *  src->nrow);

  sink->nrow = sink->nrow + src->nrow;

  return 0;
}

////////////////////

int sum_v(mats_v *v)
{
  int i, sum = 0;
  for (i = 0; i < v->n; ++i) {
    sum += get_v(v, i);
  }
  return sum;
}

////////////////////
