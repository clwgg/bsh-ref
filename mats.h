#ifndef MATS_H
#define MATS_H

typedef struct {
  int *m;
  int nrow;
  int ncol;
  int **r;
  int **c;
} mats_m;

typedef struct {
  int *v;
  int n;
} mats_v;

typedef struct {
  int **v;
  int n;
} mats_s;


mats_m* init_m(int ncol, int nrow);
void free_m(mats_m *m);
void set_m(mats_m *m, int r, int c, int e);
int get_m(mats_m *m, int r, int c);


mats_v* init_v(int n);
void free_v(mats_v *v);
void set_v(mats_v *v, int n, int e);
int get_v(mats_v *v, int n);


mats_s col_s(mats_m *m, int c);
mats_s row_s(mats_m *m, int r);
void set_s(mats_s s, int n, int e);
int get_s(mats_s s, int n);
mats_v* s_write_v(mats_s s);


void zeroc_m(mats_m *m);
void zeroi_m(mats_m *m);
void zeroc_v(mats_v *v);
void zeroi_v(mats_v *v);
void randomc_m(mats_m *m);
void randomi_m(mats_m *m);
void printc_m(mats_m *m);
void printc_v(mats_v *v);
void printi_m(mats_m *m);
void printi_v(mats_v *v);


int rbind_m(mats_m *sink, mats_m *src);


int sum_v(mats_v *v);

#endif
