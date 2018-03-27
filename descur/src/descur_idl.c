#include <stdio.h>

void curve3d_c(argc,argv)
     int argc;
     void *argv[];

{
  void curve3d_();

  int *i0, *i1, *i2, *i3, *i4, *i5, *i11;
  double *f6, *f7, *f8, *f9, *f10, *f12, *f13, *f14, *f15, *f16, *f17; 

  i0 = (int *) argv[0];
  i1 = (int *) argv[1];
  i2 = (int *) argv[2];
  i3 = (int *) argv[3];
  i4 = (int *) argv[4];
  i5 = (int *) argv[5];
  f6 = (double *) argv[6];
  f7 = (double *) argv[7];
  f8 = (double *) argv[8];
  f9 = (double *) argv[9];
  f10 = (double *) argv[10];
  i11 = (int *) argv[11];
  f12 = (double *) argv[12];
  f13 = (double *) argv[13];
  f14 = (double *) argv[14];
  f15 = (double *) argv[15];
  f16 = (double *) argv[16];
  f17 = (double *) argv[17];


  curve3d_(i0, i1, i2, i3, i4, i5, f6, f7, f8, f9, f10, i11, f12, f13, f14, f15, f16, f17);
}
