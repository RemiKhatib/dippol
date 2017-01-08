/*
  Basic mathematical operations
*/

#include <stdio.h>
#include <math.h>
#include <main.h>
#include <my_math.h>

/* Cross product */
void cross_p(double *v1, double *v2, double *v3){
  v3[0]=v1[1]*v2[2]-v2[1]*v1[2];
  v3[1]=v1[2]*v2[0]-v2[2]*v1[0];
  v3[2]=v1[0]*v2[1]-v2[0]*v1[1];

  return;
}



/* Vector normalization */
void norm(double *v1, double *v2){
  double d=0;

  d=sqrt(pow(v1[0],2.)+pow(v1[1],2.)+pow(v1[2],2.));
  v2[0]=v1[0]/d;
  v2[1]=v1[1]/d;
  v2[2]=v1[2]/d;
  
  return;
}

/* Rewrite the molecule coordinate to have the mass center
   between -param/2 and param/2 */
void pbc_mol(vect_3d cell, mol_info *mol){
  /*Calculation of the mass center*/
  (*mol).x= (MO* (*mol).xO + MH*(*mol).xH1 + MH* (*mol).xH2 )/(MO+2*MH);
  (*mol).y= (MO* (*mol).yO + MH*(*mol).yH1 + MH* (*mol).yH2 )/(MO+2*MH);
  (*mol).z= (MO* (*mol).zO + MH*(*mol).zH1 + MH* (*mol).zH2 )/(MO+2*MH);

  pbc(cell,&((*mol).x),&((*mol).y),&((*mol).z));
}

/* Rewrite some coordinates respecting the pbc conditions
   The coordinates will be between -param/2 and param/2 */
void pbc(vect_3d cell, double *x, double *y, double *z){
  *x = *x - rint(*x / cell.x) * cell.x;
  *y = *y - rint(*y / cell.y) * cell.y;
  *z = *z - rint(*z / cell.z) * cell.z;
}


/*Multiplication of a symmetric 3D matrix and a 3D vector*/
void mult_msym_v_3d(mat_sym_3d *m, vect_3d *vin, vect_3d *vout){
  (*vout).x = (*m).xx*(*vin).x + (*m).yx*(*vin).y + (*m).zx*(*vin).z;
  (*vout).y = (*m).yx*(*vin).x + (*m).yy*(*vin).y + (*m).zy*(*vin).z;
  (*vout).z = (*m).zx*(*vin).x + (*m).zy*(*vin).y + (*m).zz*(*vin).z;
}


/*Multiplication of a symmetric 3D matrix and a 3D matrix*/
void mult_msym_m_3d(mat_sym_3d *m1, mat_3d *m2, mat_3d *mout){
  (*mout).xx = (*m1).xx*(*m2).xx + (*m1).yx*(*m2).yx + (*m1).zx*(*m2).zx;
  (*mout).xy = (*m1).xx*(*m2).xy + (*m1).yx*(*m2).yy + (*m1).zx*(*m2).zy;
  (*mout).xz = (*m1).xx*(*m2).xz + (*m1).yx*(*m2).yz + (*m1).zx*(*m2).zz;

  (*mout).yx = (*m1).yx*(*m2).xx + (*m1).yy*(*m2).yx + (*m1).zy*(*m2).zx;
  (*mout).yy = (*m1).yx*(*m2).xy + (*m1).yy*(*m2).yy + (*m1).zy*(*m2).zy;
  (*mout).yz = (*m1).yx*(*m2).xz + (*m1).yy*(*m2).yz + (*m1).zy*(*m2).zz;

  (*mout).zx = (*m1).zx*(*m2).xx + (*m1).zy*(*m2).yx + (*m1).zz*(*m2).zx;
  (*mout).zy = (*m1).zx*(*m2).xy + (*m1).zy*(*m2).yy + (*m1).zz*(*m2).zy;
  (*mout).zz = (*m1).zx*(*m2).xz + (*m1).zy*(*m2).yz + (*m1).zz*(*m2).zz;
}

/*Multiplication of a matrix 3D and another matrix 3D*/
void mult_m_m_3d(mat_3d *m1, mat_3d *m2, mat_3d *mout){

  (*mout).xx = (*m1).xx*(*m2).xx + (*m1).xy*(*m2).yx + (*m1).xz*(*m2).zx;
  (*mout).xy = (*m1).xx*(*m2).xy + (*m1).xy*(*m2).yy + (*m1).xz*(*m2).zy;
  (*mout).xz = (*m1).xx*(*m2).xz + (*m1).xy*(*m2).yz + (*m1).xz*(*m2).zz;

  (*mout).yx = (*m1).yx*(*m2).xx + (*m1).yy*(*m2).yx + (*m1).yz*(*m2).zx;
  (*mout).yy = (*m1).yx*(*m2).xy + (*m1).yy*(*m2).yy + (*m1).yz*(*m2).zy;
  (*mout).yz = (*m1).yx*(*m2).xz + (*m1).yy*(*m2).yz + (*m1).yz*(*m2).zz;

  (*mout).zx = (*m1).zx*(*m2).xx + (*m1).zy*(*m2).yx + (*m1).zz*(*m2).zx;
  (*mout).zy = (*m1).zx*(*m2).xy + (*m1).zy*(*m2).yy + (*m1).zz*(*m2).zy;
  (*mout).zz = (*m1).zx*(*m2).xz + (*m1).zy*(*m2).yz + (*m1).zz*(*m2).zz;

}


/*This subroutine should be removed once the 
structures with the positions will be removed
It just serves to rotate the H atomic polarizability 
from bond referential to lab referential*/
/*AH is symmetric. xz = xy =0*/
void rot_pol0_H(sys_info *sys, mat_sym_3d *pol0, double *v1, double *v2, double *v3){
  double tmp[9];
  
  /* D*AH */
  tmp[0]= v1[0]*(*sys).AH_xx + v2[0]*(*sys).AH_yx + v3[0]*(*sys).AH_zx ;/*Xx*/
  tmp[1]= v1[0]*(*sys).AH_yx + v2[0]*(*sys).AH_yy + v3[0]*(*sys).AH_zy ;/*Xy*/
  tmp[2]= v1[0]*(*sys).AH_zx + v2[0]*(*sys).AH_zy + v3[0]*(*sys).AH_zz ;/*Xz*/
  tmp[3]= v1[1]*(*sys).AH_xx + v2[1]*(*sys).AH_yx + v3[1]*(*sys).AH_zx ;/*Yx*/
  tmp[4]= v1[1]*(*sys).AH_yx + v2[1]*(*sys).AH_yy + v3[1]*(*sys).AH_zy ;/*Yy*/
  tmp[5]= v1[1]*(*sys).AH_zx + v2[1]*(*sys).AH_zy + v3[1]*(*sys).AH_zz ;/*Yz*/
  tmp[6]= v1[2]*(*sys).AH_xx + v2[2]*(*sys).AH_yx + v3[2]*(*sys).AH_zx ;/*Zx*/
  tmp[7]= v1[2]*(*sys).AH_yx + v2[2]*(*sys).AH_yy + v3[2]*(*sys).AH_zy ;/*Zy*/
  tmp[8]= v1[2]*(*sys).AH_zx + v2[2]*(*sys).AH_zy + v3[2]*(*sys).AH_zz ;/*Zz*/


  /* D*AH*Dt */
  (*pol0).xx= tmp[0]*v1[0] + tmp[1]*v2[0] + tmp[2]*v3[0];
  (*pol0).yx= tmp[3]*v1[0] + tmp[4]*v2[0] + tmp[5]*v3[0];
  (*pol0).yy= tmp[3]*v1[1] + tmp[4]*v2[1] + tmp[5]*v3[1];
  (*pol0).zx= tmp[6]*v1[0] + tmp[7]*v2[0] + tmp[8]*v3[0];
  (*pol0).zy= tmp[6]*v1[1] + tmp[7]*v2[1] + tmp[8]*v3[1];
  (*pol0).zz= tmp[6]*v1[2] + tmp[7]*v2[2] + tmp[8]*v3[2];
  
}


/*This subroutine should be removed once the 
structures with the positions will be removed
It just serves to rotate the H atomic dipole moment
from bond referential to lab referential*/
void rot_dip0_H(sys_info *sys, vect_3d *dip0, double *v3){
  (*dip0).x= v3[0]*(*sys).MH_z ;/*X*/
  (*dip0).y= v3[1]*(*sys).MH_z ;/*Y*/
  (*dip0).z= v3[2]*(*sys).MH_z ;/*Z*/

}



