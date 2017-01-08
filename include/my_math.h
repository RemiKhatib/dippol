#ifndef MY_MATH_H
#define MY_MATH_H

void cross_p(double *v1, double *v2, double *v3);
void norm(double *v1, double *v2);
void pbc(vect_3d cell, double *x, double *y, double *z);
void pbc_mol(vect_3d cell, mol_info *mol);
void mult_msym_v_3d(mat_sym_3d *m, vect_3d *vin, vect_3d *vout);
void mult_msym_m_3d(mat_sym_3d *m1, mat_3d *m2, mat_3d *mout);
void mult_m_m_3d(mat_3d *m1, mat_3d *m2, mat_3d *mout);
void rot_pol0_H(sys_info *sys, mat_sym_3d *pol0, double *v1, double *v2, double *v3);
void rot_dip0_H(sys_info *sys, vect_3d *mol0, double *v3);

#endif 
