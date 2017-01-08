#ifndef DIP_POL_H
#define DIP_POL_H

void val_init(sys_info *sys, mol_info *mol, vect_3d *dip0, vect_3d *dip0_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, FILE *file_traj, FILE *file_dip0, FILE *file_pol0, int step);
void comp_dip_pol(sys_info *sys, mol_info *mol, double *at_coord, vect_3d *dip0, vect_3d *dip0_mol, vect_3d *dip, vect_3d *dip_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_3d *pol_mol, vect_3d *e_ind, mat_3d *a_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, mat_sym_3d *Tij_mc, mat_sym_3d *Tij_at, double *x, double *y, double *z, double *r, double *fthole, FILE *fileo, FILE *file_Tij_mc, FILE *file_Tij_at, FILE *file_dip, FILE *file_pol, FILE *file_dip_ind, FILE *file_pol_ind, int step);
void sumT_dip(mat_sym_3d *tjk, vect_3d *dip, vect_3d *v3_tmp);
void sumT_pol(mat_sym_3d *tjk, mat_3d *pol, mat_3d *m3_tmp);
void init_dip_pol(vect_3d *dip0, vect_3d *dip, mat_sym_3d *pol0, mat_3d *pol, sys_info *sys);
void calc_tij(sys_info *sys, mol_info *mol, double *at_coord, double *x, double *y, double *z, double *r, double *fthole, mat_sym_3d *Tij_mc, int size,FILE *file_Tij, int step);
void scf_protection(sys_info *sys, mol_info *mol, vect_3d *dip0, vect_3d *dip, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_sym_3d *Tij, vect_3d *e_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, int step);
void thole(double *r, double *fthole, int thole);

#endif 
