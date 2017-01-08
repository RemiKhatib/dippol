#ifndef MOL_ATOM_H
#define MOL_ATOM_H

void dippol0_molmol(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0);
void dippol0_molat(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0);
void dippol0_atat(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0);
void pol0_at2mol(mat_sym_3d *pol0, mat_sym_3d *pol0_mol);
void dip_at2mol(vect_3d *dip, vect_3d *dip_mol);
void pol_at2mol(mat_3d *pol, mat_3d *pol_mol);
  
#endif
