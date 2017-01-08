/* These file contain the subroutines which are different because */
/* of the molecular or atomic treatment of the polarizability. */

#include <stdio.h>
#include <stdlib.h>
#include <main.h>
#include <my_math.h>

/*=================================================================
  This subroutine calculates the permanent molecular dipole moment
  and polarizability
  =================================================================*/
void dippol0_molmol(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0){
  double v1[3]={0.}, v2[3]={0.}, v3[3]={0.};

  /*Calculation of the direction cosine matrix*/
  /*  Projection of the molecular axis on the lab referential
      All the vector are normalized
      v3=first bissector of H2O
      v2=vector out of plane (vectorial product of bissector with one OH bond)
      v1=second vector in the molecular plane (vectorial procuct of the first two)*/

  v3[0]=v_oh2[0]+v_oh1[0];
  v3[1]=v_oh2[1]+v_oh1[1];
  v3[2]=v_oh2[2]+v_oh1[2];
  norm(v3,v3);
  cross_p(v3,v_oh1,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
      
  /*
    Calculation of the dipole moment (gas phase-like)
    Initially it is on the bissector
  */
  (*dip0).x=(*sys).M_z*v3[0];
  (*dip0).y=(*sys).M_z*v3[1];
  (*dip0).z=(*sys).M_z*v3[2];

  /*
    Calculation of the polarizability tensor (gas phase-like)
    Initially it is on the bissector
  */
  (*pol0).xx = v1[0]*v1[0]*(*sys).A_xx + v2[0]*v2[0]*(*sys).A_yy + v3[0]*v3[0]*(*sys).A_zz;
  (*pol0).yx = v1[1]*v1[0]*(*sys).A_xx + v2[1]*v2[0]*(*sys).A_yy + v3[1]*v3[0]*(*sys).A_zz;
  (*pol0).yy = v1[1]*v1[1]*(*sys).A_xx + v2[1]*v2[1]*(*sys).A_yy + v3[1]*v3[1]*(*sys).A_zz;
  (*pol0).zx = v1[2]*v1[0]*(*sys).A_xx + v2[2]*v2[0]*(*sys).A_yy + v3[2]*v3[0]*(*sys).A_zz;
  (*pol0).zy = v1[2]*v1[1]*(*sys).A_xx + v2[2]*v2[1]*(*sys).A_yy + v3[2]*v3[1]*(*sys).A_zz;
  (*pol0).zz = v1[2]*v1[2]*(*sys).A_xx + v2[2]*v2[2]*(*sys).A_yy + v3[2]*v3[2]*(*sys).A_zz;
}




/*================================================================
  This subroutine calculates the permanent molecular dipole moment
  and atomic polarizability
  ================================================================*/
void dippol0_molat(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0){
  double v1[3]={0.}, v2[3]={0.}, v3[3]={0.};

  /*Calculation of the direction cosine matrix for the molecule and O*/
  /*  Projection of the molecular axis on the lab referential
      All the vector are normalized
      v3=first bissector of H2O
      v2=vector out of plane (vectorial product of bissector with one OH bond)
      v1=second vector in the molecular plane (vectorial procuct of the first two)*/

  v3[0]=v_oh2[0]+v_oh1[0];
  v3[1]=v_oh2[1]+v_oh1[1];
  v3[2]=v_oh2[2]+v_oh1[2];
  norm(v3,v3);
  cross_p(v3,v_oh1,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
      
  /*
    Calculation of the permanent dipole moment of the whole molecule
    Initially it is on the bissector
  */
  (*dip0).x=(*sys).M_z*v3[0];
  (*dip0).y=(*sys).M_z*v3[1];
  (*dip0).z=(*sys).M_z*v3[2];

  /*
    Calculation of the permanent polarizability tensor of O
    Initially it is on the bissector
  */
  (*(pol0+0)).xx = v1[0]*v1[0]*(*sys).AO_xx + v2[0]*v2[0]*(*sys).AO_yy + v3[0]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).yx = v1[1]*v1[0]*(*sys).AO_xx + v2[1]*v2[0]*(*sys).AO_yy + v3[1]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).yy = v1[1]*v1[1]*(*sys).AO_xx + v2[1]*v2[1]*(*sys).AO_yy + v3[1]*v3[1]*(*sys).AO_zz;
  (*(pol0+0)).zx = v1[2]*v1[0]*(*sys).AO_xx + v2[2]*v2[0]*(*sys).AO_yy + v3[2]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).zy = v1[2]*v1[1]*(*sys).AO_xx + v2[2]*v2[1]*(*sys).AO_yy + v3[2]*v3[1]*(*sys).AO_zz;
  (*(pol0+0)).zz = v1[2]*v1[2]*(*sys).AO_xx + v2[2]*v2[2]*(*sys).AO_yy + v3[2]*v3[2]*(*sys).AO_zz;



  /*Calculation of the direction cosine matrix for H*/
  /*  Projection of the molecular axis on the lab referential
      All the vector are normalized
      v3=O-H bond
      v2=vector out of plane (vectorial product of main OH and the other OH !! IN THIS ORDER !!)
      v1=second vector in the molecular plane (vectorial procuct of the first two)*/

  /* First OH direction cosine matrix */
  norm(v_oh1,v3);
  cross_p(v_oh1,v_oh2,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
  /* First OH Polarizability */
  rot_pol0_H(sys,(pol0+1),v1,v2,v3);


  /* Second OH direction cosine matrix */
  norm(v_oh2,v3);
  cross_p(v_oh2,v_oh1,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
  /* Second OH Polarizability */
  rot_pol0_H(sys,(pol0+2),v1,v2,v3);

}




/*=============================================================
  This subroutine calculates the permanent atomic dipole moment
  and atomic polarizability
  =============================================================*/
void dippol0_atat(sys_info *sys, double *v_oh1, double *v_oh2, vect_3d *dip0, mat_sym_3d *pol0){
  double v1[3]={0.}, v2[3]={0.}, v3[3]={0.};

  /*Calculation of the direction cosine matrix for  O*/
  /*  Projection of the molecular axis on the lab referential
      All the vector are normalized
      v3=first bissector of H2O
      v2=vector out of plane (vectorial product of bissector with one OH bond)
      v1=second vector in the molecular plane (vectorial procuct of the first two)*/

  v3[0]=v_oh2[0]+v_oh1[0];
  v3[1]=v_oh2[1]+v_oh1[1];
  v3[2]=v_oh2[2]+v_oh1[2];
  norm(v3,v3);
  cross_p(v3,v_oh1,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
      
  /* Calculation of the permanent dipole moment of the O atom */
  (*dip0).x=(*sys).MO_z*v3[0];
  (*dip0).y=(*sys).MO_z*v3[1];
  (*dip0).z=(*sys).MO_z*v3[2];

  /*
    Calculation of the permanent polarizability tensor of O
    Initially it is on the bissector
  */
  (*(pol0+0)).xx = v1[0]*v1[0]*(*sys).AO_xx + v2[0]*v2[0]*(*sys).AO_yy + v3[0]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).yx = v1[1]*v1[0]*(*sys).AO_xx + v2[1]*v2[0]*(*sys).AO_yy + v3[1]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).yy = v1[1]*v1[1]*(*sys).AO_xx + v2[1]*v2[1]*(*sys).AO_yy + v3[1]*v3[1]*(*sys).AO_zz;
  (*(pol0+0)).zx = v1[2]*v1[0]*(*sys).AO_xx + v2[2]*v2[0]*(*sys).AO_yy + v3[2]*v3[0]*(*sys).AO_zz;
  (*(pol0+0)).zy = v1[2]*v1[1]*(*sys).AO_xx + v2[2]*v2[1]*(*sys).AO_yy + v3[2]*v3[1]*(*sys).AO_zz;
  (*(pol0+0)).zz = v1[2]*v1[2]*(*sys).AO_xx + v2[2]*v2[2]*(*sys).AO_yy + v3[2]*v3[2]*(*sys).AO_zz;



  /*Calculation of the direction cosine matrix for H*/
  /*  Projection of the molecular axis on the lab referential
      All the vector are normalized
      v3=O-H bond
      v2=vector out of plane (vectorial product of main OH and the other OH !! IN THIS ORDER !!)
      v1=second vector in the molecular plane (vectorial procuct of the first two)*/

  /* First OH direction cosine matrix */
  norm(v_oh1,v3);
  cross_p(v_oh1,v_oh2,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);
  
  rot_dip0_H(sys,(dip0+1),v3);/*Dipole moment*/
  rot_pol0_H(sys,(pol0+1),v1,v2,v3);/*Polarizability*/


  /* Second OH direction cosine matrix */
  norm(v_oh2,v3);
  cross_p(v_oh2,v_oh1,v2);
  norm(v2,v2);
  cross_p(v2,v3,v1);

  rot_dip0_H(sys,(dip0+2),v3);/*Dipole moment*/
  rot_pol0_H(sys,(pol0+2),v1,v2,v3);/*Polarizability*/

}






/* This subroutine sums the atomic contribution to calculate
   the molecular contribution (dip)*/
void dip_at2mol(vect_3d *dip, vect_3d *dip_mol){
  (*dip_mol).x = (*(dip + 0)).x + (*(dip + 1)).x + (*(dip + 2)).x;
  (*dip_mol).y = (*(dip + 0)).y + (*(dip + 1)).y + (*(dip + 2)).y;
  (*dip_mol).z = (*(dip + 0)).z + (*(dip + 1)).z + (*(dip + 2)).z;
}



/* This subroutine sums the atomic contribution to calculate
   the molecular contribution (pol0)*/
void pol0_at2mol(mat_sym_3d *pol0, mat_sym_3d *pol0_mol){
  (*pol0_mol).xx = (*(pol0 + 0)).xx + (*(pol0 + 1)).xx + (*(pol0 + 2)).xx;
  (*pol0_mol).yx = (*(pol0 + 0)).yx + (*(pol0 + 1)).yx + (*(pol0 + 2)).yx;
  (*pol0_mol).yy = (*(pol0 + 0)).yy + (*(pol0 + 1)).yy + (*(pol0 + 2)).yy;
  (*pol0_mol).zx = (*(pol0 + 0)).zx + (*(pol0 + 1)).zx + (*(pol0 + 2)).zx;
  (*pol0_mol).zy = (*(pol0 + 0)).zy + (*(pol0 + 1)).zy + (*(pol0 + 2)).zy;
  (*pol0_mol).zz = (*(pol0 + 0)).zz + (*(pol0 + 1)).zz + (*(pol0 + 2)).zz;
}

/* This subroutine sums the atomic contribution to calculate
   the molecular contribution (pol)*/
void pol_at2mol(mat_3d *pol, mat_3d *pol_mol){
  (*pol_mol).xx = (*(pol + 0)).xx + (*(pol + 1)).xx + (*(pol + 2)).xx;
  (*pol_mol).xy = (*(pol + 0)).xy + (*(pol + 1)).xy + (*(pol + 2)).xy;
  (*pol_mol).xz = (*(pol + 0)).xz + (*(pol + 1)).xz + (*(pol + 2)).xz;
  (*pol_mol).yx = (*(pol + 0)).yx + (*(pol + 1)).yx + (*(pol + 2)).yx;
  (*pol_mol).yy = (*(pol + 0)).yy + (*(pol + 1)).yy + (*(pol + 2)).yy;
  (*pol_mol).yz = (*(pol + 0)).yz + (*(pol + 1)).yz + (*(pol + 2)).yz;
  (*pol_mol).zx = (*(pol + 0)).zx + (*(pol + 1)).zx + (*(pol + 2)).zx;
  (*pol_mol).zy = (*(pol + 0)).zy + (*(pol + 1)).zy + (*(pol + 2)).zy;
  (*pol_mol).zz = (*(pol + 0)).zz + (*(pol + 1)).zz + (*(pol + 2)).zz;
}


