/*
  In this file can be found the subroutines
  which are necessary to calculate the dipole moment
  and the polarizability tensor

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <main.h>
#include <dip_pol.h>
#include <mol_atom.h>
#include <my_math.h>
#include <write.h>


void val_init(sys_info *sys, mol_info *mol, vect_3d *dip0, vect_3d *dip0_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, FILE *file_traj, FILE *file_dip0, FILE *file_pol0, int step){
  int i=0;
  double v_oh1[3]={0.}, v_oh2[3]={0.};
  
  for(i=0 ; i<(*sys).nb_mol ; i++){

#ifndef YUKI
    fprintf(file_traj,"X %9.4f %9.4f %9.4f\n",mol[i].x,mol[i].y,mol[i].z);
#endif
    
    /*
      Regular dipole moment and polarizability for rigid water molecule
    */
    /* Calculation of the O-H vectors*/
    v_oh1[0]=mol[i].xH1-mol[i].xO;
    v_oh1[1]=mol[i].yH1-mol[i].yO;
    v_oh1[2]=mol[i].zH1-mol[i].zO;
    norm(v_oh1,v_oh1);
    
    v_oh2[0]=mol[i].xH2-mol[i].xO;
    v_oh2[1]=mol[i].yH2-mol[i].yO;
    v_oh2[2]=mol[i].zH2-mol[i].zO;
    norm(v_oh2,v_oh2);


    /*==============================================
    The atomic polarizability can be used even
    without atomic dipole: 1 on 3 point independents
    therefore there are 4 possibilities
    ================================================*/
    if((*sys).typ_dip){ /*Atomic dipole moment*/
      if((*sys).typ_pol){ /* Atomic polarizability */
	dippol0_atat(sys,v_oh1,v_oh2,&(dip0[3*i]),&(pol0[3*i]));
	dip_at2mol(&(dip0[3*i]),&(dip0_mol[i]));
	pol0_at2mol(&(pol0[3*i]),&(pol0_mol[i]));
      }
      else{ /* Molecular polarizability */
	printf("I thought that I had blocked this possibility.\n");
	printf("End of program.\n");
	exit(0);
      }
    }
    else{/*Molecular dipole moment */
      if((*sys).typ_pol){ /* Atomic polarizability */
	dippol0_molat(sys,v_oh1,v_oh2,&(dip0[i]),&(pol0[3*i]));
	dip0_mol[i]=dip0[i];
	pol0_at2mol(&(pol0[3*i]),&(pol0_mol[i]));
      }
      else{ /* Molecular polarizability */
	dippol0_molmol(sys,v_oh1,v_oh2,&(dip0[i]),&(pol0[i]));
	dip0_mol[i]=dip0[i];
	pol0_mol[i]=pol0[i];
      }
    }
    
  }

  

    
#ifdef DEBUG
  for(i=0 ; i<(*sys).nb_dip ; i++){
    fprintf(file_dip0,"%10d %5d %9.4f %9.4f %9.4f\n",step,i,dip0[i].x,dip0[i].y,dip0[i].z);
  }
  for(i=0 ; i<(*sys).nb_pol ; i++){
    fprintf(file_pol0,"%10d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",step,i,pol0[i].xx,pol0[i].yx,pol0[i].yy,	\
	    pol0[i].zx,pol0[i].zy,pol0[i].zz);
  }
#endif
  
  return;
}



/*
  Complete dipole (permanent + induced) and effective polarizbility
  by Morita and Hynes, J. Phys. Chem. B 2002, 106, 673-685
  and calculated only once (Yuki version)
*/
void comp_dip_pol(sys_info *sys, mol_info *mol, double *at_coord, vect_3d *dip0, vect_3d *dip0_mol, vect_3d *dip, vect_3d *dip_mol, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_3d *pol_mol, vect_3d *e_ind, mat_3d *a_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, mat_sym_3d *Tij_mc, mat_sym_3d *Tij_at, double *x, double *y, double *z, double *r, double *fthole, FILE *fileo, FILE *file_Tij_mc, FILE *file_Tij_at, FILE *file_dip, FILE *file_pol, FILE *file_dip_ind, FILE *file_pol_ind, int step){
  int i=0 , j=0 , ij=0, test_scf=0;

  /* Calculation of the dipole field tensor Tij*/
  calc_tij(sys,mol,at_coord,x,y,z,r,fthole,Tij_mc,(*sys).nb_mol,file_Tij_mc,step);/*Molecular*/
  if((*sys).typ_dip+(*sys).typ_pol){/*Atomic*/
    calc_tij(sys,mol,at_coord,x,y,z,r,fthole,Tij_at,(*sys).nb_point,file_Tij_at,step);
  }
  
  /* Initialization of the total terms*/
  init_dip_pol(dip0,dip,pol0,pol,sys);
  /* Initialization of the induced terms at the order 0*/
  init_dip_pol(dip0,e_ind,pol0,a_ind,sys);
  
  /*========================*/
  /*========SCF loop========*/
  /*========================*/
  test_scf=1;/*To know when we go in or out of the while loop*/
  while(test_scf){
    /*Initialization*/
    test_scf=0;
    for(i=0;i<(*sys).nb_dip;i++){
      v3_tmp[i].x=0.; v3_tmp[i].y=0.; v3_tmp[i].z=0.;
    }
    for(i=0;i<(*sys).nb_pol;i++){
      m3_tmp[i].xx=0.; m3_tmp[i].xy=0.; m3_tmp[i].xz=0.;
      m3_tmp[i].yx=0.; m3_tmp[i].yy=0.; m3_tmp[i].yz=0.;
      m3_tmp[i].zx=0.; m3_tmp[i].zy=0.; m3_tmp[i].zz=0.;
    }
      

    /*-Sum Tij.e_ind[j], -Sum Tij.a_ind[j]*/
    /*-Sum Tji.e_ind[i], -Sum Tji.a_ind[i]*/
    /* Dipole*/
    ij=0;
    if((*sys).typ_dip){/*Atomic*/
      for(i=0;i<(*sys).nb_dip;i++){
	for(j=i+1;j<(*sys).nb_dip;j++){
	  sumT_dip(&(Tij_at[ij]),&(e_ind[j]),&(v3_tmp[i]));
	  sumT_dip(&(Tij_at[ij]),&(e_ind[i]),&(v3_tmp[j]));
	  ij++;
	}
      }
    }
    else{/*Molecular*/
      for(i=0;i<(*sys).nb_dip;i++){
	for(j=i+1;j<(*sys).nb_dip;j++){
	  sumT_dip(&(Tij_mc[ij]),&(e_ind[j]),&(v3_tmp[i]));
	  sumT_dip(&(Tij_mc[ij]),&(e_ind[i]),&(v3_tmp[j]));
	  ij++;
	}
      }
    }
    
    /*Polarizability*/
    ij=0;
    if((*sys).typ_pol){/*Atomic*/
      for(i=0;i<(*sys).nb_pol;i++){
	for(j=i+1;j<(*sys).nb_pol;j++){
	  sumT_pol(&(Tij_at[ij]),&(a_ind[j]),&(m3_tmp[i]));
	  sumT_pol(&(Tij_at[ij]),&(a_ind[i]),&(m3_tmp[j]));
	  ij++;
	}
      }
    }
    else{/*Molecular*/
      for(i=0;i<(*sys).nb_pol;i++){
	for(j=i+1;j<(*sys).nb_pol;j++){
	  sumT_pol(&(Tij_mc[ij]),&(a_ind[j]),&(m3_tmp[i]));
	  sumT_pol(&(Tij_mc[ij]),&(a_ind[i]),&(m3_tmp[j]));
	  ij++;
	}
      }
    }

    /*==========================================================*/
    /* Calculation of the induced dipole moment at the order N+1*/
    /*==========================================================*/
    if((*sys).typ_dip){/*Atomic*/
      for(i=0;i<(*sys).nb_dip;i++){
	/*If atomic dipole, therefore there are atomic polarizabilities
	  So it is sure that pol0[i] is defined*/
	mult_msym_v_3d(&(pol0[i]),&(v3_tmp[i]),&(e_ind[i]));
      }
    }
    else{/*Molecular*/
      for(i=0;i<(*sys).nb_dip;i++){
	mult_msym_v_3d(&(pol0_mol[i]),&(v3_tmp[i]),&(e_ind[i]));
      }
    }
      
    for(i=0;i<(*sys).nb_dip;i++){
      /*Update of the complete dipole at the order N+1*/
      dip[i].x += e_ind[i].x;
      dip[i].y += e_ind[i].y;
      dip[i].z += e_ind[i].z;
      
      /*The norm of the induced dipole is not converged*/
      if(sqrt(pow(e_ind[i].x,2)+pow(e_ind[i].y,2)+pow(e_ind[i].z,2))>DIP_CONV){
	test_scf=1;
      }

#ifdef DEBUG
      fprintf(file_dip_ind,"%10d %5d %9.4f %9.4f %9.4f\n",step,i,e_ind[i].x,e_ind[i].y,e_ind[i].z);
      fprintf(file_dip,"%10d %5d %9.4f %9.4f %9.4f\n",step,i,dip[i].x,dip[i].y,dip[i].z);
#endif

    }

    /*===========================================================*/
    /* Calculation of the induced polarizability at the order N+1*/
    /*===========================================================*/
    for(i=0;i<(*sys).nb_pol;i++){
      mult_msym_m_3d(&(pol0[i]),&(m3_tmp[i]),&(a_ind[i]));

      /*Update of the complete polarizability at the order N+1*/
      pol[i].xx += a_ind[i].xx ; pol[i].xy += a_ind[i].xy ; pol[i].xz += a_ind[i].xz ;
      pol[i].yx += a_ind[i].yx ; pol[i].yy += a_ind[i].yy ; pol[i].yz += a_ind[i].yz ;
      pol[i].zx += a_ind[i].zx ; pol[i].zy += a_ind[i].zy ; pol[i].zz += a_ind[i].zz ;

      /*The trace of the induced polarizability is not converged*/
      if(sqrt(pow(a_ind[i].xx+a_ind[i].yy+a_ind[i].zz,2))>POL_CONV ){
	test_scf=1;
      }

#ifdef DEBUG
      fprintf(file_pol_ind,"%10d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",step,i, \
	      a_ind[i].xx,a_ind[i].xy,a_ind[i].xz,			\
	      a_ind[i].yx,a_ind[i].yy,a_ind[i].yz,			\
	      a_ind[i].zx,a_ind[i].zy,a_ind[i].zz);
      fprintf(file_pol,"%10d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",step,i, \
	      pol[i].xx,pol[i].xy,pol[i].xz,				\
	      pol[i].yx,pol[i].yy,pol[i].yz,				\
	      pol[i].zx,pol[i].zy,pol[i].zz);
#endif

      /*Protection against infinite loop
	If we are out of the limits, only the first order dipole
	moment and polarizability are written for the whole step*/
      if(pol[i].xx>POL_LIM ||pol[i].xy>POL_LIM ||pol[i].xz>POL_LIM ||\
	 pol[i].yx>POL_LIM ||pol[i].yy>POL_LIM ||pol[i].yz>POL_LIM ||	\
	 pol[i].zx>POL_LIM ||pol[i].zy>POL_LIM ||pol[i].zz>POL_LIM){
	scf_protection(sys,mol,dip0,dip,pol0,pol0_mol,pol,Tij_mc,e_ind,v3_tmp,m3_tmp,step);
	i=(*sys).nb_mol;/*Out of the for loop*/
	test_scf=0;/*Out of the while loop*/
      }

    }/* Induced polarizability*/
  }/*SCF*/


  /*===================================
    Calculation of the molecular values
    ===================================*/
  /* Dipole moment */
  if((*sys).typ_dip){
    for(i=0;i<(*sys).nb_mol;i++){
      dip_at2mol(&(dip[3*i]),&(dip_mol[i]));
    }
  }
  else{
    dip_mol=dip;
  }

  /* Polarizability */
  if((*sys).typ_pol){ 
    for(i=0;i<(*sys).nb_mol;i++){
      pol_at2mol(&(pol[3*i]),&(pol_mol[i]));
    }
  }
  else{
    pol_mol=pol;
  }

  /* We write the info about the dipole moment and the polarizability */
  write_dippol(*sys,dip0_mol,dip_mol,pol0_mol,pol_mol,fileo);

  return;
}



/* Calculation of the Sum involving Tjk and dip[k]
   Preliminary step before the calculation of the induced dipole

   Remark1: The sign "-" is already included 
   Remark2: tjk is diagonal (only lower part)

*/
void sumT_dip(mat_sym_3d *tjk, vect_3d *dip, vect_3d *v3_tmp){
  
  /* Contribution of k on the induced field of j */
  (*v3_tmp).x -= (*tjk).xx*(*dip).x + (*tjk).yx*(*dip).y + (*tjk).zx*(*dip).z;
  (*v3_tmp).y -= (*tjk).yx*(*dip).x + (*tjk).yy*(*dip).y + (*tjk).zy*(*dip).z;
  (*v3_tmp).z -= (*tjk).zx*(*dip).x + (*tjk).zy*(*dip).y + (*tjk).zz*(*dip).z;

  return;
}




/* Calculation of the Sum involving Tjk and pol[k]
   Preliminary step before the calculation of the induced polarizability 

   Remark1: The sign "-" is already included 
   Remark2: tjk is diagonal (only lower part)

*/
void sumT_pol(mat_sym_3d *tjk, mat_3d *pol, mat_3d *m3_tmp){
  /*Diagonal terms of j*/
  (*m3_tmp).xx -= (*tjk).xx*(*pol).xx + (*tjk).yx*(*pol).yx + (*tjk).zx*(*pol).zx;
  (*m3_tmp).yy -= (*tjk).yx*(*pol).xy + (*tjk).yy*(*pol).yy + (*tjk).zy*(*pol).zy;
  (*m3_tmp).zz -= (*tjk).zx*(*pol).xz + (*tjk).zy*(*pol).yz + (*tjk).zz*(*pol).zz;

  /*Lower part of j*/
  (*m3_tmp).yx -= (*tjk).yx*(*pol).xx + (*tjk).yy*(*pol).yx + (*tjk).zy*(*pol).zx;
  (*m3_tmp).zx -= (*tjk).zx*(*pol).xx + (*tjk).zy*(*pol).yx + (*tjk).zz*(*pol).zx;
  (*m3_tmp).zy -= (*tjk).zx*(*pol).xy + (*tjk).zy*(*pol).yy + (*tjk).zz*(*pol).zy;

  /*Upper part of j*/
  (*m3_tmp).xy -= (*tjk).xx*(*pol).xy + (*tjk).yx*(*pol).yy + (*tjk).zx*(*pol).zy;
  (*m3_tmp).xz -= (*tjk).xx*(*pol).xz + (*tjk).yx*(*pol).yz + (*tjk).zx*(*pol).zz;
  (*m3_tmp).yz -= (*tjk).yx*(*pol).xz + (*tjk).yy*(*pol).yz + (*tjk).zy*(*pol).zz;

  return;
}


/*Initialization before the calculation of 
  the dipole moment and the polarizability*/
void init_dip_pol(vect_3d *dip0, vect_3d *dip, mat_sym_3d *pol0, mat_3d *pol, sys_info *sys){
  int i=0;

  for(i=0;i<(*sys).nb_dip;i++){
    dip[i].x=dip0[i].x; dip[i].y=dip0[i].y; dip[i].z=dip0[i].z;
  }
  for(i=0;i<(*sys).nb_pol;i++){
    pol[i].xx=pol0[i].xx; pol[i].xy=pol0[i].yx; pol[i].xz=pol0[i].zx;
    pol[i].yx=pol0[i].yx; pol[i].yy=pol0[i].yy; pol[i].yz=pol0[i].zy;
    pol[i].zx=pol0[i].zx; pol[i].zy=pol0[i].zy; pol[i].zz=pol0[i].zz;
  }

  return;
}


/* Calculation of the dipole field tensor Tij
   Only for k>j because Tjk=Tkj and Tjj=0 */
void calc_tij(sys_info *sys, mol_info *mol, double *at_coord, double *x, double *y, double *z, double *r, double *fthole, mat_sym_3d *Tij, int size, FILE *file_Tij, int step){
  int i=0, j=0, ij=0;
  double r3=0., r5=0.;
  
  /*Calculation of the distance between the atoms i and j*/
  if(size==(*sys).nb_mol){ /*Mass center*/
    for(i=0;i<size;i++){
      for(j=i+1;j<size;j++){
	x[ij] = mol[i].x - mol[j].x;
	y[ij] = mol[i].y - mol[j].y;
	z[ij] = mol[i].z - mol[j].z;
	ij++;
      }
    }
  }
  else{ /*Atomic*/
    for(i=0;i<size;i++){
      for(j=i+1;j<size;j++){
	x[ij] = at_coord[3*i+0] - at_coord[3*j+0];
	y[ij] = at_coord[3*i+1] - at_coord[3*j+1];
	z[ij] = at_coord[3*i+2] - at_coord[3*j+2];
	ij++;
      }
    }
  }
  for(i=0;i<size*(size-1)/2;i++){
    pbc((*sys).cell,&x[i],&y[i],&z[i]);
    r[i]=sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
  }
  
  /* Calculation of the screening functions (Thole)
     The parameters are already set to 1.d0 if we do not
     want these screening function
  */
  if((*sys).thole){
    thole(r,fthole,size);
  }

  /* Calculation of Tij */
  ij=0;
  for(i=0;i<size;i++){
    for(j=i+1;j<size;j++){

      if(r[ij]<=(*sys).cutoff){

	r3=pow(r[ij],-3.)* fthole[2*ij  ]; /*This is 1/r3 with screening function*/
	r5=3/pow(r[ij],5.)*fthole[2*ij+1]; /*This is 3/r5 with screening function*/

	/*Tij (=Tji) calculation*/
	Tij[ij].xx=r3-r5*x[ij]*x[ij];
	Tij[ij].yx=  -r5*y[ij]*x[ij]; Tij[ij].yy=r3-r5*y[ij]*y[ij];
	Tij[ij].zx=  -r5*z[ij]*x[ij]; Tij[ij].zy=  -r5*z[ij]*y[ij]; Tij[ij].zz=r3-r5*z[ij]*z[ij];
	
      }
      else{
	Tij[ij].xx=0.0;
	Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
	Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;
      }

      ij++;
    }
  }

  /*If we do not want intra molecular DID (for 3 points model only)*/
  if(size>(*sys).nb_mol){
    if((*sys).intra==0){
      ij=0;
      for(i=0;i<size;i=i+3){/*Only the O are passed*/
	/*No interaction between O and H1*/
	Tij[ij].xx=0.0;
	Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
	Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;

	/*No interaction between O and H2*/
	ij++;
	Tij[ij].xx=0.0;
	Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
	Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;

	/*No interaction between H1 and H2*/
	ij=ij+size-i-2;
	Tij[ij].xx=0.0;
	Tij[ij].yx=0.0;	Tij[ij].yy=0.0;
	Tij[ij].zx=0.0;	Tij[ij].zy=0.0;	Tij[ij].zz=0.0;

	ij=ij+2*size-2*i-5; /*The intermolecular interactions are passed*/
      }
    }
  }
  
#ifdef DEBUG
  ij=0;
  for(i=0;i<size;i++){
    for(j=i+1;j<size;j++){
      fprintf(file_Tij,"%10d %10d %5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",step,i,j, \
	      Tij[ij].xx,Tij[ij].yx,Tij[ij].yy,Tij[ij].zx,Tij[ij].zy,Tij[ij].zz);
      ij++;
    }
  }
#endif
  
  return;
}


/*Protection against infinite loop
  If we are out of the limits, only the first order dipole
  moment and polarizability are written for the whole step
  XXX TRAVAIL ARBEIT WORK
  Check if it works
*/
void scf_protection(sys_info *sys, mol_info *mol, vect_3d *dip0, vect_3d *dip, mat_sym_3d *pol0, mat_sym_3d *pol0_mol, mat_3d *pol, mat_sym_3d *Tij_mc, vect_3d *e_ind, vect_3d *v3_tmp, mat_3d *m3_tmp, int step){
  int i=0, j=0, ij=0;
  FILE *file_pb;

  /*Initialization*/
  init_dip_pol(dip0,dip,pol0,pol,sys);
  for(i=0;i<(*sys).nb_mol;i++){
    v3_tmp[i].x=0.; v3_tmp[i].y=0.; v3_tmp[i].z=0.;
  }
  for(i=0;i<(*sys).nb_pol;i++){
    m3_tmp[i].xx=0.; m3_tmp[i].xy=0.; m3_tmp[i].xz=0.;
    m3_tmp[i].yx=0.; m3_tmp[i].yy=0.; m3_tmp[i].yz=0.;
    m3_tmp[i].zx=0.; m3_tmp[i].zy=0.; m3_tmp[i].zz=0.;
  }
  
  /*-Sum Tij.dip[j], -Sum Tij.pol[j]*/
  /*-Sum Tji.dip[i], -Sum Tji.pol[i] */
  ij=0;
  for(i=0;i<(*sys).nb_mol;i++){
    for(j=i+1;j<(*sys).nb_mol;j++){
      sumT_dip(&(Tij_mc[ij]),&(dip[j]),&(v3_tmp[i]));
      sumT_dip(&(Tij_mc[ij]),&(dip[i]),&(v3_tmp[j]));
      ij++;
    }
  }
  ij=0;
  for(i=0;i<(*sys).nb_point;i++){
    for(j=i+1;j<(*sys).nb_point;j++){
      sumT_pol(&(Tij_mc[ij]),&(pol[j]),&(m3_tmp[i]));
      sumT_pol(&(Tij_mc[ij]),&(pol[i]),&(m3_tmp[j]));
      ij++;
    }
  }

  
  /*1-Sum Tij.pol[j]*/
  for(i=0;i<(*sys).nb_point;i++){
    m3_tmp[i].xx+=1.;
    m3_tmp[i].yy+=1.;
    m3_tmp[i].zz+=1.;
  }

  for(i=0;i<(*sys).nb_mol;i++){
    mult_msym_v_3d(&(pol0_mol[i]),&(v3_tmp[i]),&(e_ind[i]));/*INDUCED dipole moment*/
    /*Total dipole*/
    dip[i].x = dip0[i].x + e_ind[i].x;
    dip[i].y = dip0[i].y + e_ind[i].y;
    dip[i].z = dip0[i].z + e_ind[i].z;
  }
  for(i=0;i<(*sys).nb_point;i++){
    mult_msym_m_3d(&(pol0[i]),&(m3_tmp[i]),&(pol[i]));/*TOTAL polarizability*/
  }


  /*Record of the positions with a problem of convergence
    for further debugging*/
  FOPEN_SAFE(file_pb,CONV_PB,"a");
  fprintf(file_pb,"%d\n",(*sys).nb_mol);
  fprintf(file_pb,"i = %d , this is the positions which had some convergence problem\n",step);
  for(i=0;i<(*sys).nb_mol;i++){
    fprintf(file_pb," X  %f %f %f\n",mol[i].x,mol[i].y,mol[i].z);
  }
  fclose(file_pb);
}



/*This subroutine returns the 2 screening functions*/
void thole(double *r, double *fthole, int size){
  int i=0, j=0, ij=0;
  double x=0;
  
  for(i=0;i<size;i++){
    for(j=i+1;j<size;j++){
      x=A_THOLE*r[ij];
      
      /* Screening of the term in 1/r3 */
      fthole[2*ij  ]=1-(1 + x + x*x/2.)*exp(-x);
      /* Screening of the term in 1/r5 */
      fthole[2*ij+1]=1-(1 + x + x*x/2. + pow(x,3.)/6)*exp(-x);

      ij++;
    }
  }


};


