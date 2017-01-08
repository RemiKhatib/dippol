/*
  In this file can be found the subroutines
  which are necessary with the input files
  (Either from the user or from CP2K)

  I suppose that the xyz-file comes from CP2K
  and therefore is formated!!!!!!
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <main.h>
#include <my_math.h>
#include <treat_inp.h>
#include <dip_pol.h>

/*==================================*/
/*Reading the input file of the user*/
/*==================================*/
void read_input(input_info *input, sys_info *sys, int argc, char *argv[]){
  int i=0, first_H=0;
  char symb=0;
  char chain[CHAIN_SIZE]="",p[CHAIN_SIZE]="",q[CHAIN_SIZE]="",r[CHAIN_SIZE]="";
  FILE *file=NULL;
  
  if(argc!=2){
    printf("You have to specify the name of the input file:\n");
    printf("%s input_file\n",argv[0]);
    printf("\nEnd of program\n");
    exit(0);
  }
  FOPEN_SAFE(file,argv[1],"r");

  /* Position file */
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%s",(*input).namei);
  
  /* Output file (dippol.dat) */
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%s",(*input).nameo);

#ifndef YUKI
  /*Output file for the mass center position*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%s",(*input).name_mass);
#endif
  
  /*Size of the box*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%lf %lf %lf",&((*sys).cell.x),&((*sys).cell.y),&((*sys).cell.z));

  /*Initial / Final step*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%d %d",&((*input).stepi),&((*input).stepf));

  /*Order of the atoms (only the first character will be read)*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%d",&((*sys).at_p_mol));/*Atoms per molecule*/
  MALLOC_SAFE((*sys).order,(*sys).at_p_mol,int);

  first_H=1;
  for(i=0;i<(*sys).at_p_mol;i++){
    symb=fgetc(file);
    fgets(chain,CHAIN_SIZE,file);

    if(symb=='O'){(*sys).order[i]=0;}
    else if(symb=='H'){
      if(first_H){
	(*sys).order[i]=3;/*H1*/
	first_H=0;
      }
      else{
	(*sys).order[i]=6;/*H2*/
      }
    }
    else{(*sys).order[i]=9;}/*Any dummy atom*/
  }

  /*
    Polarization of the beam
    For (*sys).p (*sys).q (*sys).r:
    X=0, Y=1, Z=2
  */
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%s %s %s",p,q,r);
  if(strcmp(r,"X")==0){(*sys).r=0;}
  else if(strcmp(r,"Y")==0){(*sys).r=1;}
  else if(strcmp(r,"Z")==0){(*sys).r=2;}
  else{printf("Problem with the polarization %s\nEnd of program.\n",r);exit(0);}

  if(strcmp(q,"X")==0){(*sys).q=0;}
  else if(strcmp(q,"Y")==0){(*sys).q=1;}
  else if(strcmp(q,"Z")==0){(*sys).q=2;}
  else if(strcmp(q,"S")==0){(*sys).q=-(*sys).r-1;}
  else{printf("Problem with the polarization %s\nEnd of program.\n",q);exit(0);}

  if(strcmp(p,"X")==0){(*sys).p=0;}
  else if(strcmp(p,"Y")==0){(*sys).p=1;}
  else if(strcmp(p,"Z")==0){(*sys).p=2;}
  else if(strcmp(p,"S")==0){
    (*sys).p=-(*sys).r-1;
    /* Extra control for SS {X,Y,Z} polarization */
    if(strcmp(p,q)!=0){
      printf("Problem with the polarization\nYou certainly wanted to write \"S S %s\"\nEnd of program.\n",r);exit(0);
    }
  }
  else{printf("Problem with the polarization %s\nEnd of program.\n",p);exit(0);}
  
  /*Cutoff for the induction*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%lf",&((*sys).cutoff));

  /*Do we use the Thole damping factor for Tij?*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%d",&((*sys).thole));

  /*Atomic or molecular dipole moment?*/
  symb=fgetc(file);
  fgets(chain,CHAIN_SIZE,file);
  if(symb=='a'||symb=='A'){
    (*sys).typ_dip=1;/* 1=Atomic dipole moment */
    fgets(chain,CHAIN_SIZE,file);
    sscanf(chain,"%lf %lf",&((*sys).MO_z),&((*sys).MH_z));
  }
  else if(symb=='m'||symb=='M'){
    (*sys).typ_dip=0;/* 0=Molecular dipole moment */
    fgets(chain,CHAIN_SIZE,file);
    sscanf(chain,"%lf %lf %lf",&((*sys).M_x),&((*sys).M_y),&((*sys).M_z));
  }
  else{
    printf("Problem with the kind of dipole moment. Please chose \"A(a)tomic\" or \"M(m)olecular\"\nEnd of program.\n");exit(0);
  }

  /*Atomic or molecular polarizability?*/
  symb=fgetc(file);
  fgets(chain,CHAIN_SIZE,file);
  if(symb=='a'||symb=='A'){
    (*sys).typ_pol=1;/* 1=Atomic polarizability */
    fgets(chain,CHAIN_SIZE,file);
    sscanf(chain,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",			\
	   &((*sys).AO_xx),&((*sys).AO_yy),&((*sys).AO_zz),		\
	   &((*sys).AH_xx),&((*sys).AH_yy),&((*sys).AH_zz),		\
	   &((*sys).AH_yx),&((*sys).AH_zx),&((*sys).AH_zy));
  }
  else if(symb=='m'||symb=='M'){
    if((*sys).typ_dip){
      printf("What you are doing is not possible\n");
      printf("To have the ATOMIC induced dipole moment,\n");
      printf("you need an ATOMIC polarizability\n");
      printf("alpha0i * Tij * mu_j  (with j!=i)\n");
      printf("End of program\n");
      exit(0);
    }
    (*sys).typ_pol=0;/* 0=Molecular polarizability */
    fgets(chain,CHAIN_SIZE,file);
    sscanf(chain,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",		\
	   &((*sys).A_xx),&((*sys).A_xy),&((*sys).A_xz),	\
	   &((*sys).A_yx),&((*sys).A_yy),&((*sys).A_zz),	\
	   &((*sys).A_zx),&((*sys).A_zy),&((*sys).A_zz));
  }
  else{
    printf("Problem with the kind of polarizability. Please chose \"A(a)tomic\" or \"M(m)olecular\"\nEnd of program.\n");exit(0);
  }

  /*Do we take into account the intra molecular DID?*/
  fgets(chain,CHAIN_SIZE,file);
  sscanf(chain,"%d",&((*sys).intra));
  
  fclose(file);
}




/*==============================*/
/*Reading of the trajectory file*/
/*==============================*/
void trajectory(sys_info *sys, input_info input, char *argv[]) {
  FILE *filei=NULL, *file_debug=NULL;
  int step=0, i=0, j=0, start_read=0, end_read=0;
  double time=0.;
  char chain[CHAIN_SIZE]="";
  
  double *x=NULL, *y=NULL, *z=NULL, *r=NULL, *fthole=NULL, *at_coord;
  vect_3d *dip0=NULL,*dip0_mol=NULL, *dip=NULL, *dip_mol=NULL, *e_ind=NULL, *v3_tmp=NULL;
  mol_info *mol=NULL;
  mat_sym_3d *pol0=NULL, *pol0_mol=NULL, *Tij_mc=NULL, *Tij_at=NULL;
  mat_3d *pol=NULL, *pol_mol=NULL, *a_ind=NULL, *m3_tmp=NULL;
  FILE *fileo=NULL, *file_traj=NULL,*file_dip0=NULL,*file_pol0=NULL;
  FILE *file_Tij_mc=NULL, *file_Tij_at=NULL, *file_dip=NULL,*file_pol=NULL,*file_dip_ind=NULL, *file_pol_ind=NULL;
    
  /*To facilitate the reading, I extract the info from structure input*/
  int stepi=input.stepi, stepf=input.stepf;
  char namei[CHAIN_SIZE]="";
  strcpy(namei,input.namei);
  
  FOPEN_SAFE(filei,namei,"r");
  FOPEN_SAFE(fileo,input.nameo,"w+");
#ifndef YUKI
  FOPEN_SAFE(file_traj,input.name_mass,"w+");
  fprintf(fileo,"#This file has been written by %s. Edit it at your own risk.\n",argv[0]);
  fprintf(fileo,"#Mass_center_file %s\n",input.name_mass);
  fprintf(fileo,"#Cell_parameters %f %f %f\n",(*sys).cell.x,(*sys).cell.y,(*sys).cell.z);
  fprintf(fileo,"#Cutoff_Ang %f\n",(*sys).cutoff);
  fprintf(fileo,"#Dip0 Dip_comp Pol0 Pol_comp\n");
#endif
  
  /*The trajectory file is read while we have not reached EOF
    or the maximal value authorized by the user */
  while(fscanf(filei,"%d",&((*sys).nb_at)) != EOF && !end_read){
    /*Comment line reading*/
    fscanf(filei,"%s",chain);/* Pass the i */
    fscanf(filei,"%s",chain);/* Pass the = */
    fscanf(filei,"%d",&step);
    fscanf(filei,"%s",chain);/*Pass the ,*/
    fscanf(filei,"%s",chain);/*Pass the time*/
    fscanf(filei,"%s",chain);/*Pass the =*/
    fscanf(filei,"%lf",&time);
  
    /*========================*/
    /*First step selected only*/
    /*========================*/
    if(step>=stepi && !start_read){
      start_read=1;

      /*Preparation of the different values requiered for the allocation*/
      (*sys).nb_mol  =(*sys).nb_at/(*sys).at_p_mol;
      (*sys).nb_point=(*sys).nb_at/(*sys).at_p_mol;
      (*sys).nb_dip  =(*sys).nb_mol;
      (*sys).nb_pol  =(*sys).nb_mol;
      if((*sys).typ_dip){/*Atomic dipole*/
	(*sys).nb_dip  =3*(*sys).nb_mol;
	(*sys).nb_point=3*(*sys).nb_mol;
      }
      if((*sys).typ_pol){/*Atomic polarizability*/
	(*sys).nb_pol  =3*(*sys).nb_mol;
	(*sys).nb_point=3*(*sys).nb_mol;
      }

      /*Creation of the different arrays and files*/
      if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
	MALLOC_SAFE(Tij_at,(*sys).nb_point*((*sys).nb_point-1)/2,mat_sym_3d);
      }
      MALLOC_SAFE(mol     ,(*sys).nb_mol                        ,mol_info);
      MALLOC_SAFE(at_coord,9*(*sys).nb_mol                      ,double);
      MALLOC_SAFE(dip0    ,(*sys).nb_dip                        ,vect_3d);
      MALLOC_SAFE(dip     ,(*sys).nb_dip                        ,vect_3d);
      MALLOC_SAFE(dip0_mol,(*sys).nb_mol                        ,vect_3d);
      MALLOC_SAFE(dip_mol ,(*sys).nb_mol                        ,vect_3d);
      MALLOC_SAFE(pol0    ,(*sys).nb_pol                        ,mat_sym_3d);
      MALLOC_SAFE(pol     ,(*sys).nb_pol                        ,mat_3d);
      MALLOC_SAFE(pol0_mol,(*sys).nb_mol                        ,mat_sym_3d);
      MALLOC_SAFE(pol_mol ,(*sys).nb_mol                        ,mat_3d);
      MALLOC_SAFE(e_ind   ,(*sys).nb_dip                        ,vect_3d);
      MALLOC_SAFE(v3_tmp  ,(*sys).nb_dip                        ,vect_3d);
      MALLOC_SAFE(a_ind   ,(*sys).nb_pol                        ,mat_3d);
      MALLOC_SAFE(m3_tmp  ,(*sys).nb_pol                        ,mat_3d);
      MALLOC_SAFE(Tij_mc  ,(*sys).nb_mol*((*sys).nb_mol-1)/2    ,mat_sym_3d);
      MALLOC_SAFE(x       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
      MALLOC_SAFE(y       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
      MALLOC_SAFE(z       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
      MALLOC_SAFE(r       ,(*sys).nb_point*((*sys).nb_point-1)/2,double);
      MALLOC_SAFE(fthole  ,(*sys).nb_point*((*sys).nb_point-1)  ,double);

      /* Initialization of thole screening factor to 1 (No thole factor)*/
      for(i=0;i<(*sys).nb_point*((*sys).nb_point-1);i++){
	fthole[i]=1.;
      }


#ifdef DEBUG
      FOPEN_SAFE(file_debug,DEBUG_TRAJ,"w+");
      FOPEN_SAFE(file_dip0,DEBUG_DIP0,"w+");
      FOPEN_SAFE(file_pol0,DEBUG_POL0,"w+");
      FOPEN_SAFE(file_Tij_mc ,DEBUG_TIJ_MC ,"w+");
      FOPEN_SAFE(file_dip_ind,DEBUG_DIP_IND,"w+");
      FOPEN_SAFE(file_dip    ,DEBUG_DIP    ,"w+");
      FOPEN_SAFE(file_pol    ,DEBUG_POL    ,"w+");
      FOPEN_SAFE(file_pol_ind,DEBUG_POL_IND,"w+");
      fprintf(file_debug,"#Step_number Record_number Mol_number atom_number x y z\n");
      fprintf(file_dip0,"#Record_number Mol_number dip0_x dip0_y dip0_z\n");
      fprintf(file_pol0,"#Record_number Mol_number pol0_xx pol0_yx pol0_yy pol0_zx pol0_zy pol0_zz\n");
      fprintf(file_Tij_mc,"#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n");
      fprintf(file_dip_ind,"#Record_number Mol_number e_ind_x e_ind_y e_ind_z at the order N\n");
      fprintf(file_dip,"#Record_number Mol_number dip_x dip_y dip_z\n");
      fprintf(file_pol,"#Record_number Mol_number pol_xx pol_xy pol_xz pol_yx pol_yy pol_yz pol_zx pol_zy pol_zz\n");
      fprintf(file_pol_ind,"#Record_number Mol_number a_ind_xx xy xz yx yy yz zx zy zz at the order N\n");
      if((*sys).typ_dip + (*sys).typ_pol){ /* Atomic dipole moment or polarizability */
	FOPEN_SAFE(file_Tij_at ,DEBUG_TIJ_AT ,"w+");
	fprintf(file_Tij_at,"#Record_number Mol_i Mol_j Tij_xx Tij_yx Tij_yy Tij_zx Tij_zy Tij_zz\n");
      }
#endif
    }
    if(step>=stepf){end_read=1;}
    fgets(chain,CHAIN_SIZE,filei);/* Read the end of the line */
      
    /*Body read if we are after the starting step*/
    if(start_read){
#ifndef YUKI
      fprintf(file_traj,"%d\n",(*sys).nb_mol);
      fprintf(file_traj,"i = %d , time = %f . Only the mass center position is written.\n",step,time);
#endif
      
      for(i=0;i<(*sys).nb_mol;i++){
	for(j=0;j<(*sys).at_p_mol;j++){
	  fscanf(filei,"%s %lf %lf %lf",chain,				\
		 (double *) &(mol[i])+(*sys).order[j],			\
		 (double *) &(mol[i])+(*sys).order[j]+1,		\
		 (double *) &(mol[i])+(*sys).order[j]+2);

	}
	pbc_mol((*sys).cell,&(mol[i])); /*Mass center + centering*/
	
#ifdef DEBUG
	fprintf(file_debug,"%10d %10d %5d 1 %f %f %f\n",step,(*sys).step_max,i,mol[i].xO ,mol[i].yO ,mol[i].zO );
	fprintf(file_debug,"%10d %10d %5d 2 %f %f %f\n",step,(*sys).step_max,i,mol[i].xH1,mol[i].yH1,mol[i].zH1);
	fprintf(file_debug,"%10d %10d %5d 3 %f %f %f\n",step,(*sys).step_max,i,mol[i].xH2,mol[i].yH2,mol[i].zH2);
#endif
      }


      /*=====================================================
	 Storage of the coordinates into an array of double
	 This array will have (sooner or later to replace the
	 unneficient structure "mol")
	 ==================================================*/
      for(i=0;i<(*sys).nb_mol;i++){
	at_coord[9*i+0]=mol[i].xO ; at_coord[9*i+1]=mol[i].yO ; at_coord[9*i+2]=mol[i].zO ; 
	at_coord[9*i+3]=mol[i].xH1; at_coord[9*i+4]=mol[i].yH1; at_coord[9*i+5]=mol[i].zH1;
	at_coord[9*i+6]=mol[i].xH2; at_coord[9*i+7]=mol[i].yH2; at_coord[9*i+8]=mol[i].zH2;
      }

      
      /*================================================================================
	Calculation of the complete dipole moment / polarizability (permanent + induced)
	================================================================================*/
      val_init(sys,mol,dip0,dip0_mol,pol0,pol0_mol,file_traj,file_dip0,file_pol0,step);/*dip0, pol0*/
      comp_dip_pol(sys,mol,at_coord,dip0,dip0_mol,dip,dip_mol,pol0,pol0_mol,pol,pol_mol,e_ind,a_ind,v3_tmp,m3_tmp,Tij_mc,Tij_at,x,y,z,r,fthole,fileo,file_Tij_mc,file_Tij_at,file_dip,file_pol,file_dip_ind,file_pol_ind,step);/*Induction*/
      

      
      (*sys).step_max++; /* Real number of steps recorded (can be lower than stepf-stepi+1) */
    }
    else{/*Passing the lines before stepi*/
      for(i=1;i<=(*sys).nb_at;i++){
	fgets(chain,CHAIN_SIZE,filei);
      }
    }
  }

  /* Release the memory */
  free(mol);
  free(at_coord);
  free(dip0);
  free(dip);
  free(dip0_mol);
  free(dip_mol);
  free(pol0);
  free(pol);
  free(pol0_mol);
  free(pol_mol);
  free(e_ind);
  free(v3_tmp);
  free(m3_tmp);
  free(Tij_mc);
  free(x);
  free(y);
  free(z);
  free(r);
  free(fthole);
  if((*sys).typ_dip + (*sys).typ_pol){free(Tij_at);}
    
  fclose(filei);
  fclose(fileo);
#ifndef YUKI
  fclose(file_traj);
#endif
#ifdef DEBUG
  fclose(file_debug);
  fclose(file_dip0);
  fclose(file_pol0);
  fclose(file_Tij_mc);
  fclose(file_dip);
  fclose(file_pol);
  fclose(file_dip_ind);
  if((*sys).typ_dip + (*sys).typ_pol){fclose(file_Tij_at);}
#endif

  

  return;
}



