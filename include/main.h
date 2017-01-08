/*Contains all:
    -the constants
    -Debug mode
    -the structures
    -the macros
*/


/*============================================*/
/*============================================*/
#ifndef CONSTANTS
#define CHAIN_SIZE 200 /* Maximal lenght of a chain of characters */

/* Molar weight */
#define MO 16. /*Oxygen*/
#define MH 1. /*Hydrogen*/

/*Convergence reached*/
#define DIP_CONV 0.001
#define POL_CONV 0.001

/*Limits to avoid Infinite SCF*/
#define A_THOLE 2.4380/0.529177249 /*Screening length (a) in Ang-1 from van Duijnen and Swart, JPCA, 102, 14, 2399, 1998 (Table 7, ab initio, expon)*/
#define POL_LIM 100
#define CONV_PB "convergence_pb.xyz"/*File where the problematic positions are written*/

#define DEBUG_TRAJ "traj.dat" /*File with trajectories info*/
#define DEBUG_DIP0 "dip0.dat" /*File with permanent dipole*/
#define DEBUG_POL0 "pol0.dat" /*File with permanent polarizability*/
#define DEBUG_DIP_IND "dip_ind.dat" /*File with the induced dipole*/
#define DEBUG_POL_IND "pol_ind.dat" /*File with the induced polarizability*/
#define DEBUG_DIP "dip.dat" /*File with complete dipole*/
#define DEBUG_POL "pol.dat" /*File with effective polarizability*/
#define DEBUG_TIJ_MC "tij.dat" /*File which contains the different Tij matrixes*/
#define DEBUG_TIJ_AT "tij_at.dat" /*File which contains the different Tij matrixes*/

#endif


/*============================================*/
/*============================================*/
#ifndef STRUCTURES
/* Structure which contains the positions
   the three different atoms composing a
   water molecule.
*/
typedef struct mol_info
{
  /*DO NOT CHANGE THE ORDER OF O, H1, H2 and dummy*/
  /*This is important for the record of atomic positions
    and the calculation of the mass center*/
  /*DO NOT CHANGE THE ORDER OF O, H1, H2 and dummy*/
  double xO,yO,zO;/*O coordinates*/
  double xH1,yH1,zH1;/*First H*/
  double xH2,yH2,zH2;/*Second H*/
  double xdummy,ydummy,zdummy;/*Second H*/
  /*DO NOT CHANGE THE ORDER OF O, H1, H2 and dummy*/

  double x,y,z;/*mass center*/
} mol_info;

/* Structure wich contains the dipole moment
   for each water molecule*/
typedef struct vect_3d
{
  double x,y,z;
} vect_3d;

/*
  Structure wich contains the symmetric
  3D-matrixes. Only the lower part is taken
  into account
*/
typedef struct mat_sym_3d
{
  double xx,yx,yy,zx,zy,zz;
} mat_sym_3d;

/*
  Structure wich contains the 3D-matrixes.
*/
typedef struct mat_3d
{
  double xx,xy,xz,yx,yy,yz,zx,zy,zz;
} mat_3d;


/*
  Information from the input file
  Starting step to record, Final step to record
  Distance between 2 t0, Window to calculate the 
  correlation function
  Name of the input and the output files
*/
typedef struct input_info
{
  int stepi,stepf; 
  char namei[CHAIN_SIZE], nameo[CHAIN_SIZE], name_mass[CHAIN_SIZE];
} input_info;

/*
  Information about the system:
  cell size
  number of atoms, of molecules
  Number of steps really recorded (can differ from
  the expectation)
  Do we use Thole damping factor
*/
typedef struct sys_info
{
  vect_3d cell;
  int at_p_mol,nb_at,nb_mol,nb_point,nb_dip,nb_pol,step_max;
  double cutoff;
  int p,q,r;
  int *order;
  int thole, intra, typ_dip, typ_pol;
  double M_x,M_y,M_z,A_xx,A_xy,A_xz,A_yx,A_yy,A_yz,A_zx,A_zy,A_zz;
  double MO_z,AO_xx,AO_yy,AO_zz;
  double MH_z,AH_xx,AH_yy,AH_zz,AH_yx,AH_zx,AH_zy;
} sys_info;
#endif


/*============================================*/
/*=================MACROS=====================*/
/*============================================*/
#ifndef MACROS
/*Secured opening of a file 
  The program is stoped if the file cannot be opened*/
#define FOPEN_SAFE(ptr,name,mode) ptr=fopen(name,mode);	\
  if (ptr == NULL){ \
    printf("Impossible to open the file %s (line %d of file %s)\n",name,__LINE__,__FILE__); \
    printf("End of program.\n"); \
    exit(0); \
  };

/*
  Secured memory allocation
  The program is stopped if the allocation can not be done
  All the blocks are initialized
*/
#define MALLOC_SAFE(ptr,nb_el,type_el) ptr = (type_el*) malloc(nb_el*sizeof(type_el)); \
  if(ptr == NULL){ \
    printf("Impossible to allocate at line %d of file %s.\n",__LINE__,__FILE__); \
    printf("End of program.\n"); \
    exit(0); \
  }; \


#endif
