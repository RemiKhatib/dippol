/*
  In this file can be found the subroutines
  which are necessary to calculate correlations
  fucntions

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <main.h>

/* Correlation function <a(t0+t).m(t0)> */
void write_dippol(sys_info sys, vect_3d *dip0, vect_3d *dip, mat_sym_3d *pol0, mat_3d *pol, FILE *fileo){
  int i=0, inc_dip=0, inc_pol=0, inc_pol0=0;
  int v3_shift=sizeof(vect_3d)/sizeof(double), m3_shift=sizeof(mat_3d)/sizeof(double), m3_sym_shift=sizeof(mat_sym_3d)/sizeof(double);
  double *p_dip=NULL, *p_pol=NULL, *p_polb=NULL, *p_pol0=NULL, *p_dip0=NULL, *p_pol0b=NULL;
  
  /* Selection of the polarization p_dip and 
     p_pol which will point exactly on the good 
     element of the structure*/

  /* Normal calculation method (XXZ, YYZ, XZX, ...)*/
  if(sys.p>=0){
    
    if(sys.p==0&&sys.q==0){p_pol0=&(pol0[0].xx);p_pol=&(pol[0].xx);}
    else if(sys.p==0&&sys.q==1){p_pol0=&(pol0[0].yx);p_pol=&(pol[0].xy);}
    else if(sys.p==0&&sys.q==2){p_pol0=&(pol0[0].zx);p_pol=&(pol[0].xz);}
    else if(sys.p==1&&sys.q==0){p_pol0=&(pol0[0].yx);p_pol=&(pol[0].yx);}
    else if(sys.p==1&&sys.q==1){p_pol0=&(pol0[0].yy);p_pol=&(pol[0].yy);}
    else if(sys.p==1&&sys.q==2){p_pol0=&(pol0[0].zy);p_pol=&(pol[0].yz);}
    else if(sys.p==2&&sys.q==0){p_pol0=&(pol0[0].zx);p_pol=&(pol[0].zx);}
    else if(sys.p==2&&sys.q==1){p_pol0=&(pol0[0].zy);p_pol=&(pol[0].zy);}
    else{p_pol0=&(pol0[0].zz);p_pol=&(pol[0].zz);}

    switch(sys.r){
    case 0 :
      p_dip=&(dip[0].x);
      p_dip0=&(dip0[0].x);
      break;
    case 1 :
      p_dip=&(dip[0].z);
      p_dip0=&(dip0[0].y);
      break;
    default :
      p_dip=&(dip[0].z);
      p_dip0=&(dip0[0].z);
      break;
    }

    for(i=0;i<sys.nb_mol;i++){
      inc_dip =v3_shift*i;
      inc_pol =m3_shift*i;
      inc_pol0=m3_sym_shift*i;
      
#ifndef YUKI
      fprintf(fileo,"%12.7f %12.7f %12.7f %12.7f\n",p_dip0[inc_dip],p_dip[inc_dip], \
	      p_pol0[inc_pol0],p_pol[inc_pol]);
#else
      /*Conversion from Ang^3 to bohr^3
       Then conversion from real to integers with Yuki style */
      fprintf(fileo,"%5.0f%5.0f%5.0f\n",\
	      (p_dip0[inc_dip])*1000+5000,				\
	      (p_dip[inc_dip]-p_dip0[inc_dip])*1000+5000, \
	      (p_pol[inc_pol]/pow(0.529177249,3)*10000-20000)*0.3				\
	      );
#endif

    }
  }
  else{  /* An average will be done between 2 polarization (eg XXZ and YYZ) */
    if(sys.p==-1){/*YYX+ZZX*/
      p_pol=&(pol[0].yy);
      p_polb=&(pol[0].zz);
      p_pol0=&(pol0[0].yy);
      p_pol0b=&(pol0[0].zz);
      p_dip=&(dip[0].x);
      p_dip0=&(dip0[0].x);
    }
    else if(sys.p==-2){/*XXY+ZZY*/
      p_pol=&(pol[0].xx);
      p_polb=&(pol[0].zz);
      p_pol0=&(pol0[0].xx);
      p_pol0b=&(pol0[0].zz);
      p_dip=&(dip[0].y);
      p_dip0=&(dip0[0].y);
    }
    else{/*XXZ+YYZ*/
      p_pol=&(pol[0].xx);
      p_polb=&(pol[0].yy);
      p_pol0=&(pol0[0].xx);
      p_pol0b=&(pol0[0].yy);
      p_dip=&(dip[0].z);
      p_dip0=&(dip0[0].z);
    }

    for(i=0;i<sys.nb_mol;i++){
      inc_dip =v3_shift*i;
      inc_pol =m3_shift*i;
      inc_pol0=m3_sym_shift*i;
      
#ifndef YUKI
      fprintf(fileo,"%12.7f %12.7f %12.7f %12.7f\n",p_dip0[inc_dip],p_dip[inc_dip],(p_pol0[inc_pol0]+p_pol0b[inc_pol0])/2, \
	      (p_pol[inc_pol]+p_polb[inc_pol])/2);
#else
      /*Conversion from Ang^3 to bohr^3
       Then conversion from real to integers with Yuki style */
      fprintf(fileo,"%5.0f%5.0f%5.0f\n",				\
	      (p_dip0[inc_dip])*1000+5000,				\
	      (p_dip[inc_dip]-p_dip0[inc_dip])*1000+5000,		\
	      (((p_pol[inc_pol]+p_polb[inc_pol])/2/pow(0.529177249,3))*10000-20000)*0.3	\
	      );
#endif
    }

  }
  return;
}
