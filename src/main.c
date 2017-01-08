/*
  The aim of this program is to simulate SFG spectra.
  The dipole moment and the polarizability are calculated
  thanks to first order approximation

  BE CAREFUL
  To write the code I supposed that:
  1)The position file is constituted of water only (dummy atoms for TIP4P can be included)
  2)The water model is not flexible (easy calculation of the bissector)
  3)The O and H are never splitted because of pbc
*/

#include <stdio.h>
#include <stdlib.h>
#include <main.h>
#include <treat_inp.h>

int main(int argc, char *argv[]){
  input_info input={0,0,"","",""};
  sys_info sys={{0.,0.,0.},0,0,0,0,0,0,0,0.0,				\
		0,0,0,NULL,0,0,0,0,					\
		0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,	\
		0.0,0.0,0.0,0.0,					\
		0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#ifdef DEBUG
  printf("You are running the DEBUG-mode.\n");
  printf("I hope you know what you are doing.\n");
  printf("If you want the release mode: change the makefile!\n\n");
#endif

  /*The file where the convergence problems are reported is removed*/
  remove(CONV_PB);

  /* Read the input file */
  read_input(&input,&sys,argc,argv);
  
  /*
    Open/Read/Close the position file
    Each step, the polarizability and
    the dipole moment are calculated
   */
  trajectory(&sys,input,argv);
  
  return 0;
}


