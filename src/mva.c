#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>




SEXP mva(SEXP R_x, SEXP R_mvg) {
     int i;
     int xlen = length(R_x); 
     double *x0 = REAL(R_x);
     int mvg = INTEGER(R_mvg)[0]; 
     
     /* result */
     SEXP res;
     PROTECT(res = allocVector(REALSXP, xlen));
     double *x1 = REAL(res); 

     double *sx = (double *) R_alloc(xlen, sizeof(double));
     double temp=0;
     
     /* first compute cumsum */
     for(i=0;i<xlen;i++){
	  temp += x0[i];
	  sx[i]=temp;
     }
	
     /* compute moving average. Be careful at the boundaries. */
     for(i=0;i<xlen;i++){
	  if(i<(mvg-1)){
	       x1[i] = sx[mvg-2]/(mvg-1);
	  } else if(i==(mvg-1)){
	       x1[i] = sx[i+mvg-1]/(2*mvg-1);
	  } else if(i>(xlen-mvg)){
	       x1[i] = (sx[xlen-1]-sx[xlen-mvg])/(mvg-1);
	  } else if(i==(xlen-mvg)){
	       x1[i] = (sx[xlen-1]-sx[xlen-2*mvg])/(2*mvg-1);
	  } else{
	       x1[i] = (sx[i+mvg-1]-sx[i-mvg])/(2*mvg-1);
	  }
     }

     UNPROTECT(1);
     return(res);
}


static R_CallMethodDef callMethods[]  = {
  {"mva", (DL_FUNC) &mva, 2},
  {NULL, NULL,0}
};

void R_init_ChIPComp( DllInfo *info )
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}








