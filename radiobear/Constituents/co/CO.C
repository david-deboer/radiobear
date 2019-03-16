/******************************CO.c****************************************/
/*    Returns absorption in cm^-1                                         */

#include "model.h"
/*extern struct ATM_LAYER layer[MAXLAYERS];*/
extern struct ATM_LAYER *layer;
extern struct useConstituents use;

/*     constants     */
#define coef      7.244E+21      /* coefficient from GEISA_PP.TEX eq. 14    */
#define NUM_LINES 30             /* number of CO lines used for calculation*/
#define T0        296.0          /* reference temperature in K              */
#define hck       1.438396       /* hc/k  [K cm]                            */
#define GHz       29.9792458     /* conversion from cm^-1 to GHz            */

double CO(int r, int j, double f)
{
	int i, NLIN;
    double shape, g2, f2, z2, num, den, alpha=0.0, ITG;
	double zeta=0.0, delta=0.0;
    double PH2, PHe, PCO, gamma, T, P;
    double GH2, GHe, GCO;
    double n_dvl, n_int;
    static double *f0, *I0, *E;
    FILE *ifp;

    NLIN = 26;
    if (!use.co) return(0.0);
    if (r==1)
    {
		printf("Allocating array memory for CO absorption.");
        f0   = (double *) calloc(NUM_LINES,sizeof(double));
        I0   = (double *) calloc(NUM_LINES,sizeof(double));
        E    = (double *) calloc(NUM_LINES,sizeof(double));
        printf("\rReading in  CO line data...                  ");   /* read co.lin */
        ifp = fopen("co.lin","r");
        for (i=0; i<NLIN; ++i)
              fscanf(ifp,"%f %f %f\n",f0+i,I0+i,E+i);
        printf("\rReading in  CO line data...  Done.  (%d lines)\n",NLIN);
        fclose(ifp);
        return ( 0 );
      }

      P = layer[j].P;
      T = layer[j].T;
      PH2 = layer[j].XH2*P;
      PCO = layer[j].XCO*P;
      PHe = layer[j].XHe*P;

/*******************************************************/
/* Set these parameters for the different constituents */
      GH2 = 1.960;
      GHe = 1.200;
      GCO = 6.000;
      n_dvl = 0.7;
      n_int = 3.0/2.0;
/*******************************************************/

      /* Ben-Reuven parameters */
      gamma = pow((T0/T),n_dvl)*(GH2*PH2 + GHe*PHe + GCO*PCO);
	  g2 = SQ(gamma);
      f2 = SQ(f);
	  z2 = SQ(zeta);

      for (i=0; i < NLIN; ++i)
      {
          ITG = I0[i]*exp(-((1.0/T)-(1.0/T0))*E[i]*hck);
          num = (gamma-zeta)*f2 + (gamma+zeta)*( SQ(f0[i]+delta) + g2 - z2);
          den = SQ(f2 - SQ(f0[i]+delta) - g2 + z2) + 4.0*f2*g2;
          shape = GHz*2.0*SQ(f/f0[i])*num/(PI*den);
          alpha += shape*ITG;
      }
      return( coef*(PCO/T0)*pow((T0/T),n_int+2)*alpha );
}
#undef NUM_LINES
#undef T0
