/************************************H2.c***********************************/
#include "model.h"
/*extern struct ATM_LAYER layer[MAXLAYERS];*/
extern struct ATM_LAYER *layer;
extern struct useConstituents use;

double H2(int r, int j, double f)  /*from Joanna's thesis*/
{
	double PH2, PHe, PCH4, alpha, cf, T;

    if (!use.h2) return (0.0);
    if (r==1)
    {
		printf("H2 absorption included.\n");
        return ( 0 );
	}

    PH2  = layer[j].XH2*layer[j].P;
    PHe  = layer[j].XHe*layer[j].P;
    PCH4 = layer[j].XCH4*layer[j].P;
    cf = 3.9522e-14*f*f*PH2;
    T = 273.0/layer[j].T;
    alpha =cf*(PH2*pow(T,3.12)+1.382*PHe*pow(T,2.24)+9.322*PCH4*pow(T,3.34));
    return( alpha );
}
