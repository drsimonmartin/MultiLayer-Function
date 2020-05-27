// multiSlab neutron reflectivity simulator

// 30 Jan 2006 remove GSL error function - use proFit version
#include "proFit_interface.h"
#ifndef __MACH__
#include <fp.h>
#endif
#include <math.h>
#include <complex.h>
// #include <gsl/gsl_sf_erf.h>

#define SLABS 8
// have six layers + top and bottoms
_Complex double g; /* a test var */
/* some globals */
double sld_array[SLABS]; /*array of sld */
double thick_array[SLABS]; /*array of thickness */
double rough_array[SLABS];/*array of roughness */
int nlayers; /*total layers in system */
int nfunclayers; /* no of layers in functional part of model */


/***************************************************************************************/

void SetUp (	short* const moduleKind,		/* set moduleKind to isFunction or isProgram */
Str255 name,					/* the name of the program or function (pascal string) */
SInt32* const requiredGlobals,	/* the number of bytes to be allocated in ExtModulesParamBlock.globals */
/* set requiredGlobals to 0 if you don't use this feature */
ExtModulesParamBlock* pb)		/* the complete parameter block passed by pro Fit to the */
/* routines defined in this file. In most cases it can be ignored */
/* SetUp is called once when the external module is linked to proFit */
{
	*moduleKind=isFunction;					/* we define a function */
	SetPascalStr(name,"\pSlab_fit",255);		/* with the name "Slab_fit" */
				 *requiredGlobals=0;					/* we define no globals */
}

/***************************************************************************************/

void InitializeFunc (
					 Boolean* const hasDerivatives,	/* set this to true if and only if you define the function */
					 /* Derivatives to calculate the partial derivatives of the parameters */
					 Str255 descr1stLine,			/* first line of the text in the parameter window */
					 Str255 descr2ndLine,			/* second line of the text in the parameter window */
					 short* const numberOfParams,	/* the number of parameters of the function */
					 DefaultParamInfo* const a0,		/* the default names, values etc. of the parameters */
					 ExtModulesParamBlock* pb)		/* the complete parameter block passed by pro Fit to the */
					 /* routines defined in this file. In most cases it can be ignored */
					 /* InitializeFunc is called once (after SetUp has been called) when the external module is linked to proFit */
					 /* Used to set all the information needed to describe a function */
					 {	
	*hasDerivatives=false;
	SetPascalStr(descr1stLine,"\pA test",255);
	SetPascalStr(descr2ndLine,"\p",255);
	
	*numberOfParams=31;							/* we have 30 parameters */
	
	/* The following is to set parameter names, fitting modes, etc. */
	
	(*a0->value)[0] = 0.0;						/* set their names and defaults	*/
	(*a0->lowest)[0]= 0.;
	(*a0->mode)[0] = inactive;
	SetPascalStr((*a0->name)[0],"\pTopSLD",  maxParamNameLength);
	
	(*a0->value)[1]= 0.0;
	(*a0->lowest)[1]= 0.;
	(*a0->mode)[1] = inactive;
	SetPascalStr((*a0->name)[1],"\psld1", maxParamNameLength);
	
	(*a0->value)[2]= 0.0;
	(*a0->lowest)[2]= 0.;
	(*a0->mode)[2] = inactive;
	SetPascalStr((*a0->name)[2],"\prough1", maxParamNameLength);
	
	(*a0->value)[3]= 0.0;
	(*a0->lowest)[3]= 0.0;
	(*a0->mode)[3] = inactive;
	SetPascalStr((*a0->name)[3],"\pTK_1", maxParamNameLength);
	
	(*a0->value)[4]= 0.0; /* Sld of layer 2 */
	(*a0->lowest)[4]= 0.;
	(*a0->mode)[4] = inactive;
	SetPascalStr((*a0->name)[4],"\psld2", maxParamNameLength);
	
	(*a0->value)[5]= 0.0; /* roughness of layer 2*/
	(*a0->lowest)[5]= 0.;
	(*a0->mode)[5] = inactive;
	SetPascalStr((*a0->name)[5],"\prough2", maxParamNameLength);
	
	(*a0->value)[6]= 0.0; /* thickness of layer 2*/
	(*a0->lowest)[6]= 0.;
	(*a0->mode)[6] = inactive;
	SetPascalStr((*a0->name)[6],"\pTK_2", maxParamNameLength);
	
	(*a0->value)[7]= 0.0; /* SLd of component of layer 3*/
	(*a0->lowest)[7]= 0.;
	(*a0->mode)[7] = inactive;
	SetPascalStr((*a0->name)[7],"\psld_3", maxParamNameLength);
	
	(*a0->value)[8]= 1.0; /* roughness 3 */
	(*a0->lowest)[8]= 0.;
	//(*a0->highest)[8] = 1.;
	(*a0->mode)[8] = inactive;
	SetPascalStr((*a0->name)[8],"\prough3", maxParamNameLength);
	
	(*a0->value)[9]= 0.0; /* thickness of layer 3*/
	(*a0->lowest)[9]= 0.;
	(*a0->mode)[9] = inactive;
	SetPascalStr((*a0->name)[9],"\pTK_3", maxParamNameLength);
	
	(*a0->value)[10]= 0.0; /* SLD 4 */
	(*a0->lowest)[10]= 0.;
	(*a0->mode)[10] = inactive;
	SetPascalStr((*a0->name)[10],"\pSLD 4", maxParamNameLength);
	
	(*a0->value)[11]= 0.0; /* rough 4 */
	(*a0->mode)[11] = inactive;
	(*a0->lowest)[11]= 0.;
	SetPascalStr((*a0->name)[11],"\prough4", maxParamNameLength);
	
	(*a0->value)[12]= 0.0; /* thick 4*/
	(*a0->lowest)[12]= 0.;
	(*a0->mode)[12] = inactive;
	SetPascalStr((*a0->name)[12],"\pTK4", maxParamNameLength);
	
	(*a0->value)[13]= 0.0; /* SLD 5 */
	(*a0->lowest)[13]= 0.;
	(*a0->mode)[13] = inactive;
	SetPascalStr((*a0->name)[13],"\pSLD5", maxParamNameLength);
	
	(*a0->value)[14]= 0.0; /* roughness 5 */
	(*a0->lowest)[14]= 0.;
	(*a0->mode)[14] = inactive;
	SetPascalStr((*a0->name)[14],"\prough5", maxParamNameLength);
	
	(*a0->value)[15]= 0.0; /* Thickness 5 */
	(*a0->lowest)[15]= 0.;
	(*a0->mode)[15] = inactive;
	SetPascalStr((*a0->name)[15],"\pTK 5", maxParamNameLength);
	
	(*a0->value)[16]= 0.; /* SLD 6 */
	(*a0->lowest)[16]= 0.;
	(*a0->mode)[16] = inactive;
	SetPascalStr((*a0->name)[16],"\pSLD 6", maxParamNameLength);
	
	(*a0->value)[17]= 0.0; /*roughness 6*/
	(*a0->lowest)[17]= 0.;
	(*a0->mode)[17] = inactive;
	SetPascalStr((*a0->name)[17],"\prough 6", maxParamNameLength);
	
	(*a0->value)[18]= 0.0; /* Thick 6 */
	(*a0->lowest)[18]= 0.;
	(*a0->mode)[18] = inactive;
	SetPascalStr((*a0->name)[18],"\pTK 6", maxParamNameLength);
	
	(*a0->value)[19]= 2.07e-6; /* SLD substrate */
	(*a0->lowest)[19]= 0.;
	(*a0->highest)[19] = 1.;
	(*a0->mode)[19] = inactive;
	SetPascalStr((*a0->name)[19],"\pSLDsubs", maxParamNameLength);
	
	(*a0->value)[20]= 4.0; /* rough subs */
	(*a0->lowest)[20]= 0.;
	(*a0->mode)[20] = inactive;
	SetPascalStr((*a0->name)[20],"\pRough subs", maxParamNameLength);
	
	(*a0->value)[21]= 5e-7; /* background */
	(*a0->lowest)[21]= 0.;
	(*a0->mode)[21] = inactive;
	SetPascalStr((*a0->name)[21],"\pBground", maxParamNameLength);
	
	(*a0->value)[22]= 1.0; /* scalefactor */
	(*a0->lowest)[22]= 0.;
	(*a0->mode)[22] = inactive;
	SetPascalStr((*a0->name)[22],"\pScaleF", maxParamNameLength);
	
	(*a0->value)[23]= 0.04; /* resolution */
	(*a0->lowest)[23]= 0.;
	(*a0->mode)[23] = inactive;
	SetPascalStr((*a0->name)[23],"\pResn", maxParamNameLength);
	
	(*a0->value)[24]= 4.4; /* Wavelength */
	(*a0->lowest)[24]= 0.;
	 (*a0->mode)[24] = inactive;
	 SetPascalStr((*a0->name)[24],"\plambda", maxParamNameLength);
						 
    (*a0->value)[25]= 0; /* number of multilayers */
    (*a0->lowest)[25]= 0.;
    (*a0->mode)[25] = inactive;
    SetPascalStr((*a0->name)[25],"\prepeats", maxParamNameLength);

    (*a0->value)[26]= 0; /* Where multilayers are positioned */
    (*a0->lowest)[26]= 0.;
    (*a0->mode)[26] = inactive;
    SetPascalStr((*a0->name)[26],"\pMlStart", maxParamNameLength);

    (*a0->value)[27]= 10; /* thickness layer 1 */
    (*a0->lowest)[27]= 0.1;
     (*a0->mode)[27] = inactive;
     SetPascalStr((*a0->name)[27],"\pt1", maxParamNameLength);

    (*a0->value)[28]= 1e-6; /* sld layer 1 */
    (*a0->lowest)[28]= 0.1e-6;
    (*a0->mode)[28] = inactive;
    SetPascalStr((*a0->name)[28],"\pSLDl1", maxParamNameLength);
                         
    (*a0->value)[29]= 10.; /* thickness layer 2 */
    (*a0->lowest)[29]= 0.1;
    (*a0->mode)[29] = inactive;
    SetPascalStr((*a0->name)[29],"\pt2", maxParamNameLength);

    (*a0->value)[30]= 1e-6; /* sld layer 2 */
    (*a0->lowest)[30]= 0.1e-6;
    (*a0->mode)[30] = inactive;
    SetPascalStr((*a0->name)[30],"\pSLDl2", maxParamNameLength);
                                              

	
	
	pb->numberOfAdditionalReturnVals = 5;
	(*pb->y0.value)[0] = 0;
	SetPascalStr((*pb->y0.name)[0],"\pR(Q)",maxParamNameLength);
	(*pb->y0.value)[1] = 1;
	SetPascalStr((*pb->y0.name)[1],"\pSLD profile",maxParamNameLength);
	(*pb->y0.value)[2] = 0.;
	SetPascalStr((*pb->y0.name)[2],"\pRQ^4",maxParamNameLength);
	(*pb->y0.value)[3] = 0.;
	SetPascalStr((*pb->y0.name)[3],"\pphi[z]",maxParamNameLength); // irrelevant here
	(*pb->y0.value)[4] = 0.;
	SetPascalStr((*pb->y0.name)[4],"\pthickness",maxParamNameLength); // the total thickness of the sample
	(*pb->y0.value)[5] = 0.;
	SetPascalStr((*pb->y0.name)[5],"\pRTheta",maxParamNameLength); // the data is a function of theta
						 
}
					 
/***************************************************************************************/
 
short Check(short paramNo,						/* the parameter that was changed */
DefaultParamInfo* const a0,		/* the default names, values etc of the paramters */
ExtModulesParamBlock* pb)		/* the complete parameter block passed by pro Fit to the */
/* routines defined in this file. In most cases it can be ignored */
/* Can be left emtpy (returning good) if not needed. */
/* called when the user has changed a value in the parameters window. This routine */
/* can then check if this parameters is fine. It can also change some of the */
/* other entries in a0. The returned values can be: */
/*	good:		return this value if you agree with the new parameter value */
/*	update:		return this value if you want the parameters window */
/*				to be updated because you changed some of the values in a0 */
/*	bad:		return this value if you want the new parameter value to be refused */
{
	return good;
}

/***************************************************************************************/

void First (	ParamArray a,				/* the new parameters */
ExtModulesParamBlock* pb)	/* the complete parameter block passed by pro Fit to the */
/* routines defined in this file. In most cases it can be ignored */
/* Can be left emtpy if not needed. */
/* Called whenever the parameters are changed. Can be used to accelerate */
/* some calculations. See manual for more info */

{
	//#ifdef __leavein__
	int i=0.;
	/* set up rho with a sld 'profile' of the system */
	sld_array[0]=a[0];
	rough_array[0]=0.;
	thick_array[0]=0.;
    int total_slabs=SLABS+2*(int)a[25]; /* account for any mulitlayers in layer count*/
    int j=1; // counting variable used to keep track of slabs
    
    /* check to see if need to insert a multilayer section */
    if (total_slabs>SLABS)
    {
        // have to insert multilayer at correct point
        for(i=1;i<SLABS-1;i++)
        {
            if (i!=(int)a[26])
            {
                // have ordinary layer
                sld_array[j]=a[3*i-2];
                rough_array[j]=a[3*i - 1];
                thick_array[j]=a[3*i];
                j++;
            } else // have to insert multilayer structure
            {
                /* now insert multilayer */
                for (int k=1;k<=2*(int)a[25];k++)
                {
                    // check if k is even/odd, insert appropriate thickness, sld.
                    if(k % 2 == 0)
                    {
                        // is even
                        sld_array[j]  =a[30];
                        rough_array[j]=a[3*i - 1];
                        thick_array[j]=a[29];
                        j++;
                    }else
                    {
                        sld_array[j]  =a[28];
                        rough_array[j]=a[3*i - 1];
                        thick_array[j]=a[27];
                        j++;
                    }
                    
                }
            }
            
        }
        //WriteNumber(j);
    }else
    {
        // have ordinary system
        for(i=1;i<SLABS-1;i++)
        {
            sld_array[j]=a[3*i-2];
            rough_array[j]=a[3*i - 1];
            thick_array[j]=a[3*i];
            j++;
        }
        
    }
    
    /* now set up substrate values */
	sld_array[j]=a[19];
	rough_array[j]=a[20];
	thick_array[j]=0.;
	// #endif
    
    /* for (int l=0;l<=SLABS+2*(int)a[25];l++)
    {
        WriteNumber(sld_array[l]); Write(" "); WriteNumber(thick_array[l]); Write(" ") ;WriteNumber(rough_array[l]);Writeln(" ");
    }
    Writeln("");*/
}
/***************************************************************************************/



/***************************************************************************************/
/* some functions for the reflectivity calcultation */
_Complex double fresnel(_Complex double q1,_Complex double q2,double rough)
{
	return((q1-q2)/(q1+q2))*cexp(-0.5*q1*q2*rough*rough);
}
_Complex double  qmedium(_Complex double q0,_Complex double sldmedium)
{
	
	return csqrt(q0*q0-(16.*M_PI + 0.I)*sldmedium);
}

double CReflectivity( _Complex double Q,int nlayers)  // calcs reflectivity
{
	
	int i;
	_Complex double q0;
	_Complex double q1;
	_Complex double q2;
	_Complex double R;
	_Complex double numer=0;
	_Complex double denom=0;
	_Complex double rn;
	
	// start by calculating q in capping medium
	  q0 = cabs(csqrt((Q)*(Q)+(16.*M_PI + 0.I)*(sld_array[0])));
	R=0.; // make sure start from correct place
	
	// calc first term seperately
	// saves some floating point maths
	
//	q1=qmedium(q0,sld_array[nlayers]); // + 0.I);
//	q2=qmedium(q0,sld_array[nlayers+1] ); //+ 0.I);
//	R=fresnel(q1,q2,rough_array[nlayers+1]);
	
	
	
//	for (i=nlayers-1;i>=0;i--) // loop through the layers
for (i=nlayers;i>=0;i--) // loop through the layers
	{
		
		q1=qmedium(q0,sld_array[i]); // + 0.I);
		q2=qmedium(q0,sld_array[i+1]); //+ 0.I);
		rn=fresnel(q1,q2,rough_array[i+1]);
		
		numer=rn + R*cexp(1.I*q2*thick_array[i+1]);
		denom=1. + (rn*R*cexp(1.I*q2*thick_array[i+1]));
		
		R=numer/denom;
		
		}
	
	//R=numer/denom;
	
	return creal(R*conj(R));	
	
}

// Function to return value of scattering length density at z
double calcProfile(double *z,int nSlabs)
{
	double profSum;   //{sum variable used to calculate the profile}
	int i;            // {a counting variable}
	double zi;		  // distance from top of sample to ith layer}
	double sigma;     // {roughness of ith layer}
	
	profSum = sld_array[0];   //{ initialize the sum variable}
	zi=0;             //{start from z=0}
	for (i=1;i<=nSlabs-1;i++)
	{
		if (rough_array[i]==0.)  
		{ 
			sigma=1.e-3;
		}else sigma=rough_array[i];
		// WriteNumber(sld_array[i]);
		//profSum = profSum + (sld_array[i] - sld_array[i-1])*(2. - gsl_sf_erfc((*z-zi)*sqrt(2.)/(2.*sigma))/2.);
		
		profSum=profSum + (sld_array[i] - sld_array[i-1])*(1. + Erf((*z-zi)/sqrt(2.)/sigma))/2.;
		zi = zi + thick_array[i];
			}
	
	return(profSum);
}


/***************************************************************/

void Func (		double x,						/* the x-value */
ParamArray a,					/* the parameters */
double* const y,				/* the y-value to be returned */
ExtModulesParamBlock* pb)		/* the complete parameter block passed by pro Fit to the */
/* routines defined in this file. In most cases it can be ignored */

{
	double sum=0.;
	double backg; /* background*/
	double scalef; /* scaling factor */
    int total_slabs=SLABS+2*(int)a[25];
	backg=a[21];
	scalef=a[22];
	sum=0.;
    			/************************************************************************************************/
				nlayers=total_slabs-2;
    Write("Nlayers=");WriteNumber(nlayers);
if (Output(0)) // ordinary calculation of reflectivity
{
	sum=sum+ 0.135*CReflectivity(x-x*a[23],nlayers);
	sum=sum+ 0.135*CReflectivity(x+x*a[23],nlayers);  
	sum=sum+ 0.198*CReflectivity(x-x*a[23]*0.9,nlayers); 
	sum=sum+ 0.198*CReflectivity(x+x*a[23]*0.9,nlayers); 
	sum=sum+ 0.278*CReflectivity(x-x*a[23]*0.8,nlayers); 
	sum=sum+ 0.278*CReflectivity(x+x*a[23]*0.8,nlayers); 
	sum=sum+ 0.375*CReflectivity(x-x*a[23]*.7,nlayers); 
	sum=sum+ 0.375*CReflectivity(x+x*a[23]*.7,nlayers); 
	sum=sum+ 0.487*CReflectivity(x-x*a[23]*.6,nlayers); 
	sum=sum+ 0.487*CReflectivity(x+x*a[23]*.6,nlayers); 
	sum=sum+ 0.606*CReflectivity(x-x*a[23]*.5,nlayers); 
	sum=sum+ 0.606*CReflectivity(x+x*a[23]*.5,nlayers); 
	sum=sum+ 0.726*CReflectivity(x-x*a[23]*.4,nlayers); 
	sum=sum+ 0.726*CReflectivity(x+x*a[23]*.4,nlayers); 
	sum=sum+ 0.835*CReflectivity(x-x*a[23]*.3,nlayers); 
	sum=sum+ 0.835*CReflectivity(x+x*a[23]*.3,nlayers); 
	sum=sum+ 0.923*CReflectivity(x-x*a[23]*.2,nlayers); 
	sum=sum+ 0.923*CReflectivity(x+x*a[23]*.2,nlayers); 
	sum=sum+ 0.98*CReflectivity(x-x*a[23]*.1,nlayers); 
	sum=sum+ 0.98*CReflectivity(x+x*a[23]*.1,nlayers); 
	sum=sum + CReflectivity(x,nlayers);
	
	
	//*y=(backg + sum/12.086);
	y[0]=(backg + scalef*sum/12.086);
}
if(Output(1))
{
	// profile has been requested
	y[1]=calcProfile(&x,total_slabs);
}
if (Output(2)) // Calculate RQ^4
{
	sum=sum+ 0.135*CReflectivity(x-x*a[23],nlayers);
	sum=sum+ 0.135*CReflectivity(x+x*a[23],nlayers);  
	sum=sum+ 0.198*CReflectivity(x-x*a[23]*0.9,nlayers); 
	sum=sum+ 0.198*CReflectivity(x+x*a[23]*0.9,nlayers); 
	sum=sum+ 0.278*CReflectivity(x-x*a[23]*0.8,nlayers); 
	sum=sum+ 0.278*CReflectivity(x+x*a[23]*0.8,nlayers); 
	sum=sum+ 0.375*CReflectivity(x-x*a[23]*.7,nlayers); 
	sum=sum+ 0.375*CReflectivity(x+x*a[23]*.7,nlayers); 
	sum=sum+ 0.487*CReflectivity(x-x*a[23]*.6,nlayers); 
	sum=sum+ 0.487*CReflectivity(x+x*a[23]*.6,nlayers); 
	sum=sum+ 0.606*CReflectivity(x-x*a[23]*.5,nlayers); 
	sum=sum+ 0.606*CReflectivity(x+x*a[23]*.5,nlayers); 
	sum=sum+ 0.726*CReflectivity(x-x*a[23]*.4,nlayers); 
	sum=sum+ 0.726*CReflectivity(x+x*a[23]*.4,nlayers); 
	sum=sum+ 0.835*CReflectivity(x-x*a[23]*.3,nlayers); 
	sum=sum+ 0.835*CReflectivity(x+x*a[23]*.3,nlayers); 
	sum=sum+ 0.923*CReflectivity(x-x*a[23]*.2,nlayers); 
	sum=sum+ 0.923*CReflectivity(x+x*a[23]*.2,nlayers); 
	sum=sum+ 0.98*CReflectivity(x-x*a[23]*.1,nlayers); 
	sum=sum+ 0.98*CReflectivity(x+x*a[23]*.1,nlayers); 
	sum=sum + CReflectivity(x,nlayers);
	
	
	//*y=(backg + sum/12.086);
	y[2]=x*x*x*x*(backg + scalef*sum/12.086);
}

if (Output(4))
{
	// return total thickness of sample - useful for SLD plots
	y[4]=(a[3] + a[6] + a[9] + a[12] + a[15] + a[18]);
}

//*y=gsl_sf_erfc(-1000.);	
}

/***************************************************************************************/

void Derivatives(double x,						/* the x-value */
ParamArray a,					/* the parameters */
ParamArray dyda,				/* the derivatives to be returned */
ExtModulesParamBlock* pb)		/* the complete parameter block passed by pro Fit to the */
/* routines defined in this file. In most cases it can be ignored */
/* Can be left empty if InitializeFunc sets hasDerivatives to false */
/* called to calculate the partial derivatives of the function with respect to */
/* its parameters. If you leave this function empty and set hasDerivatives to false in */
/* FuncInitialize, the derivatives will be calcuated numerically, otherwise pro Fit */
/* calls this function to obtain the values of ALL derivatives. */
/* As a result of the numerical calculation fitting will be slower */
{
}



/***************************************************************************************/

void Last (ExtModulesParamBlock* pb)
/* Can be left emtpy if not needed. */
/* Called when calculating is through. See manual for more info */
{
}

/***************************************************************************************/

void CleanUp (ExtModulesParamBlock* pb)
/* called when the external module is removed from pro Fit's menus */
/* in most cases, this function can be empty */
{
}





/***************************************************************************************/
/* for programs, not used here: */
/***************************************************************************************/

void InitializeProg (ExtModulesParamBlock* pb)
{}


void Run(ExtModulesParamBlock* pb)
{
}

