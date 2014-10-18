/*
interp1_SAscan - 1-D interpolation routine for beamforming a synthetic aperture scan (creates subimages)

Loops through 3-D rf data using different delays for each transmit event and receive channel
 
CAUTION: NO TYPE CHECKING IS PERFORMED
All inputs are "single" type variables except mask (logical)
			   
Vout = interp1_DAS(x1,dx,v,tt,tr,mask,nan,spline,channel);

    x1 - Input vector start
    dx - Input vector spacing
    v - Input values
    tt - Transmit distances
    tr - Receive distances
    mask - Masking function as a function of transmit event
    nan - NaN value
    spline - 1 for spline, 0 for linear
    channel - 1 for channel data, 0 for summed data
*/

#include "mex.h"
#include <omp.h>

/* Input arguments */
#define I_X1        prhs[0]
#define I_DX	    prhs[1]
#define I_V         prhs[2]
#define I_XQ1       prhs[3]
#define I_XQ2       prhs[4]
#define I_MASK      prhs[5]
#define I_NAN       prhs[6]
#define I_SPLINE    prhs[7]
#define I_CHANNEL   prhs[8]

/* Output arguments */
#define O_VOUT        plhs[0]

/*Linear interpolation routine*/
static void linterp(float *vout, float x1, float dx, float *v, float *tt, float *tr, bool *mask, int ax, int tx, int rx, int ninterp, float nan, int channel)
{
    int i,j,k, whole;
    float cur;
    float *dv;
    dv = (float*) mxMalloc(sizeof(float)*ax);
	
    for (k=0;k<tx;k++)
    {
        for (j=0;j<rx;j++)
        {
            /*Calculate slopes*/
            #pragma omp parallel for private(i) shared(ax,rx,j,k,v,dv)
            for (i=0;i<ax-1;i++)
            {
                dv[i]=v[k*rx*ax+j*ax+i+1]-v[k*rx*ax+j*ax+i];
            }
            dv[ax-1]=0;
	
            /*Rescale interp vector to indices, calculate output value*/
            #pragma omp parallel for private(i,cur,whole) shared(ax,rx,ninterp,nan,j,k,v,dv,x1,dx,vout,mask)
            for (i=0;i<ninterp;i++)
            {
                
                if((channel==1&&(mask[k*ninterp+i]==0||vout[j*ninterp+i]==nan))||(channel==0&&(mask[k*ninterp+i]==0||vout[i]==nan))){}
                else
                {
                    cur=(tt[k*ninterp+i]+tr[j*ninterp+i]-x1)/dx;
                    if(cur > ax-1 || cur < 0)
                    { 
                        if(channel==0)
                            vout[i]=nan;
                        else
                            vout[j*ninterp+i]=nan;
                    }
                    else
                    {
                        whole=(int)cur; /*Cast instead of floor, fine for positive numbers*/
                        if(channel==0)
                            vout[i]=vout[i]+v[k*rx*ax+j*ax+whole]+dv[whole]*(cur-whole);
                        else
                            vout[j*ninterp+i]=vout[j*ninterp+i]+v[k*rx*ax+j*ax+whole]+dv[whole]*(cur-whole);
                    }
                }
            }
        }
    }
	
   /*Clean up dynamic array*/
   mxFree(dv);
}

/*Cubic spline interpolation routine*/
static void sinterp(float *vout, float x1, float dx, float *v, float *tt, float *tr, bool *mask, int ax, int tx, int rx, int ninterp, float nan, int channel)
{
    int i,j,k, whole;
	float cur, t;
	float *d, *cp, *dp, *kp, *a, *b;
	d = (float*) mxMalloc(sizeof(float)*ax);
	cp = (float*) mxMalloc(sizeof(float)*ax);
    dp = (float*) mxMalloc(sizeof(float)*ax);
    kp = (float*) mxMalloc(sizeof(float)*ax);
    a = (float*) mxMalloc(sizeof(float)*(ax-1));
    b = (float*) mxMalloc(sizeof(float)*(ax-1));
       
    /*Calculate data-independent parameter*/
    cp[0]=.5;
    for (i=1;i<ax-1;i++)
    {
        cp[i]=1/(4-cp[i-1]);
    }
    cp[ax-1]=1/(2-cp[ax-2]);
    
    for (k=0;k<tx;k++)
    {
        for (j=0;j<rx;j++)
        {
            /*Calculate (modified) slopes*/
            d[0]=v[k*rx*ax+j*ax+1]-v[k*rx*ax+j*ax+0];
            for (i=1;i<ax-1;i++)
            {
                d[i]=v[k*rx*ax+j*ax+i+1]-v[k*rx*ax+j*ax+i-1];
            }
            d[ax-1]=v[k*rx*ax+(j+1)*ax-1]-v[k*rx*ax+(j+1)*ax-2];

            /*Calculate data-dependent parameter*/
            dp[0]=d[0]/2;
            for (i=1;i<ax-1;i++)
            {
                dp[i]=(d[i]-dp[i-1])/(4-cp[i-1]);
            }
            dp[ax-1]=(d[ax-1]-dp[ax-2])/(2-cp[ax-2]);

            /*Calculate intermediate parameters*/
            kp[ax-1]=dp[ax-1];
            for (i=ax-2;i>=0;i--)
            {
                kp[i]=dp[i]-cp[i]*kp[i+1];
            }

            /*Calculate a and b parameters*/
            for (i=0;i<ax-1;i++)
            {
                a[i]=3*kp[i]-v[k*rx*ax+j*ax+i+1]+v[k*rx*ax+j*ax+i];
                b[i]=-3*kp[i+1]+v[k*rx*ax+j*ax+i+1]-v[k*rx*ax+j*ax+i];
            }

            /*Rescale interp vector to indices, calculate output value*/
            for (i=0;i<ninterp;i++)
            {
                if((channel==1&&(mask[k*ninterp+i]==0||vout[j*ninterp+i]==nan))||(channel==0&&(mask[k*ninterp+i]==0||vout[i]==nan))){}
                else
                {
                    cur=(tt[k*ninterp+i]+tr[j*ninterp+i]-x1)/dx;
                    if(cur > ax-1 || cur < 0)
                    { 
                        if(channel==0)
                            vout[i]=nan;
                        else
                            vout[j*ninterp+i]=nan;
                    }
                    else
                    {
                        whole=(int)cur; /*Cast instead of floor, fine for positive numbers*/
                        t=cur-whole;
                        if(channel==0)
                            vout[i]=vout[i]+(1-t)*v[k*rx*ax+j*ax+whole]+t*v[k*rx*ax+j*ax+whole+1]+t*(1-t)*(a[whole]*(1-t)+b[whole]*t);
                        else
                            vout[j*ninterp+i]=vout[j*ninterp+i]+(1-t)*v[k*rx*ax+j*ax+whole]+t*v[k*rx*ax+j*ax+whole+1]+t*(1-t)*(a[whole]*(1-t)+b[whole]*t);
                    }
                }
            }
        }
    }
	
	/*Clean up dynamic array*/
	mxFree(d);
    mxFree(cp);
    mxFree(dp);
	mxFree(kp);
    mxFree(a);
    mxFree(b);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /*Input variables*/
   float x1, dx, nan;
   float *vout, *v, *tt, *tr;
   int spline;
   bool *mask;
  
   /*Local variables*/
   int ax, rx, tx, ninterp, ndims, channel;
   const mwSize *dim_array;
  
   /*Get input sizes*/
   dim_array=mxGetDimensions(I_V);
   ndims=mxGetNumberOfDimensions(I_V);
   ax=dim_array[0];
   rx=dim_array[1];
   if(ndims==2){tx=1;}
   else{tx=dim_array[2];}
   ninterp = mxGetM(I_XQ1);
   
   /*Get input data*/
   x1       = (float)mxGetScalar(I_X1);
   dx       = (float)mxGetScalar(I_DX);
   nan      = (float)mxGetScalar(I_NAN);
   v        = (float*)mxGetData(I_V);
   tt       = (float*)mxGetData(I_XQ1);
   tr       = (float*)mxGetData(I_XQ2);
   spline   = (int)mxGetScalar(I_SPLINE);
   mask     = (bool*)mxGetLogicals(I_MASK);
   channel  = (int)mxGetScalar(I_CHANNEL);
   
   /*Set up the output*/
   if(channel==0)
       O_VOUT = mxCreateNumericMatrix(ninterp,1,mxSINGLE_CLASS,mxREAL);
   else
       O_VOUT = mxCreateNumericMatrix(ninterp,rx,mxSINGLE_CLASS,mxREAL);
   vout  = (float*)mxGetData(O_VOUT);
   
   if(spline==0)
   {
       linterp(vout,x1,dx,v,tt,tr,mask,ax,tx,rx,ninterp,nan,channel);
   }
   else
   {
       sinterp(vout,x1,dx,v,tt,tr,mask,ax,tx,rx,ninterp,nan,channel);
   }
   return;
}
