#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "functions.h"
#include "nr.h"
#include "tree.h"


/*  this function wraps the distance x to the closest image
 *  for the given box size
 */

/* matrix multiplication of two matrices C=AB*/
void m_multiply(float **a, float **b,float **c,int n1, int n2,int transpose)
{
    int i,j,k;

    for (i=n1;i<=n2;i++)
	for (j=n1;j<=n2;j++)
	    for (k=n1,c[i][j]=0.0;k<=n2;k++)
	    {
		switch(transpose)
		{
		    case 0: c[i][j]+=a[i][k]*b[k][j]; break;
		    case 1: c[i][j]+=a[k][i]*b[k][j]; break;
		    case 2: c[i][j]+=a[i][k]*b[j][k]; break;
		    case 3: c[i][j]+=a[k][i]*b[j][k]; break;
		    default: cout<<"give m_multiply options correctly"<<endl;endrun(10); break;
		}
	    }

}

/* matrix operation on a vector C=AB*/
void m_multiply_a(float **a, float *b,float *c,int n1, int n2,int transpose)
{
    int i,k;

    for (i=n1;i<=n2;i++)
	for (k=n1,c[i]=0.0;k<=n2;k++)
	{

	    if(transpose)
		c[i]+=a[k][i]*b[k];
	    else
		c[i]+=a[i][k]*b[k];
//       switch(transpose)
//       {
//       case 0: c[i]+=a[i][k]*b[k]; break;
//       case 1: c[i]+=a[k][i]*b[k]; break;
//       default: cout<<"give m_multiply options correctly"<<endl;endrun(10); break;
//       }
	}

}

/* print matrix */
void m_print(float **a,int n1 ,int n2)
{
    int i,j;
    for (i=n1;i<=n2;i++)
	for (j=n1;j<=n2;j++)
	{
	    printf("%12.6f ",a[i][j]);
	    if(j==n2) cout<<endl;
	}
}

/* copy matrix */
void m_copy(float **a,float **b,int n1, int n2,int transpose)
{
    int i,j;
    for (i=n1;i<=n2;i++)
	for (j=n1;j<=n2;j++)
	{
	    if (transpose)      b[i][j]=a[j][i];
	    else                    b[i][j]=a[i][j];
	}
}


/* copy array */
void a_copy(float *a,float *b,int n1, int n2)
{
    int i;
    for (i=n1;i<=n2;i++)
	b[i]=a[i];

}



/* print array */
void a_print(float *a,int n1,int n2)
{
    int i;
    for (i=n1;i<=n2;i++)
    {
	printf("%12.6e ",a[i]);
	if(i==n2) cout<<endl;
    }
}

/* max element array */
float a_max(float *a,int n1 , int n2)
{
    int i;
    float amax=a[n1];
    for (i=n1;i<=n2;i++)
    {
	if(a[i] > amax ) amax=a[i];
    }
    return amax;
}

/* min element array */
float a_min(float *a,int n1, int n2)
{
    int i;
    float amin=a[n1];
    for (i=n1;i<=n2;i++)
    {
	if(a[i] < amin ) amin=a[i];
    }
    return amin;
}





double second(void)
{
    return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

    /* note: on AIX and presumably many other 32bit systems,
     * clock() has only a resolution of 10ms=0.01sec
     */
}


/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0,double t1)
{
    double dt;

    dt=t1-t0;

    if(dt<0)  /* overflow has occured */
    {
	dt=t1 + pow(2.0,32.0)/CLOCKS_PER_SEC - t0;
    }

    return dt;
}



/* returns the maximum of two double
 */
double dmax(double x,double y)
{
    if(x>y)
	return x;
    else
	return y;
}

/* returns the minimum of two double
 */
double dmin(double x,double y)
{
    if(x<y)
	return x;
    else
	return y;
}

/* returns the maximum of two integers
 */
int imax(int x,int y)
{
    if(x>y)
	return x;
    else
	return y;
}


/* returns the minimum of two integers
 */
int imin(int x,int y)
{
    if(x<y)
	return x;
    else
	return y;
}



/* end the run */
void endrun(int ierr)
{
    if(ierr)
    {
	fprintf(stdout,"endrun called with an error level of %d\n\n\n", ierr);
	exit(1);
    }
    exit(0);
}



void checkrun(string s1 ,string s2, string s3)
{
    int temp;
    cout<<s1<<endl;
    cout<<"Type 0 to "<<s2<<" or 1 to "<<s3<<" :";
    cin>>temp;
    if(temp!=1)
	exit(1);
}


void endrunm(string s )
{
    cout<<s<<endl;
    exit(1);
}



/* calculate the gamma function */
float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

