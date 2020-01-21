#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "functions.h"

/*
 *  Does the initial set-up.
 *  Reading the parameterfile, setting units,
 *  getting IC's  etc.
 */
void begrun(void)
{
    t11=second();
    read_parameter_file(ParameterFile);   /* ... read in parameters for this run */

    cout<<"Dimensions = "<<ND<<endl;

#ifdef MEDIAN
    if(All.MedianSplittingOn==0)
    {
	cout<<"MedianSplittingOn should be set to 1"<<endl;
	endrun(10);
    }
#endif


    set_sph_kernel(All.TypeOfKernel,ND);


    srand48(42);

    t22=second();  t00+=timediff(t11,t22);

    init();    /* ... read in initial model */
}



/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile. (Taken from GADGET)
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

    FILE *fd,*fdout;

    char buf[200],buf1[200],buf2[200],buf3[200];
    int  i,j,nt;
    int  id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int  errorFlag=0;

    nt=0;


    strcpy(tag[nt],"SpatialScale"); 
    addr[nt]=&All.SpatialScale;
    id[nt++]=DOUBLE;

    strcpy(tag[nt],"Anisotropy"); 
    addr[nt]=&All.Anisotropy;
    id[nt++]=DOUBLE;


//     strcpy(tag[nt],"InputDir"); // for batch processing input dir
//     addr[nt]=All.InputDir;
//     id[nt++]=STRING;

//     strcpy(tag[nt],"InitFileBase");  // for batch processing initial file base
//     addr[nt]=All.InitFileBase;
//     id[nt++]=STRING;



    strcpy(tag[nt],"TypeOfSmoothing"); 
    addr[nt]=&All.TypeOfSmoothing;
    id[nt++]=INT;

    strcpy(tag[nt],"AnisotropicKernel"); 
    addr[nt]=&All.AnisotropicKernel;
    id[nt++]=INT;

    strcpy(tag[nt],"CubicCells"); 
    addr[nt]=&All.CubicCells;
    id[nt++]=INT;


    strcpy(tag[nt],"VolCorr"); 
    addr[nt]=&All.VolCorr;
    id[nt++]=INT;

    strcpy(tag[nt],"KernelBiasCorrection"); 
    addr[nt]=&All.KernelBiasCorrection;
    id[nt++]=INT;

    strcpy(tag[nt],"TypeOfKernel");    
    addr[nt]=&All.TypeOfKernel;
    id[nt++]=INT;

    strcpy(tag[nt],"PartBoundary"); 
    addr[nt]=&All.PartBoundary;
    id[nt++]=INT;

    strcpy(tag[nt],"NodeSplittingCriterion"); 
    addr[nt]=&All.NodeSplittingCriterion;
    id[nt++]=INT;

    strcpy(tag[nt],"MedianSplittingOn"); 
    addr[nt]=&All.MedianSplittingOn;
    id[nt++]=INT;


    strcpy(tag[nt],"InitCondFile");
    addr[nt]=All.InitCondFile;
    id[nt++]=STRING;

    strcpy(tag[nt],"TypeListOn"); 
    addr[nt]=&All.TypeListOn;
    id[nt++]=INT;

    strcpy(tag[nt],"SnapshotFileBase");
    addr[nt]=All.SnapshotFileBase;
    id[nt++]=STRING;


    strcpy(tag[nt],"DesNumNgb");
    addr[nt]=&All.DesNumNgb;
    id[nt++]=INT;


    strcpy(tag[nt],"DesNumNgbA");
    addr[nt]=&All.DesNumNgbA;
    id[nt++]=INT;


    strcpy(tag[nt],"ICFormat");
    addr[nt]=&All.ICFormat;
    id[nt++]=INT;


    strcpy(tag[nt],"PeriodicBoundaryOn"); // read from file periodic_lengths.txt
    addr[nt]=&All.PeriodicBoundaryOn;
    id[nt++]=INT;





    if((fd=fopen(fname,"r")))
    {
	sprintf(buf,"%s%s",fname,"_enbid-usedvalues");
	if(!(fdout=fopen(buf,"w")))
	{
	    fprintf(stdout,"error opening file '%s' \n",buf);
	    errorFlag=1;
	}
	else
	{
	    while(!feof(fd))
	    {
		fgets(buf,200,fd);
		if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
		    continue;

		if(buf1[0]=='%')
		    continue;

		for(i=0,j=-1;i<nt;i++)
		    if(strcmp(buf1,tag[i])==0)
		    {
			j=i;
			tag[i][0]=0;
			break;
		    }

		if(j>=0)
		{
		    switch(id[j])
		    {
			case DOUBLE:
			    *((double*)addr[j])=atof(buf2);
			    fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
			    break;
			case STRING:
			    strcpy((char *)addr[j],buf2);
			    fprintf(fdout,"%-35s%s\n",buf1,buf2);
			    break;
			case INT:
			    *((int*)addr[j])=atoi(buf2);
			    fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
			    break;
		    }
		}
		else
		{
		    fprintf(stdout,"Error in file %s:   Tag '%s' not allowed or multiple defined.\n",fname,buf1);
		    errorFlag=1;
		}
	    }
	}
	fclose(fd);
	fclose(fdout);

//	sprintf(buf1, "%s%s", fname, "_enbid-usedvalues");
    }
    else
    {
	fprintf(stdout,"Parameter file %s not found.\n", fname);
	errorFlag=1;
	endrun(1);
    }


    for(i=0;i<nt;i++)
    {
	if(*tag[i])
	{
	    fprintf(stdout,"Error. I miss a value for tag '%s' in parameter file '%s'.\n",tag[i],fname);
	    errorFlag=1;
	}
    }






    if(errorFlag)
	endrun(1);


#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS
}


/* Here the lookup table for the kernel of the SPH calculation
 * is initialized.
 */
void set_sph_kernel(int TypeOfKernel, int dim)
{
    int i,nd1;
    double vd,f1=0.0,f2;
// more accurate 11/16/05
    float fsp[]={1.3333369,1.8189136,2.5464790,3.6606359,5.4037953,8.1913803,12.748839,20.366416,33.380983,56.102186,96.621159,170.39909,307.49826,567.37865,1069.6362,2058.8172,4043.1010,8095.3490,16515.921,34312.457};
    float fbw[]={0.93750176,0.95492964,1.0444543,1.2158542,1.4960706,1.9350925,2.6191784,3.6957561,5.4191207,8.2347774,12.937213,20.969717,35.003422,60.073937,105.84864,191.22182,353.77415,669.54824,1295.0185,2557.4981};

    float fcic[]={1.0000020,0.95492963,0.95492963,1.0132118,1.1398631,1.3545643,1.6932052,2.2174514,3.0316711,4.3134446,6.3690677,9.7358945,15.373955,25.030601,41.945765,72.238385,127.67576,231.29376,428.98102,813.72530};

    float
	fep[]={0.75000113,0.63661975,0.59683102,0.60792705,0.66492015,0.77403670,0.95242788,1.2319173,1.6674189,2.3527875,3.4499109,5.2424031,8.2360444,13.349647,22.283674,38.243824,67.384374,121.73344,225.21478,426.23651};

    float
	ftsc[]={1.6875038,1.9833154,2.4171656,3.0521959,3.9950015,5.4207369,7.6205882,11.087829,16.674494,25.880863,41.399134,68.151573,115.30531,200.24626,356.54234,650.15322,1212.9379,2312.9804,4504.4024,8951.2220};



    if((All.TypeOfSmoothing==4)||(All.TypeOfSmoothing==5)) nd1=1;  else nd1=dim;

    vd=2.0*pow(PI,ND/2.0)/(ND*exp(gammln(ND/2.0)));
    f2=1/4.0;
   
    switch (TypeOfKernel)
    {
	case 0:f1=fsp[nd1-1];break;          //Spline
	case 1:f1=1.0/vd;break;             //Top Hat
	case 2:f1=fbw[nd1-1];break;          //Bi-weight
	case 3:f1=fep[nd1-1];break;          //Epanechikov
	case 4:f1=fcic[nd1-1];break;          //CIC
	case 5:f1=ftsc[nd1-1];break;          //TSC
	case 6:f1=1.0/pow(f2*sqrt(2*PI),nd1);break;          //Gaussian
	default:cout<<"specify the Kernel Type"<<endl;endrun(10);
    }

    if((ND<1)||(ND>20))
    {
	cout<<"Specify the Normalization constant for dimensions > 20"<<endl;
	;endrun(10);
    }
    else
    {
	cout<<"Normalization constant of Kernel type "<< TypeOfKernel<<": "<<f1<<endl;
    }


    for(i=0;i<=KERNEL_TABLE+1;i++)
	KernelRad[i] = ((double)i)/KERNEL_TABLE;

    Kernel[KERNEL_TABLE+1] = KernelDer[KERNEL_TABLE+1]= 0;
    for(i=0;i<=KERNEL_TABLE;i++)
    {

	if(TypeOfKernel==0)
	{
	    if(KernelRad[i]<=0.5)
	    {
		Kernel[i] = f1 *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
		KernelDer[i] = f1 *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
		KernelDer2[i] = f1 *( -12 + 36*KernelRad[i]);
	    }
	    else
	    {
		Kernel[i] = f1 * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
		KernelDer[i] = f1 *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
		KernelDer2[i] = f1 *( 12*(1-KernelRad[i]));
	    }
	}

	if(TypeOfKernel==1)
	{
	    Kernel[i] = f1;
	}

	if(TypeOfKernel==2)
	{
	    Kernel[i] = f1 *(1-KernelRad[i]*KernelRad[i])*(1-KernelRad[i]*KernelRad[i]);
	}
	if(TypeOfKernel==3)
	{
	    Kernel[i] = f1 *(1-KernelRad[i]*KernelRad[i]);
	}
	if(TypeOfKernel==4)
	{
	    Kernel[i] = f1 *(1-KernelRad[i]);
	}
	if(TypeOfKernel==5)
	{

	    if(KernelRad[i]<1.0/3.0)
	    {
		Kernel[i] = f1 *(2.0/3.0-2*KernelRad[i]*KernelRad[i]);
	    }
	    else
	    {
		Kernel[i] = f1 *(1-KernelRad[i])*(1-KernelRad[i]);
	    }
	}
	if(TypeOfKernel==6)
	{
	    Kernel[i] = f1 *exp(-KernelRad[i]*KernelRad[i]/(2*f2*f2));
	}
    }

}



