#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include<fstream>
#include "functions.h"


/* This function reads initial conditions that are in the default file format
 * of Gadget or in ASCII format x[i] y[i] z[i].
 */
/* Read GADGET format file */
void read_ic1(char *fname)
{
#define SKIP my_fread1(&blklen,patterns,1,All.flag_swap,fd);
    FILE *fd;
    int   i,k,massflag,count;
    float dummy[3];
    int   pc,type ;
    int patternd[]={4,3,0},patterns[]={4,1,0},patternh[]={4,6,8,8,4,10,8,4,4,3,4,21,0};
    int4byte  intdummy, blklen;
    All.flag_swap=0;
  
  
     if((ND!=3)&&(ND!=6))
 	endrunm("\nFor Gadget format file number of dimensions should be 3 or 6\n");

    cout<<"filename is: " << fname << "\n";
  
    if((fd=fopen(fname,"r")))
    {
	fprintf(stdout,"Reading GADGET format file: %s \n",fname); fflush(stdout);
    
	SKIP;
	if(blklen!=256)
	{
	    All.flag_swap=1;
	    SwapEndian(&blklen,patterns);
	    if(blklen!=256)
	    {
		printf("incorrect header format (2)\n");
		endrun(889);
	    }
	}

	my_fread1(&header1,patternh,1,All.flag_swap,fd);


	SKIP;
	if(blklen!=256)
	{
	    printf("incorrect header format (2)\n");
	    endrun(889);
	}
    
//	All.BoxSize=header1.BoxSize;
	for(i=0, massflag=0;i<6;i++)
	{
	    if(header1.mass[i]==0 && header1.npart[i]>0)
		massflag=1;
	}
    
    
	for(i=0,NumPart=0;i<6;i++)
	    NumPart+=header1.npart[i];

	All.MaxPart = NumPart;    // sets the maximum number of particles 
				    
	allocate_memory();
    
	SKIP;
	for(i=1;i<=NumPart;i++)
	{
	    my_fread1(&dummy[0],patternd,1,All.flag_swap,fd);
		for(k=0;k<3;k++)
		    P[i].Pos[k]=dummy[k];	
	}
	SKIP;
    
    
	SKIP;
	for(i=1;i<=NumPart;i++)
	{
	    my_fread1(&dummy[0],patternd,1,All.flag_swap,fd);
#ifdef DIM3
	    for(k=0;k<3;k++)
		P[i].Vel[k]=dummy[k];
#else
	    if(ND>3)
	    for(k=0;k<(ND-3);k++)
		P[i].Pos[k+3]=dummy[k];
#endif
	}
	SKIP;
    
    
	SKIP;
	for(i=1;i<=NumPart;i++)
	{
	    my_fread1(&intdummy, patterns, 1,All.flag_swap, fd);
	    P[i].ID= intdummy;
	}
	SKIP;
    
    
	if(massflag)
	    SKIP;
	for(type=0, count=1; type<6; type++)
	{
	    if(header1.mass[type]==0 && header1.npart[type]>0)
	    {
		k=count;
		for(i=1;i<=header1.npart[type];i++)
		{
		    my_fread1(&dummy[0], patterns, 1,All.flag_swap, fd);
		    P[count++].Mass=dummy[0];
          
		}
        
	    }
	    else
	    {
		for(i=1;i<=header1.npart[type];i++)
		{
		    P[count++].Mass= header1.mass[type];
		}
	    }
	}
	if(massflag)
	    SKIP;
    
	fclose(fd);
	fprintf(stdout,"....done with reading.\n"); fflush(stdout);
    

	npart.clear();
	npartc.clear();
	for(i=0,k=0; i<6; i++)
	{
	    npart.push_back(header1.npart[i]);
	    npartc.push_back(k);
	    k+=npart[i];
	}
	/* set the particle types */
	for(type=0, pc=1; type<6; type++)
	    for(i=0; i<npart[type]; i++)
		P[pc++].Type = type;

    
    }
    else
    {
	fprintf(stdout,"File %s not found.\n", fname);
	endrun(7);
    }
  
    for(i=0; i<6; i++)
	if(header1.npart[i]>0) cout<<"Type = "<<i<<" Particles = "<<header1.npart[i]<<endl;
  
    fprintf(stdout,"Total particles = %d\n", NumPart);
  
#undef SKIP    
  
}











/* Read ASCII format file */
void read_ic0(char *fname)
{
#define SKIP my_fread(&blklen,sizeof(int4byte),1,fd);
    // FILE *fd;
    int   i,k,m=0;
    float temp;
    int   pc,type ;
  
    //double u_init;
    //ofstream fd (fname);
  
    cout<<"\n gets here in reading ascii file \n";    
    cout<<"\n filename is " << fname << "\n";

    ifstream fd;
    fd.open(fname);
    if (fd.is_open())
    {
	fprintf(stdout,"Reading ASCII format file: %s \n",fname); fflush(stdout);
	m=0;
	while(fd>>temp) m++;
	cout << "m: " << m << "\n";
	cout << "ND: " << ND << "\n";
	if((m%ND)==0)
	    cout<<"Read "<<m<<" records"<<endl;
	else
	{
	    cout<<"File format incorrect expecting "<<m%ND<<" more records"<<endl;
	    endrun(10);
	}
	fd.close();
    }
    else
    {
	fprintf(stdout,"File %s not found.\n", fname);
	endrun(7);
    }
  
    NumPart=m/ND;
    All.MaxPart = NumPart;    // sets the maximum number of particles 
				 
    for(i=0;i<6;i++)
    {
	header1.npart[i]=0;
	header1.npartTotal[i]=0;
	header1.mass[i]=0;
    }
    header1.npart[1]=NumPart;
    header1.mass[1]=1.0;
    header1.npartTotal[1]=NumPart;
  
  
    allocate_memory();
  
  
    /* reading data */
    ifstream fd1;
    fd1.open(fname);
    if (fd1.is_open())
    {
	for(i=1; i<=NumPart; i++) /*  start-up initialization */
	{
	    P[i].ID=i;
	    for(k=0;k<ND;k++)
		if(!(fd1>>P[i].Pos[k]))
		{
		    cout<<" file reading problem terminatiing"<<i<<" "<<k<<endl;
		    //                endrun(10);
		}
	}
	fd1.close();
    }
  
    for(i=1; i<=NumPart; i++) /*  start-up initialization */
	P[i].Mass=1.0;
    /* set the particle types */
	npart.clear();
	npartc.clear();
	for(i=0,k=0; i<6; i++)
	{
	    npart.push_back(header1.npart[i]);
	    npartc.push_back(k);
	    k+=npart[i];
	}



    for(type=0, pc=1; type<6; type++)
	for(i=0; i<header1.npart[type]; i++)
	    P[pc++].Type = type;
  
  
    for(i=0; i<6; i++)
	if(header1.npart[i]>0) cout<<"Type = "<<i<<" Particles = "<<header1.npart[i]<<endl;
  
    fprintf(stdout,"Total particles = %d\n", NumPart);
  
#undef SKIP  
}

/* write  a new reading routine */
void read_ic2(char *fname)
{
  
    int i,k,type,pc;
    ifstream fd;
    fd.open(fname);
    if (fd.is_open())
    {
	fprintf(stdout,"Reading --- format file: %s \n",fname); fflush(stdout);
	// read here

	fd.close();
    }
    else
	cout<<"Error opening file"<<endl;
  
    cout<<"Write your input routine here"<<endl;
    endrun(10);
  
    //-------------------------------------
    /* SPECIFY TOTAL NUMBER OF PARTICLES*/
    NumPart=1;
    //-------------------------------------
  
/* leave this unchanged */
    All.MaxPart = NumPart;
    for(i=0;i<6;i++)
    {
	header1.npart[i]=0;
	header1.npartTotal[i]=0;
	header1.mass[i]=0;
    }
    header1.npart[1]=NumPart;
    header1.npartTotal[1]=NumPart;
    allocate_memory();
  
    //-------------------------------------
    /* READ THE DATA HERE AND ASSIGN IT TO P[i].Pos */
    for(i=1; i<=NumPart; i++)
	for(k=0;k<ND;k++)
	    P[i].Pos[k]=1.0;
    for(i=1; i<=NumPart; i++)
    {
	    P[i].ID=i;
	P[i].Mass=1.0;
    }
    //-------------------------------------
  
  

	npart.clear();
	npartc.clear();
	for(i=0,k=0; i<6; i++)
	{
	    npart.push_back(header1.npart[i]);
	    npartc.push_back(k);
	    k+=npart[i];
	}



    /* set the particle types */
    for(type=0, pc=1; type<6; type++)
	for(i=0; i<header1.npart[type]; i++)
	    P[pc++].Type = type;
  
    for(i=0; i<6; i++)
	if(header1.npart[i]>0) cout<<"Type = "<<i<<" Particles = "<<header1.npart[i]<<endl;
    fprintf(stdout,"Total particles = %d\n", NumPart);
  
  
  
}







/* Read typelist file */
void read_typelist(char *fname,vector<int> &npart, struct particle_data P[])
{
    //  FILE *fd;
    int   i,k=0;
//    unsigned int j;
    int temp;
    int   pc,type ;
  
    ifstream fd;
    fd.open(fname);
    if (fd.is_open())
    {
	npart.clear();
	npartc.clear();
	fprintf(stdout,"Reading TypeList file: %s \n",fname); fflush(stdout);
	while(fd>>temp) 
	{
	    npart.push_back(temp);
	    npartc.push_back(k);
	    k+=temp;
	}
	fd.close();
	cout<<"Total No of particle types= "<<npart.size()<<" and Total Particles "<<k<<endl;
    }
    else
    {
	cout<<"Error opening typelist file"<<endl;
	fprintf(stdout,"%s \n",fname); fflush(stdout);
	endrun(1);
    }

//    if(k!=NumPart)  

	for(type=0, pc=1; type<int(npart.size()); type++)
	    for(i=0; i<npart[type]; i++)
		P[pc++].Type = type;
    
}



/* Read periodic lengths file */
void read_periodic_lengths(char *fname,double boxh[ND])
{
    //  FILE *fd;
    int   k=0;
//    unsigned int j;
    float temp;
  
    ifstream fd;
    fd.open(fname);
    if (fd.is_open())
    {
	fprintf(stdout,"Reading periodic lengths file: %s \n",fname); fflush(stdout);

	for(k=0;k<ND;k++)
	{
	    if(fd>>temp) 
	    {
		boxh[k]=temp/2.0;
	    }
	    else
	    {
		cout<<"Error reading periodic lengths "<<endl;
		cout<<"Read "<<k<<" out of "<<ND<<" records"<<endl;
		endrun(1);
	    }
	    
	}
	fd.close();
    }
    else
    {
	cout<<"periodic_lengths.txt file not found"<<endl;
	cout<<"should be in the curent directory"<<endl;
	endrun(1);
    }
}

