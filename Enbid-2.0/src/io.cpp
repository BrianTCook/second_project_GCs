#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include<fstream>
#include "tree.h"
#include "functions.h"


/* This wrapper function select the desired output
 * routine for snapshot files.
 */
void savepositions(void )
{
    int i,j;

    for(i=1; i<=NumPart; i++) /*  start-up initialization */
    {
	for(j=0;j<ND;j++)
	{
	    P[i].Pos[j]=P[i].Pos[j]*All.hs[j];
	}
    }

    printf("\nwriting snapshot file ....\n");

    switch(All.ICFormat)
    {
	case 0:savepositions_ioformat0( );break;
	case 1:savepositions_ioformat1( );break;
	case 2:savepositions_ioformat2( );break;
	default:cout<<"Parameter ICFormat should be between 0 to 2"<<endl;endrun(10);
    }

    printf("done with snapshot.\n");

}


//This function writes ASCII format file
void savepositions_ioformat0( )
{
    int i;
    char buf[100];
    char file_suf[5]=".est";
    sprintf(buf,"%s%s%s",All.InitCondFile,All.SnapshotFileBase,file_suf);
    fprintf(stdout,"%s\n",buf); fflush(stdout);

    ofstream fd (buf);
    if (fd.is_open())
    {
	for(i=1;i<=NumPart;i++)
	{
	    fd<<P[i].Density<<endl;
	}
	fd.close();
    }

}


/* This function writes a density in Gadget type of  format
 * Each file contains a header first, then particle ID's if flag_id is enabled
 * in header, then the density of particles.
 */
void savepositions_ioformat1( )
{
    FILE *fd;
    char buf[100];
    float dummy[3];
    int i,j;
    int   blklen;
    int patterns[]={4,1,0},patternh[]={4,6,8,8,4,10,8,4,4,3,4,21,0};

#define BLKLEN my_fwrite1(&blklen,patterns,1,All.flag_swap,fd);

    char file_suf[5]=".est";
    sprintf(buf,"%s%s%s",All.InitCondFile,All.SnapshotFileBase,file_suf);

    fprintf(stdout,"%s\n",buf); fflush(stdout);
    if((fd=fopen(buf,"w")))
    {

	blklen=sizeof(header1);
	BLKLEN;
	my_fwrite1(&header1,patternh,1,All.flag_swap,fd);
	BLKLEN;

//    blklen=NumPart*3*sizeof(float);

	header1.flag_dim=ND;

	if(header1.flag_density==1)
	{
	    /* densities of particles */
	    blklen=NumPart*sizeof(float);  /* added density  */
	    BLKLEN;
	    for(i=1;i<=NumPart;i++)
	    {
		dummy[0]=P[i].Density;
		my_fwrite1(&dummy[0],patterns,1,All.flag_swap,fd);
	    }
	    BLKLEN;
	}

	if(header1.flag_density>1)
	{
	    for(j=0;j<header1.flag_density;j++)
	    {
		blklen=NumPart*sizeof(float);  /* added density  */
		BLKLEN;
		for(i=1;i<=NumPart;i++)
		{
//		    dummy[0]=Pa[i].QntSm[j];
		    my_fwrite1(&dummy[0],patterns,1,All.flag_swap,fd);
		}
		BLKLEN;
	    }
	}

	fclose(fd);
    }
    else
    {
	fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
	endrun(10);
    }
}



// This function can be used to write a new  output file
void savepositions_ioformat2( )
{
    int i;
    char buf[100];
    char file_suf[5]=".est";
    sprintf(buf,"%s%s%s",All.InitCondFile,All.SnapshotFileBase,file_suf);
    fprintf(stdout,"%s\n",buf); fflush(stdout);

    //write your own output routine here
    // below is output in ascii format

    ofstream fd (buf);
    if (fd.is_open())
    {
	for(i=1;i<=NumPart;i++)
	{
	    fd<<P[i].Density<<endl;
	}
	fd.close();
    }
}







size_t my_fwrite1(void *ptr,int* pattern,size_t nmemb,int flag_swap,FILE *stream)
{
    size_t nwritten,size,i;
//    int i;
    void* addr;

    for(i=0,size=0;pattern[i]>0;i+=2)
	size+=pattern[i]*pattern[i+1];

    if (flag_swap==1)
    {
	for(i=0;i<nmemb;i++)
	{
	    addr=(char *)ptr+i*size;
	    SwapEndian(addr,pattern);
	}
    }


    if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
	printf("I/O error (fwrite) on has occured.\n");
	fflush(stdout);
	endrun(777);
    }

    if (flag_swap==1)
    {
	for(i=0;i<nmemb;i++)
	{
	    addr=(char *)ptr+i*size;
	    SwapEndian(addr,pattern);
	}
    }

    return nwritten;
}



size_t my_fread1(void *ptr, int *pattern, size_t nmemb,int flag_swap,FILE *stream)
{
    size_t nread,i;
    void *addr;
    int size;

    for(i=0,size=0;pattern[i]>0;i+=2)
	size+=pattern[i]*pattern[i+1];


    if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
	printf("I/O error (fread) has occured.\n");
	fflush(stdout);
	endrun(778);
    }


    if (flag_swap==1)
    {
	for(i=0;i<nmemb;i++)
	{
	    addr=(char *)ptr+i*size;
	    SwapEndian(addr,pattern);
	}
    }

    return nread;
}

void SwapEndian(void* addr, int* pattern)
{
    int i,j=0,k;
    char c;
    while (pattern[j]>0)
    {
        for(k=0;k<pattern[j+1];k++)
        {
            for(i=0;i<pattern[j]/2;i++)
            {
                c=*((char*)addr+i);
                *((char*)addr+i)=*((char*)addr+(pattern[j]-i-1));
                *((char*)addr+(pattern[j]-i-1))=c;
            }
            addr=((char *)addr)+pattern[j];
        }
        j+=2;
    }
}




