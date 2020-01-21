/***************************************************************************
                          ngb_search.cpp  -  description
                             -------------------
    begin                : Mon Jan 16 2006
    copyright            : (C) 2006 by Sanjib Sharma
    email                : ssharma@aip.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <algorithm>
#include "allvars.h"
#include "nr.h"
#include "tree.h"
#include "functions.h"


#define adjust_heap(x,size1)\
{\
     hi=0;\
    while(1)\
    {\
	hj=((hi+1)<<1);\
	if(hj>size1) break;\
	--hj;\
	if(hj<(size1-1))\
	    if(x[hj]<x[hj+1])\
		++hj;\
	if(x[hi]<x[hj])\
	{\
	    ptemp=x[hj];x[hj]=x[hi];x[hi]=ptemp;\
	    hi=hj;\
	}\
	     else\
	    break;\
    }\
    pqHead=x->pq;\
}

inline double periodic(double x,double l2)
{
    if(x > l2)
	x-=l2*2;
    if(x < -l2)
	x+=l2*2;
    return x;
}


class rmlist 
{
public:
    float r,Mass;
//    rmlist(float r1,float Mass1){r=r1; Mass=Mass1;}
};



void ngb_treesearch_sphere_metric(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool* imarkA)
{
    int k,l,i,cell,cp,ct;
    double temp1,temp2,temp3;
    bool b1=0;
    struct pqueue ptemp;
    int hj,hi;

//    cout<<nodes1<<endl;

    cell=no;
    while(cell!=ROOT)
    {
	cp=SIBLING(cell,nodes1);
	ct=cp;
	SETNEXT(ct,nodes1);

	while(1)
	{

#ifdef PERIODIC
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
	    {
		if(All.boxh[k]>0)
		{
		    temp1=(nodes1+cp)->bnd->x[k][1] - searchx[k][0];
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
	    {

		if(All.boxh[k]>0)
		{
		    temp1=searchx[k][1]-(nodes1+cp)->bnd->x[k][0]; 
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

#else
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
		goto GETNEXT;
	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
		goto GETNEXT;

	    for(k=0,temp3=0; k<ND; ++k)
	    {
		temp1=((nodes1+cp)->bnd->x[k][1]-searchcenter[k])*metric[k];
		temp2=((nodes1+cp)->bnd->x[k][0]-searchcenter[k])*metric[k];
		if(temp1>0) 
		    if(temp2<0) 
			continue;
//			break;
		if(temp1<0) {temp1=-temp1; temp2=-temp2;}
		if(temp1<temp2) temp2=temp1;
		temp3+=temp2*temp2;
		if(temp3>pqHead->r) goto GETNEXT;
	    }

#endif	    

	    if((nodes1+cp)->count>NumBucket)
	    {
		cp=LOWER(cp,nodes1);
		continue;
	    }
	    else
	    {
		b1=0;
		for(l=0; l<(nodes1+cp)->count; ++l)
		{
		    i=(nodes1+cp)->lid+l;
		    if(imarkA[i]) continue;
		    for(k=0,temp2=0; k<ND; ++k)
		    {
#ifdef PERIODIC
			temp1=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
			temp1=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
			temp2+=temp1*temp1;
		    }
		    if(temp2 < pqHead->r)
		    {
			imarkA[pqHead->p]=0; imarkA[i]=1;
			pqHead->r=temp2;	pqHead->p=i; 
			for(k=0;k<ND;++k)
			    pqHead->x[k]=Part[i].Pos[k];
 			adjust_heap(pqxA,DesNumNgb); b1=1;
		    }
		}
		if(b1)
		    for(k=0; k<ND; ++k)
		    {
			temp1=sqrt(pqHead->r)/metric[k];    
			searchx[k][0]=searchcenter[k]-temp1;
			searchx[k][1]=searchcenter[k]+temp1;
		    }
	    }

	GETNEXT:	  
	    SETNEXT(cp,nodes1);
	    if(cp==ct) break;
	}
	cell=PARENT(cell,nodes1);

    }

    return;
}


void ngb_treesearch_sphere_metric_exact_h(int no,struct NODE *nodes1, int NumBucket, vector<class rmlist> & rmlistv,bool* imarkA)
{
    int k,l,i,cell,cp,ct;
    double temp1,temp2,temp3;
    class rmlist rmliste;

//    cout<<nodes1<<endl;

    cell=no;
    while(cell!=ROOT)
    {
	cp=SIBLING(cell,nodes1);
	ct=cp;
	SETNEXT(ct,nodes1);

	while(1)
	{

#ifdef PERIODIC
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
	    {
		if(All.boxh[k]>0)
		{
		    temp1=(nodes1+cp)->bnd->x[k][1] - searchx[k][0];
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
	    {

		if(All.boxh[k]>0)
		{
		    temp1=searchx[k][1]-(nodes1+cp)->bnd->x[k][0]; 
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

#else
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
		goto GETNEXT;
	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
		goto GETNEXT;

	    for(k=0,temp3=0; k<ND; ++k)
	    {
		temp1=((nodes1+cp)->bnd->x[k][1]-searchcenter[k])*metric[k];
		temp2=((nodes1+cp)->bnd->x[k][0]-searchcenter[k])*metric[k];
		if(temp1>0) 
		    if(temp2<0) 
			continue;
//			break;
		if(temp1<0) {temp1=-temp1; temp2=-temp2;}
		if(temp1<temp2) temp2=temp1;
		temp3+=temp2*temp2;
		if(temp3>pqHead->r) goto GETNEXT;
	    }

#endif	    


	    if((nodes1+cp)->count>NumBucket)
	    {
		cp=LOWER(cp,nodes1);
		continue;
	    }
	    else
	    {
		for(l=0; l<(nodes1+cp)->count; ++l)
		{
		    i=(nodes1+cp)->lid+l;
		    if(imarkA[i]) continue;
		    for(k=0,temp2=0; k<ND; ++k)
		    {
#ifdef PERIODIC
			temp1=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
			temp1=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
			temp2+=temp1*temp1;
		    }
		    if(temp2 < 1.0)
		    {
			rmliste.r=temp2;
			rmliste.Mass=Part[i].Mass;
 			rmlistv.push_back(rmliste);
		    }
		}
	    }

	GETNEXT:	  
	    SETNEXT(cp,nodes1);
	    if(cp==ct) break;
	}
	cell=PARENT(cell,nodes1);

    }

    return;
}









void ngb_treesearch_sphere_gmatrix(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool * imarkA)
{


    int k,l,i,j,cell,cp,ct;
    double temp1,temp2,temp3;
    bool b1=0;
    struct pqueue ptemp;
    int hj,hi;

//    cout<<nodes1<<endl;

    cell=no;
    while(cell!=ROOT)
    {
	cp=SIBLING(cell,nodes1);
	ct=cp;
	SETNEXT(ct,nodes1);
	while(1)
	{

#ifdef PERIODIC
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
	    {
		if(All.boxh[k]>0)
		{
		    temp1=(nodes1+cp)->bnd->x[k][1] - searchx[k][0];
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
	    {

		if(All.boxh[k]>0)
		{
		    temp1=searchx[k][1]-(nodes1+cp)->bnd->x[k][0]; 
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

#else
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
		goto GETNEXT;
	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
		goto GETNEXT;

	    for(k=0,temp3=0; k<ND; ++k)
	    {
 		temp1=((nodes1+cp)->bnd->x[k][1]-searchcenter[k])*metric1[k];
 		temp2=((nodes1+cp)->bnd->x[k][0]-searchcenter[k])*metric1[k];
 		if(temp1>0) 
 		    if(temp2<0) 
 			continue;
// 			break;
		if(temp1<0) {temp1=-temp1; temp2=-temp2;}
		if(temp1<temp2) temp2=temp1;
		temp3+=temp2*temp2;
		if(temp3>pqHead->r) goto GETNEXT;
	    }

#endif

	    if((nodes1+cp)->count>NumBucket)
	    {
		cp=LOWER(cp,nodes1);
		continue;
	    }
	    else
	    {
		b1=0;
		for(l=0; l<(nodes1+cp)->count; ++l)
		{
		    i=(nodes1+cp)->lid+l;
		    if(imarkA[i]) continue;
		    for(k=0,temp2=0; k<ND; ++k)
		    {

#ifdef PERIODIC
 			for(j=0,temp1=0; j<ND; ++j)
			    temp1+=periodic((Part[i].Pos[j]-searchcenter[j]),All.boxh[j])*gmatrix[j][k]*metric1[j];
#else
 			for(j=0,temp1=0; j<ND; ++j)
 			    temp1+=(Part[i].Pos[j]-searchcenter[j])*gmatrix[j][k]*metric1[j];
#endif
			temp1=temp1*metric[k];
			temp2+=temp1*temp1;
		    }
		    if(temp2 < pqHead->r)
		    {
			imarkA[pqHead->p]=0; imarkA[i]=1;
			pqHead->r=temp2;	pqHead->p=i; 
			for(k=0;k<ND;++k)
			    pqHead->x[k]=Part[i].Pos[k];
 			adjust_heap(pqxA,DesNumNgb); b1=1;
		    }
		}
		if(b1)
		{
		    temp2=sqrt(pqHead->r);
		    for(k=0; k<ND; ++k)
		    {
			temp1=temp2/metric1[k];
			searchx[k][0]=searchcenter[k]-temp1;
			searchx[k][1]=searchcenter[k]+temp1;
		    }
		}
	    }

	GETNEXT:	  
	    SETNEXT(cp,nodes1);
	    if(cp==ct) break;
	}
	cell=PARENT(cell,nodes1);
    }

    return;
}


void ngb_treesearch_sphere_nometric(int no, struct NODE* nodes1,int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool * imarkA)
{
    int k,l,i,cell,cp,ct;
    double temp1,temp2,temp3;
    bool b1=0;
    struct pqueue ptemp;
    int hj,hi;

    cell=no;
    while(cell!=ROOT)
    {
	cp=SIBLING(cell,nodes1);
	ct=cp;
	SETNEXT(ct,nodes1);
	while(1)
	{
#ifdef PERIODIC
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
	    {
		if(All.boxh[k]>0)
		{
		    temp1=(nodes1+cp)->bnd->x[k][1] - searchx[k][0];
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
	    {

		if(All.boxh[k]>0)
		{
		    temp1=searchx[k][1]-(nodes1+cp)->bnd->x[k][0]; 
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

#else
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
		goto GETNEXT;
	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
		goto GETNEXT;

	    for(k=0,temp3=0; k<ND; ++k)
	    {
 		temp1=(nodes1[cp].bnd->x[k][1]-searchcenter[k]);
 		temp2=(nodes1[cp].bnd->x[k][0]-searchcenter[k]);
		if(temp1>0) 
		    if(temp2<0) 
			continue;
//			break;
		if(temp1<0) {temp1=-temp1; temp2=-temp2;}
		if(temp1<temp2) temp2=temp1;
		temp3+=temp2*temp2;
		if(temp3>pqHead->r) goto GETNEXT;
	    }
#endif
	    if((nodes1+cp)->count>NumBucket)
	    {
		cp=LOWER(cp,nodes1);
		continue;
	    }
	    else
	    {
		b1=0;
		for(l=0; l<(nodes1+cp)->count; ++l)
		{
		    i=(nodes1+cp)->lid+l;
		    if(imarkA[i]) continue;
		    for(k=0,temp2=0; k<ND; ++k)
		    {

#ifdef PERIODIC
			temp1=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k]);
#else
			temp1=(Part[i].Pos[k]-searchcenter[k]);
#endif
			temp2+=temp1*temp1;
		    }
		    if(temp2 < pqHead->r)
		    {
			imarkA[pqHead->p]=0; imarkA[i]=1;
			pqHead->r=temp2;	pqHead->p=i; 
			for(k=0;k<ND;++k)
			    pqHead->x[k]=Part[i].Pos[k];
 			adjust_heap(pqxA,DesNumNgb); b1=1;
		    }
		}
		if(b1)
		    for(k=0; k<ND; ++k)
		    {
			temp1=sqrt(pqHead->r);    
			searchx[k][0]=searchcenter[k]-temp1;
			searchx[k][1]=searchcenter[k]+temp1;
		    }
	    }

	GETNEXT:	  
	    SETNEXT(cp,nodes1);
	    if(cp==ct) break;
	}
	cell=PARENT(cell,nodes1);
    }

    return;
}



void ngb_treesearch_box_metric(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA, bool * imarkA)
{
    int k,l,i,cell,cp,ct;
    double temp1,temp2;
    bool b1=0;
    struct pqueue ptemp;
    int hj,hi;

//    cout<<nodes1<<endl;

    cell=no;
    while(cell!=ROOT)
    {
	cp=SIBLING(cell,nodes1);
	ct=cp;
	SETNEXT(ct,nodes1);
	while(1)
	{

#ifdef PERIODIC
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
	    {
		if(All.boxh[k]>0)
		{
		    temp1=(nodes1+cp)->bnd->x[k][1] - searchx[k][0];
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }

	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
	    {

		if(All.boxh[k]>0)
		{
		    temp1=searchx[k][1]-(nodes1+cp)->bnd->x[k][0]; 
		    if(temp1<All.boxh[k])
			goto GETNEXT;
		    else
			if((temp1-(All.boxh[k]*2))<0)
			    goto GETNEXT;
		}
		else
		    goto GETNEXT;
	    }
#else
	    k=(nodes1+(PARENT(cp,nodes1)))->k1;
	    if((nodes1+cp)->bnd->x[k][0] > searchx[k][1])
		goto GETNEXT;
	    if((nodes1+cp)->bnd->x[k][1] < searchx[k][0])
		goto GETNEXT;
#endif

	    if((nodes1+cp)->count>NumBucket)
	    {
		cp=LOWER(cp,nodes1);
		continue;
	    }
	    else
	    {
		b1=0;
		for(l=0; l<(nodes1+cp)->count; ++l)
		{
		    i=(nodes1+cp)->lid+l;
		    if(imarkA[i]) continue;
		    for(k=0,temp2=0; k<ND; ++k)
		    {
#ifdef PERIODIC
			temp1=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
			temp1=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
			if(temp1<0) temp1=-temp1;
			if(temp2<temp1) temp2=temp1;
		    }
		    if(temp2 < pqHead->r)
		    {
			imarkA[pqHead->p]=0; imarkA[i]=1;
			pqHead->r=temp2;	pqHead->p=i; 
			for(k=0;k<ND;++k)
			    pqHead->x[k]=Part[i].Pos[k];
 			adjust_heap(pqxA,DesNumNgb); b1=1;
		    }
		}
		if(b1)
		    for(k=0; k<ND; ++k)
		    {
			temp1=pqHead->r/metric[k];    
			searchx[k][0]=searchcenter[k]-temp1;
			searchx[k][1]=searchcenter[k]+temp1;
		    }
	    }

	GETNEXT:	  
	    SETNEXT(cp,nodes1);
	    if(cp==ct) break;
	}
	cell=PARENT(cell,nodes1);
    }

    return;
}




float  ngb_treedensity_sphere_metric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA)
{
    real u,wk,r;  /* search radius */
    int   i,j,k,th;
    float rhoxyz;
    double xmin[ND],xmax[ND],temp,temp1;
    struct linklist *pq;
    bool b1,b2=1;
    int b3=1;
    struct pqueue ptemp;
    int hj,hi;

    /* calculate metric */
    for(j=0; j<ND; j++)
	metric[j]=1.0;


    if(All.CubicCells==1)
    {

	if(ND==6)
	{

	    for(j=0; j<ND; j++)
		searchcenter[j]=0.0;
	    ngb_calculate_metric(xyz,nodes1);
	    for(j=0; j<ND; j++)
	    {
		metric[j]=searchcenter[j];
	    }


	    for(j=0,temp=1.0; j<3; j++)
		temp*=(metric[j]);
	    temp=pow(temp,0.3333);
	    for(j=0; j<3; j++)
		metric[j]=1/temp;


	    for(j=3,temp=1.0; j<6; j++)
		temp*=(metric[j]);
	    temp=pow(temp,0.3333);
	    for(j=3; j<6; j++)
		metric[j]=1/temp;
	}
    }
    else
    {
	for(j=0; j<ND; j++)
	    searchcenter[j]=0.0;
	ngb_calculate_metric(xyz,nodes1);
	for(j=0; j<ND; j++)
	{
	    metric[j]=1.0/searchcenter[j];
	}
    }

    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];




    while(b2)
    {
	/* identify the node containing the cell*/
	th=1;
	while(nodes1[th].count>desngb)
	{
 	    j=nodes1[th].k1;
	    if(xyz[j] < nodes1[(LOWER(th,nodes1))].bnd->x[j][1])
		th=(LOWER(th,nodes1));
	    else
		th=(UPPER(th,nodes1));
	}
 	while(nodes1[th].count<desngb)
	{
 	    th=PARENT(th,nodes1);

	}



	if(pnew)
	{

	    for(pq=pqStartA;pq<(pqStartA+desngb);pq++)
		imarkA[pq->p]=0;

	    for(pq=pqStartA,j=0;pq<(pqStartA+desngb);j++,pq++)
	    {
		i=nodes1[th].lid+j;
		pq->p=i;
		for(k=0; k<ND; k++)
		    pq->x[k]=Part[i].Pos[k];
		imarkA[i]=1;
	    }
	}



	for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
	{
	    for(j=0,pq->r=0.0; j<ND; j++)
	    {
#ifdef PERIODIC
		temp=periodic((pq->x[j]-searchcenter[j]),All.boxh[j])*metric[j];
#else
		temp=(pq->x[j]-searchcenter[j])*metric[j];
#endif

		pq->r+=temp*temp;
	    }
	    pqxA[desngb-i-1].pq=pq;
	}


	make_heap(pqxA,(pqxA+desngb));      pqHead=pqxA->pq;

	for(k=0; k<ND; k++)
	{
	    temp=sqrt(pqHead->r)/metric[k];    
	    searchx[k][0]=searchcenter[k]-temp;
	    searchx[k][1]=searchcenter[k]+temp;
	}


	b1=0;
	if (pnew) k=desngb; else k=0;
	for(j=0; j<nodes1[th].count; j++)
	{

	    i=nodes1[th].lid+j;
	    if(imarkA[i]) continue;
	    for(k=0,r=0; k<ND; k++)
	    {
#ifdef PERIODIC
		temp=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
		temp=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
		r+=temp*temp;
	    }


	    if(r < pqHead->r)
	    {
		imarkA[pqHead->p]=0; imarkA[i]=1;
		pqHead->r=r;	pqHead->p=i; 
		for(k=0;k<ND;k++)
		    pqHead->x[k]=Part[i].Pos[k];
		adjust_heap(pqxA,desngb); b1=1;
	    }

	}


	if(b1)
	    for(k=0; k<ND; k++)
	    {
		temp=sqrt(pqHead->r)/metric[k];    
		searchx[k][0]=searchcenter[k]-temp;
		searchx[k][1]=searchcenter[k]+temp;
	    }




	ngb_treesearch_sphere_metric(th,nodes1,desngb/2,desngb,pqxA,imarkA);


//----------------------------------------
	b2=0;
	if(b3<2)
	    if(All.VolCorr==1)
	    {

#ifdef PERIODIC
		for(i=0; i<desngb; i++)
		{

		    for(k=0; k<ND; k++)
		    {
			temp=periodic((pqStartA[i].x[k]-searchcenter[k]),All.boxh[k]);
			if((i==0)||(temp<xmin[k])) xmin[k]=temp;
			if((i==0)||(temp>xmax[k])) xmax[k]=temp;
		    }
		}

		    for(k=0; k<ND; k++)
		    {
			xmin[k]+=searchcenter[k];
			xmax[k]+=searchcenter[k];
		    
		    }
#else
		for(i=0; i<desngb; i++)
		{

		    for(k=0; k<ND; k++)
		    {
			if((i==0)||(pqStartA[i].x[k]<xmin[k])) xmin[k]=pqStartA[i].x[k];
			if((i==0)||(pqStartA[i].x[k]>xmax[k])) xmax[k]=pqStartA[i].x[k];
		    }
		}
#endif





		for(k=0; k<ND; k++)
		{
		    temp=sqrt(pqHead->r)/metric[k];    
		    searchx[k][0]=searchcenter[k]-temp;
		    searchx[k][1]=searchcenter[k]+temp;
		}

		for(k=0;k<ND;k++)
		{

		    if(All.CubicCells==1) 		    
			temp1=(xmax[k]-xmin[k])*25.0/(All.DesNumNgb);
		    else
			temp1=(xmax[k]-xmin[k])*25.0/(All.DesNumNgb);



		    if(((searchx[k][1]-xmax[k])>temp1)&&((xmin[k]-searchx[k][0])>temp1))
		    {			
			b2=1;
		    }
		}

		if(b2)
		    for(k=0;k<ND;k++)
		    {
			temp=(xmax[k]-xmin[k])*0.5;
			if(xmax[k]!=xmin[k])
			{			
			    metric[k]=1/temp;
			}
			else
			    metric[k]=1/(searchx[k][1]-searchx[k][0]);


		    }
		b3++;
	    }


    }

//----------------------------------------


    pnext=pqHead->p;
    r=pqHead->r;

    for(i=0,rhoxyz=0; i<desngb; i++)
    {
	u = sqrt(pqStartA[i].r/pqHead->r);
	if(u<1)
	{
	    if ((u==0) && (All.KernelBiasCorrection==1))
		u=ND/(1.0+ND);
	    k = (int)(u*KERNEL_TABLE);
	    wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);      rhoxyz+=wk*Part[pqStartA[i].p].Mass;
	}

	if(idoneA[pqStartA[i].p]==0)
	    if(pqStartA[i].r<r)
	    {pnext=pqStartA[i].p;r=pqStartA[i].r;}

//	    cout<<pqStartA[i].p<<" b "<<u<<endl;


    }
    if(idoneA[pnext]==1) pnew=1; else pnew=0;


    temp=sqrt(pqHead->r);
    for(k=0,r=1.0,wk=1.0;k<ND;k++)
    {
	r*=metric[k];
	wk*=temp;
    }
    temp=r/wk;
    //   cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;
//    cout<<rhoxyz<<" "<<All.hsv<<" "<<r/wk<<endl;

    return rhoxyz;
}






float  ngb_treedensity_bruteforce(float xyz[ND], int parttype,int desngb, struct pqueue* pqxA, struct linklist* pqStartA)
{
    real u,wk=0.0,r=0.0;  /* search radius */
    int   i,j,k,th,l;
    float rhoxyz=0.0;
    double xmin[ND],xmax[ND],temp,temp1,temp2;
    struct linklist *pq;
//     bool b1,b2=1;
//     int b3=1;
    struct pqueue ptemp;
    int hj,hi;


    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];

    th=1;
    for(pq=pqStartA,j=0;pq<(pqStartA+desngb);j++,pq++)
    {
	i=npartc[parttype]+j;
	pq->p=i;
	for(k=0; k<ND; k++)
	    pq->x[k]=Part[i].Pos[k];
    }

    for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
    {

	if(All.AnisotropicKernel==1)
	{
	    for(k=0,pq->r=0.0; k<ND; k++)
	    {
		for(l=0,temp=0; l<ND; ++l)
		{
#ifdef PERIODIC
		    temp+=periodic((pq->x[l]-searchcenter[l]),All.boxh[l])*gmatrix[l][k]*metric1[l];
#else
		    temp+=(pq->x[l]-searchcenter[l])*gmatrix[l][k]*metric1[l];
#endif
		}
		temp=temp*metric[k];
		pq->r+=temp*temp;

	    }
	}
	else
	{

	    for(j=0,pq->r=0.0; j<ND; j++)
	    {
#ifdef PERIODIC
		temp=periodic((pq->x[j]-searchcenter[j]),All.boxh[j])*metric[j];
#else
		temp=(pq->x[j]-searchcenter[j])*metric[j];
#endif
		if((All.TypeOfSmoothing==2)||(All.TypeOfSmoothing==3))
		    pq->r+=temp*temp;
		if((All.TypeOfSmoothing==4)||(All.TypeOfSmoothing==5))
		{
		    if(temp<0) temp=-temp;
		    if( pq->r<temp) pq->r=temp;
		}
	    }
	}

	pqxA[desngb-i-1].pq=pq;
    }

    make_heap(pqxA,(pqxA+desngb));      pqHead=pqxA->pq;


    for(j=desngb; j<npart[parttype]; ++j)
    {
	i=npartc[parttype]+j;

	if((All.TypeOfSmoothing==2)||(All.TypeOfSmoothing==3))
	{
	    if(All.AnisotropicKernel==1)
	    {

		for(k=0,temp2=0; k<ND; k++)
		{
		    for(l=0,temp=0; l<ND; ++l)
		    {
#ifdef PERIODIC
			temp+=periodic((Part[i].Pos[l]-searchcenter[l]),All.boxh[l])*gmatrix[l][k]*metric1[l];
#else
			temp+=(Part[i].Pos[l]-searchcenter[l])*gmatrix[l][k]*metric1[l];
#endif
		    }
		    temp=temp*metric[k];
		    temp2+=temp*temp;

		}

	    }
	    else
	    {


		for(k=0,temp2=0; k<ND; ++k)
		{
#ifdef PERIODIC
		    temp1=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
		    temp1=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
		    temp2+=temp1*temp1;
		}

	    }
	    if(temp2 < pqHead->r)
	    {
		pqHead->r=temp2;	pqHead->p=i; 
		adjust_heap(pqxA,desngb);
	    }
	}

	if((All.TypeOfSmoothing==4)||(All.TypeOfSmoothing==5))
	{
	    for(k=0,r=0; k<ND; k++)
	    {
#ifdef PERIODIC
		temp=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
		temp=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif

		if(temp<0) temp=-temp;
		if(r<temp) r=temp;
	    }

	    if(r < pqHead->r)
	    {
		pqHead->r=r;	pqHead->p=i; 
		for(k=0;k<ND;k++)
		    pqHead->x[k]=Part[i].Pos[k];
		adjust_heap(pqxA,desngb);
	    }

	}

    }
    

    

    if((All.TypeOfSmoothing==2)||(All.TypeOfSmoothing==3))
    {

//     vector<float> rt;
//     for(i=0,rhoxyz=0; i<desngb; i++)
// 	rt.push_back(sqrt(pqStartA[i].r/pqHead->r));
//     sort(rt.begin(),rt.end());


	for(i=0,rhoxyz=0; i<desngb; i++)
	{
	    u = sqrt(pqStartA[i].r/pqHead->r);
//	    u=rt[i];
	    if(u<1)
	    {
		if ((u==0) && (All.KernelBiasCorrection==1))
		    u=ND/(1.0+ND);
		k = (int)(u*KERNEL_TABLE);
		wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);      rhoxyz+=wk*(Part[pqStartA[i].p].Mass);
	    }
//	    cout<<i<<" a "<<pqStartA[i].p<<" a "<<u<<endl;


	}
	temp=sqrt(pqHead->r);

	for(k=0,r=1.0,wk=1.0;k<ND;k++)
	{
	    r*=metric[k];
	    if(All.AnisotropicKernel==1)
		r*=metric1[k];
	    wk*=temp;
	}

    }

    if((All.TypeOfSmoothing==4)||(All.TypeOfSmoothing==5))
    {
	for(i=0,rhoxyz=0; i<desngb; i++)
	{

	    for(k=0,r=0; k<ND; k++)
	    {
		xmax[k]=pqHead->r/metric[k];
#ifdef PERIODIC
		xmin[k]=periodic((pqStartA[i].x[k] - searchcenter[k]),All.boxh[k]);
#else
		xmin[k]=pqStartA[i].x[k] - searchcenter[k];
#endif
		temp=(xmin[k])/xmax[k];
		if(temp<0) temp=-temp;
		xmin[k]=temp;	       
		r+=temp;
	    }

	    for(l=0,wk=1.0;l<ND;l++)
	    {
		u=xmin[l];

		if ((r==0.0) && (All.KernelBiasCorrection==1))
		    u=0.5;
		if(u<1)
		{
		    k = (int)(u*KERNEL_TABLE);
		    wk*=( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);

		}
		else
		    wk=0.0;
	    }

	    rhoxyz+=Part[pqStartA[i].p].Mass*wk;
	}
	temp=(pqHead->r);

	for(k=0,r=1.0,wk=1.0;k<ND;k++)
	{
	    r*=metric[k];
	    wk*=temp;
	}

    }




    temp=r/wk;
//      cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;

//    cout<<rhoxyz/All.hsv<<" "<<r<<" "<<wk<<endl;

    return (rhoxyz/All.hsv);
}








float  ngb_treedensity_sphere_nometric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA)
{
    real u,wk,r;  /* search radius */
    int   i,j,k,th;
    float rhoxyz;
    double temp;
    struct linklist *pq;
    bool b1;
    struct pqueue ptemp;
    int hj,hi;


    /* calculate metric */

    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];


    /* identify the node containing the cell*/
    th=1;
    while(nodes1[th].count>desngb)
    {
	j=nodes1[th].k1;
	if(xyz[j] < nodes1[(LOWER(th,nodes1))].bnd->x[j][1])
	    th=(LOWER(th,nodes1));
	else
	    th=(UPPER(th,nodes1));
    }
    while(nodes1[th].count<desngb)
    {
	th=PARENT(th,nodes1);
    }


    if(pnew)
    {

	for(pq=pqStartA;pq<(pqStartA+desngb);pq++)
	    imarkA[pq->p]=0;

	for(pq=pqStartA,j=0;pq<(pqStartA+desngb);j++,pq++)
	{
	    i=nodes1[th].lid+j;
	    pq->p=i;
	    for(k=0; k<ND; k++)
		pq->x[k]=Part[i].Pos[k];
	    imarkA[i]=1;
	}
    }


    for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
    {
	for(j=0,pq->r=0.0; j<ND; j++)
	{

#ifdef PERIODIC
	    temp=periodic((pq->x[j]-searchcenter[j]),All.boxh[j]);
#else
	    temp=pq->x[j]-searchcenter[j];
#endif
	    pq->r+=temp*temp;
	}
	pqxA[desngb-i-1].pq=pq;
    }


    make_heap(pqxA,(pqxA+desngb));      pqHead=pqxA->pq;


    for(k=0; k<ND; k++)
    {
	temp=sqrt(pqHead->r);    
	searchx[k][0]=searchcenter[k]-temp;
	searchx[k][1]=searchcenter[k]+temp;
    }

    b1=0;
    if (pnew) k=desngb; else k=0;
    for(j=k; j<nodes1[th].count; j++)
    {
	i=nodes1[th].lid+j;
	if(imarkA[i]) continue;
	for(k=0,r=0; k<ND; k++)
	{
#ifdef PERIODIC
	    temp=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k]);
#else
	    temp=(Part[i].Pos[k]-searchcenter[k]);
#endif

	    r+=temp*temp;
	}

	if(r < pqHead->r)
	{
	    imarkA[pqHead->p]=0; imarkA[i]=1;
	    pqHead->r=r;	pqHead->p=i; 
	    for(k=0;k<ND;k++)
		pqHead->x[k]=Part[i].Pos[k];
	    adjust_heap(pqxA,desngb); b1=1;
	}

    }
    if(b1)
	for(k=0; k<ND; k++)
	{
	    temp=sqrt(pqHead->r);    
	    searchx[k][0]=searchcenter[k]-temp;
	    searchx[k][1]=searchcenter[k]+temp;
	}


    ngb_treesearch_sphere_nometric(th,nodes1,desngb/2,desngb,pqxA,imarkA);



//----------------------------------------


//----------------------------------------



    pnext=pqHead->p;
    r=pqHead->r;

//     vector<float> rt;
//     for(i=0,rhoxyz=0; i<desngb; i++)
// 	rt.push_back(sqrt(pqStartA[i].r/pqHead->r));
//     sort(rt.begin(),rt.end());


    for(i=0,rhoxyz=0; i<desngb; i++)
    {

	u = sqrt(pqStartA[i].r/pqHead->r);
//	u=rt[i];
	if(u<1)
	{
	    if ((u==0) && (All.KernelBiasCorrection==1))
		u=ND/(1.0+ND);
	    k = (int)(u*KERNEL_TABLE);
	    wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);  rhoxyz+=wk*(Part[pqStartA[i].p].Mass);
	}

//	    cout<<i<<" b "<<pqStartA[i].p<<" a "<<u<<endl;

	if(idoneA[pqStartA[i].p]==0)
	    if(pqStartA[i].r<r)
	    {pnext=pqStartA[i].p;r=pqStartA[i].r;}

    }
    if(idoneA[pnext]==1) pnew=1; else pnew=0;

    temp=sqrt(pqHead->r);
    for(k=0,r=1.0,wk=1.0;k<ND;k++)
    {
//	 r*=metric[k];
	wk*=temp;
    }
    temp=r/wk;
//    cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;
//    cout<<rhoxyz<<" "<<All.hsv<<" "<<r/wk<<endl;


    return rhoxyz;
}




float  ngb_treedensity_sphere_metric_exact_h(float xyz[ND],struct NODE *nodes1, int desngb, struct linklist* pqStartA,bool *imarkA)
{
    real u,wk,r;  /* search radius */
    int   i,j,k,th;
    float rhoxyz;
    double temp,temp1;
    struct linklist *pq;
    vector<class rmlist> rmlistv;
    class rmlist rmliste;

    /* calculate metric */

    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];


//    rmlistv.clear();

    for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
    {
	for(j=0,temp1=0.0; j<ND; j++)
	{

#ifdef PERIODIC
	    temp=periodic((pq->x[j]-searchcenter[j]),All.boxh[j])*metric[j];
#else
	    temp=(pq->x[j]-searchcenter[j])*metric[j];
#endif
	    temp1+=temp*temp;
	}
	imarkA[pq->p]=1;
	rmliste.r=temp1;
	rmliste.Mass=Part[pq->p].Mass;
	rmlistv.push_back(rmliste);
    }


    /* identify the node containing the cell*/
    th=1;
    while(nodes1[th].count>desngb)
    {
	j=nodes1[th].k1;
	if(xyz[j] < nodes1[(LOWER(th,nodes1))].bnd->x[j][1])
	    th=(LOWER(th,nodes1));
	else
	    th=(UPPER(th,nodes1));
    }
    while(nodes1[th].count<desngb)
    {
	th=PARENT(th,nodes1);
    }



    for(k=0; k<ND; k++)
    {
	temp=1/metric[k];    
	searchx[k][0]=searchcenter[k]-temp;
	searchx[k][1]=searchcenter[k]+temp;
    }


    k=0;
    for(j=k; j<nodes1[th].count; j++)
    {
	i=nodes1[th].lid+j;
	if(imarkA[i]) continue;
	for(k=0,r=0; k<ND; k++)
	{
#ifdef PERIODIC
	    temp=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
	    temp=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif
	    r+=temp*temp;
	}
	if(r < 1.0)
	{
	    rmliste.r=r;rmliste.Mass=Part[i].Mass;
	    rmlistv.push_back(rmliste);
	}
    }


    ngb_treesearch_sphere_metric_exact_h(th,nodes1,desngb/2,rmlistv,imarkA);



//----------------------------------------


//----------------------------------------



    for(i=0,rhoxyz=0; i<int(rmlistv.size()); i++)
    {
	u = sqrt(rmlistv[i].r);
	if(u<1)
	{
	    if ((u==0) && (All.KernelBiasCorrection==1))
		u=ND/(1.0+ND);
	    k = (int)(u*KERNEL_TABLE);
	    wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);  rhoxyz+=wk*rmlistv[i].Mass;
	}

    }

    for(k=0,r=1.0,wk=1.0;k<ND;k++)
    {
	r*=metric[k];
    }
    temp=r/wk;
//    cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;
//    cout<<rhoxyz<<" "<<All.hsv<<" "<<r/wk<<endl;


    for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
	imarkA[pq->p]=0;


    return rhoxyz;
}





float  ngb_treedensity_box_metric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA)
{
    real u,wk,r;  /* search radius */
    int   i,j,k,th,l;
    float rhoxyz;
    real xmin[ND],xmax[ND],temp,temp1;
    struct linklist *pq;
    bool b1,b2=1;
    int b3=1;
    struct pqueue ptemp;
    int hj,hi;

    /* calculate metric */
    for(j=0; j<ND; j++)
	metric[j]=1.0;

    if(All.TypeOfSmoothing==5)
    {
	if(All.CubicCells==1)
	{
	    if(ND==6)
	    {


// 	    th=1;
// 	    while(nodes1[th].count>1)
// 	    {
// 		j=nodes1[th].k1;
// 		if(xyz[j] < nodes1[LOWER(th,nodes1)].bnd->x[j][1])
// 		    th=LOWER(th,nodes1);
// 		else
// 		    th=UPPER(th,nodes1);
// 	    }
// 	    for(j=0; j<ND; j++)
// 	    {
// 		metric[j]=(nodes1[th].bnd->x[j][1]-nodes1[th].bnd->x[j][0]);
// 	    }


		for(j=0; j<ND; j++)
		    searchcenter[j]=0.0;
		ngb_calculate_metric(xyz,nodes1);
		for(j=0; j<ND; j++)
		{
		    metric[j]=searchcenter[j];
		}


		for(j=0,temp=1.0; j<3; j++)
		    temp*=(metric[j]);
		temp=pow(temp,0.3333);
		for(j=0; j<3; j++)
		    metric[j]=1/temp;

		for(j=3,temp=1.0; j<6; j++)
		    temp*=(metric[j]);
		temp=pow(temp,0.3333);
		for(j=3; j<6; j++)
		    metric[j]=1/temp;
	    }
	}
	else
	{
	    for(j=0; j<ND; j++)
		searchcenter[j]=0.0;
	    ngb_calculate_metric(xyz,nodes1);
	    for(j=0; j<ND; j++)
	    {
		metric[j]=1.0/searchcenter[j];
	    }
	}
    }


    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];


    while(b2)
    {
	/* identify the node containing the cell*/
	th=1;
	while(nodes1[th].count>desngb)
	{
	    j=nodes1[th].k1;
	    if(xyz[j] < nodes1[(LOWER(th,nodes1))].bnd->x[j][1])
		th=(LOWER(th,nodes1));
	    else
		th=(UPPER(th,nodes1));
	}
	while(nodes1[th].count<desngb)
	{
	    th=PARENT(th,nodes1);
	}


	if(pnew)
	{
	    for(pq=pqStartA;pq<(pqStartA+desngb);pq++)
		imarkA[pq->p]=0;

	    for(pq=pqStartA,j=0;pq<(pqStartA+desngb);j++,pq++)
	    {
		i=nodes1[th].lid+j;
		pq->p=i;
		for(k=0; k<ND; k++)
		    pq->x[k]=Part[i].Pos[k];
		imarkA[i]=1;
	    }
	}

	for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
	{
	    for(j=0,pq->r=0.0; j<ND; j++)
	    {
#ifdef PERIODIC
		temp=periodic((pq->x[j]-searchcenter[j]),All.boxh[j])*metric[j];
#else
		temp=(pq->x[j]-searchcenter[j])*metric[j];
#endif
		if(temp<0) temp=-temp;
		if(pq->r<temp) pq->r=temp;
	    }
	    pqxA[desngb-i-1].pq=pq;
	}


	make_heap(pqxA,(pqxA+desngb));      pqHead=pqxA->pq;

	for(k=0; k<ND; k++)
	{
	    temp=pqHead->r/metric[k];    
	    searchx[k][0]=searchcenter[k]-temp;
	    searchx[k][1]=searchcenter[k]+temp;
	}


	b1=0;
	if (pnew) k=desngb; else k=0;
	for(j=0; j<nodes1[th].count; j++)
	{
	    i=nodes1[th].lid+j;
	    if(imarkA[i]) continue;
	    for(k=0,r=0; k<ND; k++)
	    {
#ifdef PERIODIC
		temp=periodic((Part[i].Pos[k]-searchcenter[k]),All.boxh[k])*metric[k];
#else
		temp=(Part[i].Pos[k]-searchcenter[k])*metric[k];
#endif

		if(temp<0) temp=-temp;
		if(r<temp) r=temp;
	    }

	    if(r < pqHead->r)
	    {
		imarkA[pqHead->p]=0; imarkA[i]=1;
		pqHead->r=r;	pqHead->p=i; 
		for(k=0;k<ND;k++)
		    pqHead->x[k]=Part[i].Pos[k];
		adjust_heap(pqxA,desngb); b1=1;
	    }

	}
	if(b1)
	    for(k=0; k<ND; k++)
	    {
		temp=pqHead->r/metric[k];    
		searchx[k][0]=searchcenter[k]-temp;
		searchx[k][1]=searchcenter[k]+temp;
	    }




	ngb_treesearch_box_metric(th,nodes1,desngb/2,desngb,pqxA,imarkA);

//----------------------------------------


	if(All.TypeOfSmoothing==4) b3=2;
	b2=0;
	if(b3<2)
	    if(All.VolCorr==1)
	    {

#ifdef PERIODIC

		for(i=0; i<desngb; i++)
		{

		    for(k=0; k<ND; k++)
		    {
			temp=periodic((pqStartA[i].x[k]-searchcenter[k]),All.boxh[k]);
			if((i==0)||(temp<xmin[k])) xmin[k]=temp;
			if((i==0)||(temp>xmax[k])) xmax[k]=temp;
		    }
		}

		    for(k=0; k<ND; k++)
		    {
			xmin[k]+=searchcenter[k];
			xmax[k]+=searchcenter[k];
		    
		    }
#else
		for(i=0; i<desngb; i++)
		{

		    for(k=0; k<ND; k++)
		    {
			if((i==0)||(pqStartA[i].x[k]<xmin[k])) xmin[k]=pqStartA[i].x[k];
			if((i==0)||(pqStartA[i].x[k]>xmax[k])) xmax[k]=pqStartA[i].x[k];
		    }
		}
#endif


		for(k=0; k<ND; k++)
		{
		    temp=sqrt(pqHead->r)/metric[k];    
		    searchx[k][0]=searchcenter[k]-temp;
		    searchx[k][1]=searchcenter[k]+temp;
		}

		for(k=0;k<ND;k++)
		{

		    if(All.CubicCells==1) 		    
			temp1=(xmax[k]-xmin[k])*15.0/(All.DesNumNgb);
		    else
			temp1=(xmax[k]-xmin[k])*15.0/(All.DesNumNgb);

		    if(((searchx[k][1]-xmax[k])>temp1)&&((xmin[k]-searchx[k][0])>temp1))
		    {			
			b2=1;
		    }

		}

		if(b2)
		    for(k=0;k<ND;k++)
		    {
// 			temp=(xmax[k]-xmin[k])*0.5;
// 			metric[k]=1/temp;
			temp=(xmax[k]-xmin[k])*0.5;
			if(xmax[k]!=xmin[k])
			{			
			    metric[k]=1/temp;
			}
			else
			    metric[k]=1/(searchx[k][1]-searchx[k][0]);


		    }

//   	    if(b3>5) 
//  	    cout<<b3<<endl;
		b3++;
	    }
    }

//----------------------------------------




    pnext=pqHead->p;
    r=pqHead->r;


    for(i=0,rhoxyz=0; i<desngb; i++)
    {

	for(k=0,r=0; k<ND; k++)
	{
	    xmax[k]=pqHead->r/metric[k];
#ifdef PERIODIC
	    xmin[k]=periodic((pqStartA[i].x[k] - searchcenter[k]),All.boxh[k]);
#else
	    xmin[k]=pqStartA[i].x[k] - searchcenter[k];
#endif
	    temp=(xmin[k])/xmax[k];
	    if(temp<0) temp=-temp;
	    xmin[k]=temp;	       
	    r+=temp;
	}

	for(l=0,wk=1.0;l<ND;l++)
	{
	    u=xmin[l];

	    if ((r==0.0) && (All.KernelBiasCorrection==1))
		u=0.5;
	    if(u<1)
	    {
		k = (int)(u*KERNEL_TABLE);
		wk*=( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);

	    }
	    else
		wk=0.0;
	}



	rhoxyz+=Part[pqStartA[i].p].Mass*wk;



	if(idoneA[pqStartA[i].p]==0)
	    if(pqStartA[i].r<r)
	    {pnext=pqStartA[i].p;r=pqStartA[i].r;}



    }
    if(idoneA[pnext]==1) pnew=1; else pnew=0;


    temp=pqHead->r;
    for(k=0,r=1.0,wk=1.0;k<ND;k++)
    {
	r*=metric[k];
	wk*=temp;
    }
    temp=r/wk;
//    cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;


//     cout<<"density"<<rhoxyz<<endl;
//     endrun(10);

    return rhoxyz;
}








float  ngb_treedensity_sphere_gmatrix(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA)
{
    real u,wk,r;  /* search radius */
    int   i,j,k,th,l;
    float rhoxyz;
    double temp;
    struct linklist *pq;
    bool b1;


    struct pqueue ptemp;
    int hj,hi;

    /* calculate metric */

    for(j=0; j<ND; j++)
	searchcenter[j]=xyz[j];


    /* identify the node containing the cell*/
    th=1;
    while(nodes1[th].count>desngb)
    {
	j=nodes1[th].k1;
	if(xyz[j] < nodes1[(LOWER(th,nodes1))].bnd->x[j][1])
	    th=(LOWER(th,nodes1));
	else
	    th=(UPPER(th,nodes1));
    }
    while(nodes1[th].count<desngb)
    {
	th=PARENT(th,nodes1);
    }

//     pnew=1;
    if(pnew)
    {
	for(pq=pqStartA,j=0;pq<(pqStartA+desngb);j++,pq++)
	{
	    i=nodes1[th].lid+j;
	    if(imarkA[i]==1) continue;
	    imarkA[pq->p]=0;
	    pq->p=i;
	    for(k=0; k<ND; k++)
		pq->x[k]=Part[i].Pos[k];
	    imarkA[i]=1;
	}
    }

    for(pq=pqStartA,i=0;pq<(pqStartA+desngb);++i,pq++)
    {
	for(j=0,pq->r=0.0; j<ND; j++)
	{
	    for(l=0,temp=0; l<ND; ++l)
	    {
#ifdef PERIODIC
		temp+=periodic((pq->x[l]-searchcenter[l]),All.boxh[l])*gmatrix[l][j]*metric1[l];
#else
		temp+=(pq->x[l]-searchcenter[l])*gmatrix[l][j]*metric1[l];
#endif
	    }

	    temp=temp*metric[j];
	    pq->r+=temp*temp;
	}
	pqxA[desngb-i-1].pq=pq;
    }


    make_heap(pqxA,(pqxA+desngb));      pqHead=pqxA->pq;

    for(k=0; k<ND; k++)
    {
	temp=sqrt(pqHead->r)/metric1[k];    
	searchx[k][0]=searchcenter[k]-temp;
	searchx[k][1]=searchcenter[k]+temp;
    }


    b1=0;
    if (pnew) k=desngb; else k=0;
    for(j=0; j<nodes1[th].count; j++)
    {
	i=nodes1[th].lid+j;
	if(imarkA[i]) continue;
	for(k=0,r=0; k<ND; k++)
	{
	    for(l=0,temp=0; l<ND; ++l)
	    {
#ifdef PERIODIC
		temp+=periodic((Part[i].Pos[l]-searchcenter[l]),All.boxh[l])*gmatrix[l][k]*metric1[l];
#else
		temp+=(Part[i].Pos[l]-searchcenter[l])*gmatrix[l][k]*metric1[l];
#endif
	    }
	    temp=temp*metric[k];
	    r+=temp*temp;

	}

	if(r < pqHead->r)
	{
	    imarkA[pqHead->p]=0; imarkA[i]=1;
	    pqHead->r=r;	pqHead->p=i; 
	    for(k=0;k<ND;k++)
		pqHead->x[k]=Part[i].Pos[k];
	    adjust_heap(pqxA,desngb); b1=1;

	}

    }
    if(b1)
	for(k=0; k<ND; k++)
	{
	    temp=sqrt(pqHead->r)/metric1[k];    
	    searchx[k][0]=searchcenter[k]-temp;
	    searchx[k][1]=searchcenter[k]+temp;
	}

//	cout<<" c "<<pqHead->r<<endl;


    ngb_treesearch_sphere_gmatrix(th,nodes1,desngb/2,desngb,pqxA,imarkA);
//----------------------------------------

    pnext=pqHead->p;
    r=pqHead->r;


    for(k=0,temp=metric[0];k<ND;k++)
	if(temp>metric[k]) temp=metric[k];


//    int ii=0;
    for(i=0,rhoxyz=0; i<desngb; i++)
    {


	u = sqrt(pqStartA[i].r/pqHead->r);


	if(u<1)
	{
	    if ((u==0) && (All.KernelBiasCorrection==1))
		u=ND/(1.0+ND);
	    k = (int)(u*KERNEL_TABLE);
	    wk =( Kernel[k]  + (Kernel[k+1]-Kernel[k])*(u-KernelRad[k])*KERNEL_TABLE);      rhoxyz+=wk*Part[pqStartA[i].p].Mass;
	}

//	    cout<<i<<" "<<pqStartA[i].p<<" a "<<u<<endl;


	if(idoneA[pqStartA[i].p]==0)
	    if(pqStartA[i].r<r)
	    {pnext=pqStartA[i].p;r=pqStartA[i].r;}
//	cout<<pqStartA[i].p<<" b "<<u<<endl;

    }
    if(idoneA[pnext]==1) pnew=1; else pnew=0;


    temp=sqrt(pqHead->r);
    for(k=0,r=1.0,wk=1.0;k<ND;k++)
    {
	r*=metric[k]*metric1[k];
	wk*=temp;
    }
    temp=r/wk;
//    cout<<rhoxyz<<" "<<r<<" "<<wk<<endl;
    rhoxyz*=temp;
//    cout<<rhoxyz<<endl;


    return rhoxyz;
}







float anisokernel_density(float xyz[ND],struct NODE *nodes1, int desngbA, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA, int desngbB, struct pqueue* pqxB, struct linklist* pqStartB, bool *idoneB,bool *imarkB)
{
    int n,k,l,nrot;
    double dx[ND],dxcm[ND],mtot=0.0,ht1,rcm;
    float rhoxyz;

    for(k=0;k<ND;k++)
	dxcm[k]=0.0;

    if(All.TypeOfSmoothing==2)
    {
	rhoxyz=ngb_treedensity_sphere_nometric(xyz,nodes1,desngbA,pqxA,pqStartA,idoneA,imarkA);
	for(k=0;k<ND;k++)
	    metric1[k]=1.0;
    }
    if(All.TypeOfSmoothing==3)
    {
	rhoxyz=ngb_treedensity_sphere_metric(xyz,nodes1,desngbA,pqxA,pqStartA,idoneA,imarkA);
	for(k=0;k<ND;k++)
	    metric1[k]=metric[k];
    }


//    for(k=0;k<ND;k++)
//         cout<<sqrt(pqHead->r)/metric[k]<<" "<<endl;

    for(k=1;k<=ND;k++)
	for(l=1;l<=ND;l++)
	    mrho[k][l]=0.0;

    // Calculate the center of mass dxcm[k]
    for(n=0; n<desngbA; n++)
    {
	for(k=0;k<ND;k++)
	{
#ifdef PERIODIC
	    dx[k] = periodic((xyz[k] - pqStartA[n].x[k]),All.boxh[k]);
#else
	    dx[k] = xyz[k] - pqStartA[n].x[k];
#endif

	    dxcm[k]+=dx[k]*Part[pqStartA[n].p].Mass;
	}
	mtot +=Part[pqStartA[n].p].Mass;
    }

    for(k=0,rcm=0.0;k<ND;k++)
    {
	dxcm[k]=dxcm[k]/mtot;
	rcm+=dxcm[k]*dxcm[k];
    }
    rcm=sqrt(rcm);




    // Calculate the covariance matrix
    for(n=0; n<desngbA; n++)
    {
	for(k=0;k<ND;k++)
	{
#ifdef PERIODIC
	    dx[k] = periodic((xyz[k] - pqStartA[n].x[k]),All.boxh[k]);
#else
	    dx[k] = xyz[k] - pqStartA[n].x[k];
#endif
	    dx[k] =dx[k]-dxcm[k];
	}

	for(k=1;k<=ND;k++)
	    for(l=1;l<=ND;l++)
	    {
		mrho[k][l]+=dx[k-1]*dx[l-1]*Part[pqStartA[n].p].Mass;
	    }
    }


    for(k=1;k<=ND;k++)
	for(l=1;l<=ND;l++)
	{
	    mrho[k][l]=mrho[k][l]*metric1[k-1]*metric1[l-1]/(mtot*sqrt(pqHead->r*pqHead->r));
//	    cout<<k<<" "<<l<<" "<<mrho[k][l]<<endl;

	}


    // Calculate the transoformation matrix "ve" of unit covariance and update scaling "d"
    jacobi(mrho,ND,d,ve,&nrot);

    for(k=1;k<=ND;k++)
    {
	d[k]=fabs(sqrt(d[k]));
    }


    ht1=a_max(d,1,ND);
    for(k=1;k<=ND;k++)
    {
	d[k]=d[k]/ht1;
    }


//     for(k=1;k<=ND;k++)
//       for(l=1;l<=ND;l++)
// 	    cout<<k<<" "<<l<<" ve "<<ve[k][l]<<endl;



    if(a_min(d,1,ND)<All.Anisotropy)
    {
	for(k=1;k<=ND;k++)
	{
	    if(d[k]<All.Anisotropy)
		d[k]=All.Anisotropy;
	}
	for(k=1;k<=ND;k++)
	    for(l=1;l<=ND;l++)
	    {
		ve[k][l]=0.0;
		if(k==l) ve[k][l]=1.0;
	    }
    }


    for(k=1;k<=ND;k++)
	for(l=1;l<=ND;l++)
	{
	    gmatrix[k-1][l-1]=ve[k][l];
//  	    gmatrix[k-1][l-1]=0.0;
//  	    if(k==l)gmatrix[k-1][l-1]=1.0;
	}

    for(k=1,ht1=1.0/d[1];k<=ND;k++)
    {
	metric[k-1]=1.0/d[k];
	if(ht1>metric[k-1])
	    ht1=metric[k-1];
    }

    for(k=0;k<ND;k++)
    {
	metric[k]=metric[k]/ht1;
//	metric[k]=1.0;
//	cout<<" a "<<1/metric[k]<<" "<<metric[k]<<endl;
    }



    rhoxyz=ngb_treedensity_sphere_gmatrix(xyz,nodes1,desngbB,pqxB,pqStartB,idoneB,imarkB);



//    cout<<rhoxyz<<endl;
    return rhoxyz;
}


void create_linklist(struct pqueue* &pqx,struct linklist* &pqStart,bool* &imark,bool* &idone,int desngb,int numpart)
{
    int i;
    pqStart=new linklist[desngb];
    pqx=new pqueue[desngb];
    for(i=0;i<desngb;++i)
	pqx[i].pq=pqStart+i;
    imark= new bool[numpart];
    idone= new bool[numpart];
    for(i=0; i<numpart; i++)
    {idone[i]=0; }
}



void initialize_linklist(struct linklist* pqStart,bool* imark,int iStart,int desngb,int numpart)
{
    int k,ji;
    struct linklist *pq;
    for(k=0; k<numpart; k++)
	imark[k]=0;
//     for(j=0; j<parttype; j++)
// 	i+=npart[j];

    for(pq=pqStart,ji=iStart;pq<(pqStart+desngb);ji++,pq++)
    {
//		cout<<i<<" "<<ji<<" "<<nodes1<<" "<<trees[j]<<endl;
	pq->p=ji;
	for(k=0; k<ND; k++)
	    pq->x[k]=Part[ji].Pos[k];
	imark[ji]=1;
    }
    
}


float density_general(float xyz[ND],struct NODE *nodes1, int desngbA, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA, int desngbB, struct pqueue* pqxB, struct linklist* pqStartB, bool *idoneB,bool *imarkB)
{
    float Density=0;
    if(All.AnisotropicKernel==1)
    {
	Density= anisokernel_density(xyz, nodes1,desngbA,pqxA,pqStartA,idoneA,imarkA,desngbB,pqxB,pqStartB,idoneB,imarkB)/All.hsv;
    }
    else
    {
	switch (All.TypeOfSmoothing)
	{
	    case 0:
		Density=ngb_treedensity_raw(xyz,nodes1)/All.hsv;
		break;
	    case 1:
		Density=ngb_treedensity_fiestas(xyz,nodes1)/All.hsv;
		break;
	    case 2:
		Density= ngb_treedensity_sphere_nometric(xyz, nodes1,desngbB,pqxB,pqStartB,idoneB,imarkB)/All.hsv;
		break;
	    case 3:
		Density= ngb_treedensity_sphere_metric(xyz, nodes1,desngbB,pqxB,pqStartB,idoneB,imarkB)/All.hsv;
		break;
	    case 4:
		Density= ngb_treedensity_box_metric(xyz, nodes1,desngbB,pqxB,pqStartB,idoneB,imarkB)/All.hsv;
		break;
	    case 5:
		Density= ngb_treedensity_box_metric(xyz, nodes1,desngbB,pqxB,pqStartB,idoneB,imarkB)/All.hsv;
		break;
	    default:
		cout<<"Usage: TypeOfSmoothing  0 to 6"<<endl;
		endrun(10);
		break;
	}
    }
    return Density;
}




void print_list(struct linklist *pq1,int size1)
{
    struct linklist *pq;
    int i;
    cout<<" List "<<endl;
    pq=pq1;
    for(i=0;i<size1;i++)
    {
	cout<<i<<" "<<pq<<"  "<<"  "<<pq->r<<" "<<pq->p<<endl;
	if(i!=(size1-1))
	    pq++;
    }
    cout<<" complete "<<endl;
}



