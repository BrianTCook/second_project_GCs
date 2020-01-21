/***************************************************************************
                          tree_search.cpp  -  description
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
void ngb_treesearch_sphere_metric(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool* imarkA);
void ngb_treesearch_sphere_metric_periodic(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool* imarkA);
void ngb_treesearch_sphere_gmatrix(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool * imarkA);
void ngb_treesearch_sphere_nometric(int no, struct NODE* nodes1,int NumBucket, int DesNumNgb,struct pqueue* pqxA,bool * imarkA);
void ngb_treesearch_box_metric(int no,struct NODE *nodes1, int NumBucket, int DesNumNgb,struct pqueue* pqxA, bool * imarkA);
float  ngb_treedensity_sphere_metric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA);
float  ngb_treedensity_sphere_metric_periodic(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA);
float  ngb_treedensity_sphere_nometric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA);
float  ngb_treedensity_box_metric(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA);
float  ngb_treedensity_sphere_gmatrix(float xyz[ND],struct NODE *nodes1, int desngb, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA);
float anisokernel_density(float xyz[ND],struct NODE *nodes1, int desngbA, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA, int desngbB, struct pqueue* pqxB, struct linklist* pqStartB, bool *idoneB,bool *imarkB);
void create_linklist(struct pqueue* &pqx,struct linklist* &pqStart,bool* &imark,bool* &idone,int desngb,int numpart);
void initialize_linklist(struct linklist* pqStart,bool* imark,int iStart,int desngb,int numpart);
float density_general(float xyz[ND],struct NODE *nodes1, int desngbA, struct pqueue* pqxA, struct linklist* pqStartA, bool *idoneA,bool *imarkA, int desngbB, struct pqueue* pqxB, struct linklist* pqStartB, bool *idoneB,bool *imarkB);
void print_list(struct linklist *pq1,int size1);
float  ngb_treedensity_bruteforce(float xyz[ND], int parttype,int desngb, struct pqueue* pqxA, struct linklist* pqStartA);



