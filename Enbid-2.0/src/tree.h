/***************************************************************************
                          tree.h  -  description
                             -------------------
    begin                : Mon Jan 16 2006
    copyright            : (C) 2006 by Sanjib Sharma
    email                : sharma@wesleyan.edu
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
// tree building routines
void   treeallocate(int maxpart);
void   treefree(void);
int    treebuild_all(void);
int    treebuild_single(int startnode, vector<int> &npart, int type,int *creatednodes);
void   treebuild(int);
void treeprint(int no);
void sort2(unsigned long n, float arr[], int brr[]);
//* search in a box routines
void   ngb_treesearch_kdts(int no);
void   ngb_treesearch_kdts_vc(int no);
void   ngb_treesearch_kdts_h(int no);
float ngb_treedensity_fiestas(float xyz[ND], struct NODE* nodes1);
float ngb_treedensity_raw(float xyz[ND], struct NODE* nodes1);
void ngb_calculate_metric(float xyz[ND], struct NODE* nodes1);



