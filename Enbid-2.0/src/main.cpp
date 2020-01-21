/***************************************************************************
                          main.cpp  -  description
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
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "tree.h"
#include "functions.h"

/* Calls begrun() */
int main(int argc, char **argv)
{
    if(argc!=2)
    {
	fprintf(stdout,"Parameter file  is missing. Usage:\n");
	cout<<"./enbid <ParameterFile> [<begin_no> <end_no>] "<<endl;
	exit(1);
    }
    strcpy(ParameterFile,argv[1]);
  
	t00=0.0;
	begrun();     /* set-up run  */
	cout<<"\nTotal Time = "<<t00<<" s \n"<<endl;
	treefree();
	free_memory();


    return 0;
}









// int main(int argc, char **argv)
// {
//     if(argc!=2)
//     {
// 	fprintf(stdout,"Parameter file  is missing. Usage:\n");
// 	cout<<"./enbid <ParameterFile> [<begin_no> <end_no>] "<<endl;
// 	exit(1);
//     }
//     strcpy(ParameterFile,argv[1]);
  
//     Var.BeginFileNum=0;
//     Var.EndFileNum=0;
//     BatchFlag=0;
  
//     if(argc>2)
//     {
// 	BatchFlag=1;
// 	if(argc==4)
// 	{
// 	    Var.BeginFileNum=atoi(argv[2]);
// 	    Var.EndFileNum=atoi(argv[3]);
// 	    if(Var.BeginFileNum>Var.EndFileNum)
// 	    {
// 		cout<<"Final file no should be greater than inital"<<endl;
// 		endrun(10);
// 	    }
// 	}
// 	else
// 	{
// 	    cout<<"File numbers missing. Usage:"<<endl;
// 	    cout<<endl;
// 	    cout<<"./enbid <ParameterFile> [<begin_no> <end_no>] "<<endl;
// 	    cout<<endl;
// 	    endrun(10);
// 	}
//     }
  
//     for(Var.FileNum=Var.EndFileNum; Var.FileNum >= Var.BeginFileNum;Var.FileNum--)
//     {
// 	t00=0.0;
// 	begrun();     /* set-up run  */
// 	cout<<"\nTotal Time = "<<t00<<" s \n"<<endl;
// 	treefree();
// 	free_memory();
//     }
  


//     return 0;
// }




