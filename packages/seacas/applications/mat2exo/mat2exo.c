/*
 * Copyright(C) 2012 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *                                                 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/**********************************************************************/
/* matlab mat file to exodus II.  This takes .mat files exactly as
   generated by the tool exo2mat written by mrtabba and converts them
   back to exodus II format

   rmnaeth. August 8, 2003

   modified by D. Todd Griffith on 12/09/2005
   * modifcations include:
   1) writes global, nodal and element variable names
   2) writes global, nodal and elemnent variable results
   3) writes complete set of time steps (previous version
          skipped first step) 
   4) writes complete node set information (node set numbers, 
          dist. factors, etc)
   5) writes complete side set information (side set numbers, 
          dist. factors, etc)

   modified by D. Todd Griffith on 12/16/2005
   * side set distribution factors now written as double (not int)

   modified by Greg Sjaardema, 07/05/2012 to use matio instead of matlab libraries.
*/

#include <exodusII.h>                   // for ex_inquire_int, ex_put_var, etc
#include <stddef.h>                     // for size_t
#include <stdio.h>                      // for sprintf, NULL, printf, etc
#include <stdlib.h>                     // for calloc, free, exit
#include <string.h>                     // for strtok, memcpy, strcat, etc
#include "add_to_log.h"                 // for add_to_log
#include "matio.h"                      // for matvar_t, Mat_VarFree, etc

/**********************************************************************/
mat_t *mat_file=0;  /* file for binary .mat input */

/**********************************************************************/
static char *qainfo[] =
{
  "mat2exo",
  "2014/01/14",
  "2.02",
};

/**********************************************************************/
int matGetStr  (char *name,char *str);
int matGetDbl  (char *name,int n1,int n2, double *data);
int matGetInt  (char *name,int n1,int n2, int *data);
int matArrNRow (char *name);
int matArrNCol (char *name);
void del_arg(int *argc, char* argv[], int j);

/**********************************************************************/
int main (int argc, char *argv[]){

  char **str2,*line,*curr;
    
  const char* ext=".exo";

  int   
    i,j,k,n,n1,cpu_word_size,io_word_size,exo_file,
    num_axes,num_nodes,num_elements,num_blocks,
    num_side_sets,num_node_sets,num_time_steps,
    num_global_vars,
    num_nodal_vars,num_element_vars,*ids,*iscr,
    *nsssides,*nssdfac,*elem_list,*side_list,
    *nnsnodes,*nnsdfac,*node_list;

  double
    *scr,*x,*y,*z,
    *escr;

  char * blknames = NULL;
  int *num_elem_in_block = NULL;

  /* QA Info */
  printf("%s: %s, %s\n", qainfo[0], qainfo[2], qainfo[1]);

  /* usage message*/
  if(argc != 2){
    printf("%s matlab_file_name.\n",argv[0]);
    printf("   the matlab_file_name is required\n");
    printf("%d", argc);
    exit(1);
  }
  
  /*open input file*/
  mat_file = Mat_Open(argv[1], MAT_ACC_RDONLY);
  if (mat_file == NULL) {
    printf("Error opening matlab file %s\n", argv[1]);
    return(1);
  }

  /*open output file*/
  cpu_word_size=sizeof(double);
  io_word_size=sizeof(double);
  /* QA records */
  ext=".exo";
  line = (char *) calloc (2049,sizeof(char));
  strcpy(line,argv[1]);
  strtok(line,".");  
  strcat(line,ext);
  exo_file = ex_create(line,EX_CLOBBER,&cpu_word_size,&io_word_size);
  if (exo_file < 0){
    printf("error creating %s\n",line);
    exit(1);
  }

  /* print */
  fprintf(stderr,"translating %s to %s ... ",argv[1],line);
  
  /* read database parameters */
  matGetInt("naxes",  1, 1,&num_axes);
  matGetInt("nnodes", 1, 1,&num_nodes);
  matGetInt("nelems", 1, 1,&num_elements);
  matGetInt("nblks",  1, 1,&num_blocks);
  matGetInt("nnsets", 1, 1,&num_node_sets);
  matGetInt("nssets", 1, 1,&num_side_sets);
  matGetInt("nsteps", 1, 1,&num_time_steps);
  matGetInt("ngvars", 1, 1,&num_global_vars);
  matGetInt("nnvars", 1, 1,&num_nodal_vars);
  matGetInt("nevars", 1, 1,&num_element_vars);

  /*export parameters */
  ex_put_init(exo_file,line,
	      num_axes,num_nodes,num_elements,num_blocks,
	      num_node_sets,num_side_sets);
  free(line);
  
  if ( num_global_vars > 0 ){
    ex_put_variable_param(exo_file,EX_GLOBAL,num_global_vars);
  }
  
  if ( num_nodal_vars > 0 ){
    ex_put_variable_param(exo_file,EX_NODAL,num_nodal_vars);
  }
  
  if ( num_element_vars > 0 ){
    ex_put_variable_param(exo_file,EX_ELEM_BLOCK,num_element_vars);
  }

  /* nodal coordinates */
  x = (double *) calloc(num_nodes,sizeof(double));
  y = (double *) calloc(num_nodes,sizeof(double));
  if (num_axes == 3) 
    z = (double *) calloc(num_nodes,sizeof(double));
  else 
    z = NULL;
  matGetDbl("x0", num_nodes, 1, x);
  matGetDbl("y0", num_nodes, 1, y);
  if (num_axes == 3)
    matGetDbl("z0", num_nodes,1,z);
  ex_put_coord(exo_file,x,y,z);
  free(x);
  free(y);
  if (num_axes == 3){ 

    free(z);
  }
  

  /* side sets (section by dgriffi) */
  if(num_side_sets > 0){ 
     
    /* ssids */
    ids = (int *) calloc(num_side_sets,sizeof(int)); 
    matGetInt("ssids",num_side_sets, 1,ids);

    /* nsssides */
    nsssides = (int *) calloc(num_side_sets,sizeof(int));
    matGetInt("nsssides",num_side_sets,1,nsssides);

    /* nssdfac */
    nssdfac = (int *) calloc(num_side_sets,sizeof(int));
    matGetInt("nssdfac",num_side_sets,1,nssdfac);

    for(i=0;i<num_side_sets;i++){
      char name[32];
  
      ex_put_set_param(exo_file,EX_SIDE_SET,ids[i],nsssides[i],nssdfac[i]);
      elem_list = (int *) calloc(nsssides[i],sizeof(int));
      side_list = (int *) calloc(nsssides[i],sizeof(int));
      escr = (double *) calloc(nssdfac[i],sizeof(double));
           
      sprintf(name,"sselem%02d",i+1);
      matGetInt(name,nsssides[i],1,elem_list);

      sprintf(name,"ssside%02d",i+1);
      matGetInt(name,nsssides[i],1,side_list);
      ex_put_set(exo_file,EX_SIDE_SET,ids[i],elem_list,side_list);

      free(elem_list);
      free(side_list);
      sprintf(name,"ssfac%02d",i+1);
      matGetDbl(name,nssdfac[i],1,escr);
      ex_put_set_dist_fact(exo_file,EX_SIDE_SET,ids[i],escr);
      free(escr);      
    }
   
    free(nsssides);
    free(nssdfac);
    free(ids);
  }  

  /* node sets (section by dgriffi) */
  if(num_node_sets > 0){ 
     
    /* nsids */
    ids = (int *) calloc(num_node_sets,sizeof(int)); 
    matGetInt("nsids",num_node_sets, 1,ids);

    /* nnsnodes */
    nnsnodes = (int *) calloc(num_node_sets,sizeof(int));
    matGetInt("nnsnodes",num_node_sets,1,nnsnodes);

    /* nnsdfac */
    nnsdfac = (int *) calloc(num_node_sets,sizeof(int));
    matGetInt("nnsdfac",num_node_sets,1,nnsdfac);

    for(i=0;i<num_node_sets;i++){
      char name[32];

      ex_put_set_param(exo_file,EX_NODE_SET,ids[i],nnsnodes[i],nnsdfac[i]);
      node_list = (int *) calloc(nnsnodes[i],sizeof(int));
      escr = (double *) calloc(nnsdfac[i],sizeof(double));
           
      sprintf(name,"nsnod%02d",i+1);
      matGetInt(name,nnsnodes[i],1,node_list);
      ex_put_set(exo_file,EX_NODE_SET,ids[i],node_list,NULL);
      free(node_list);
      
      sprintf(name,"nsfac%02d",i+1);
      matGetDbl(name,nnsdfac[i],1,escr);
      ex_put_set_dist_fact(exo_file,EX_NODE_SET,ids[i],escr);
      free(escr);      
    }
   
    free(nnsdfac);
    free(nnsnodes);
    free(ids);
  }  


  /* element blocks */ 
  /* get elem block ids */
  ids = (int *) calloc(num_blocks,sizeof(int));
  matGetInt("blkids",num_blocks,1,ids);

  /* get elem block types */
  blknames = (char *) calloc(num_blocks*(MAX_STR_LENGTH+1),sizeof(char));
  matGetStr("blknames",blknames);
  num_elem_in_block = (int *) calloc(num_blocks,sizeof(int));
  curr = blknames;
  curr = strtok(curr,"\n");
  for(i=0;i<num_blocks;i++){
    char name[32];

    sprintf(name,"blk%02d",i+1);
    n1 = matArrNRow(name);
    n = matArrNCol(name);
    iscr = (int *) calloc(n*n1,sizeof(int));
    matGetInt(name,n1,n,iscr);
    num_elem_in_block[i]=n;
    ex_put_elem_block(exo_file,ids[i],curr,n,n1,0);
    ex_put_conn(exo_file,EX_ELEM_BLOCK,ids[i],iscr,NULL,NULL);
    free(iscr);
    curr = strtok(NULL, "\n");
  }
  free(blknames);

  /* time values */
  if (num_time_steps > 0 ) {
    scr = (double *) calloc(num_time_steps,sizeof(double));
    matGetDbl( "time", num_time_steps, 1,scr);
    for (i=0;i<num_time_steps;i++){
      ex_put_time(exo_file,i+1,&scr[i]);
    }
    free(scr); 
  }
  
  /* global variables */
  if (num_global_vars > 0 ){
    int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    char *str = (char *) calloc(num_global_vars * (max_name_length+1), sizeof(char));
    matGetStr("gnames",str);
    str2 = (char **) calloc(num_global_vars,sizeof(char*));
    curr = strtok(str,"\n");
    for(i=0;i<num_global_vars;i++){
      str2[i]=curr;
      curr = strtok(NULL,"\n");
    }
    ex_put_variable_names(exo_file, EX_GLOBAL, num_global_vars, str2);
    free(str);
    free(str2);

    {
      double * global_var_vals;
      double * temp;
      global_var_vals = (double *) calloc(num_global_vars*num_time_steps,sizeof(double));
      temp            = (double *) calloc(num_time_steps,sizeof(double));
      for (j=0;j<num_global_vars;j++) {
	char name[32];
	sprintf(name,"gvar%02d",j+1);
	matGetDbl(name,num_time_steps,1,temp);
	for (i=0; i < num_time_steps; i++) {
	  global_var_vals[num_global_vars*i+j]=temp[i];
	}
      }
      for (i=0; i<num_time_steps; i++) {
	size_t offset = num_global_vars * i;
	ex_put_var(exo_file,i+1,EX_GLOBAL,1,0,num_global_vars,&global_var_vals[offset]);
      }
      free(temp);
      free(global_var_vals);
    }
  }

  
  /* nodal variables */ /* section by dtg */

  if (num_nodal_vars > 0){
    int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    char *str = (char *) calloc(num_nodal_vars * (max_name_length+1), sizeof(char));
    matGetStr("nnames",str);
    str2 = (char **) calloc(num_nodal_vars,sizeof(char*));
    curr = strtok(str,"\n");
    for(i=0;i<num_nodal_vars;i++){
      str2[i]=curr;
      curr = strtok(NULL,"\n");
    }
    ex_put_variable_names(exo_file, EX_NODAL, num_nodal_vars, str2);	
    free(str);
    free(str2);
    {
      double * nodal_var_vals;
      for (i=0;i<num_nodal_vars;i++) {
	char name[32];
	nodal_var_vals = (double *) calloc(num_nodes*num_time_steps,sizeof(double));
	sprintf(name,"nvar%02d",i+1);
	matGetDbl(name,num_nodes,num_time_steps,nodal_var_vals);
	for (j=0;j<num_time_steps;j++) {
	  ex_put_var(exo_file,j+1,EX_NODAL,i+1,num_nodes,1,nodal_var_vals+num_nodes*j);
	}
	free(nodal_var_vals); 
      }
    }
  }

  /* elemental variables */ /* section by dtg */
  
  if (num_element_vars > 0){
    int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
    char *str = (char *) calloc(num_element_vars * (max_name_length+1), sizeof(char));
    matGetStr("enames",str);
    str2 = (char **) calloc(num_element_vars,sizeof(char*));
    curr = strtok(str,"\n");
    for(i=0;i<num_element_vars;i++){
      str2[i]=curr;
      curr = strtok(NULL,"\n");
    }
    ex_put_variable_names(exo_file, EX_ELEM_BLOCK, num_element_vars, str2);	
    free(str);
    free(str2);
    {
      double * element_var_vals;
      for (i=0;i<num_element_vars;i++) {
	char name[32];
	element_var_vals = (double *) calloc(num_elements*num_time_steps,sizeof(double));       
	sprintf(name,"evar%02d",i+1);
	matGetDbl(name,num_elements,num_time_steps,element_var_vals);
	n=0;       
	for (j=0;j<num_time_steps;j++) {
	  for (k=0;k<num_blocks;k++) {
	    ex_put_var(exo_file,j+1,EX_ELEM_BLOCK, i+1,ids[k],num_elem_in_block[k],element_var_vals+n);
	    n=n+num_elem_in_block[k];
	  }
	}
	free(element_var_vals);
      }
    }
  }
  free(ids); 

  /* node and element number maps */
  ids = (int *) calloc (num_nodes,sizeof(int));
  if ( !matGetInt("node_num_map",num_nodes,1,ids)){
    ex_put_node_num_map(exo_file,ids);
  }
  free(ids);

  ids = (int *) calloc (num_elements,sizeof(int));
  if ( !matGetInt("elem_num_map",num_elements,1,ids)){
    ex_put_elem_num_map(exo_file,ids);
  }
  free(ids);
  free(num_elem_in_block);
  
  /* close exo file */
  ex_close(exo_file);
  
  /* close mat file */
  Mat_Close(mat_file);

  /* */
  fprintf(stderr,"done.\n");

  /* exit status */
  add_to_log("mat2exo", 0);
  return(0);
}

/**********************************************************************/
int matGetStr (char *name,char *data)
{
  int strlen;

  matvar_t *matvar;
  matvar = Mat_VarRead(mat_file, name);
  if (matvar == NULL)
    return -1;

  strlen = matvar->nbytes;

  if (matvar->dims[0] != 1)
    printf("Error: Multiline string copy attempted\n");

  memcpy(data, matvar->data, strlen);

  Mat_VarFree(matvar);
  return 0;
}

/**********************************************************************/
int matGetDbl (char *name,int n1,int n2, double *data)
{
    matvar_t *matvar;
    matvar = Mat_VarRead(mat_file, name);
    if (matvar == NULL)
      return -1;

    memcpy(data, matvar->data, n1*n2*sizeof(double));
    
    Mat_VarFree(matvar);
    return 0;
}

/**********************************************************************/
int matGetInt (char *name,int n1,int n2, int *data)
{
    matvar_t *matvar;
    matvar = Mat_VarRead(mat_file, name);
    if (matvar == NULL)
      return -1;

    memcpy(data, matvar->data, n1*n2*sizeof(int));
    
    Mat_VarFree(matvar);
    return 0;
}

/**********************************************************************/
int matArrNRow (char *name)
{
  int nrow = 0;
  matvar_t *matvar = Mat_VarRead(mat_file, name);
  if (matvar == NULL)
    return -1;

  nrow = matvar->dims[0];
  Mat_VarFree(matvar);
  return nrow;
}

/**********************************************************************/
int matArrNCol (char *name)
{
  int ncol = 0;
  matvar_t *matvar = Mat_VarRead(mat_file, name);
  if (matvar == NULL)
    return -1;

  ncol = matvar->dims[1];
  Mat_VarFree(matvar);
  return ncol;
}


/**********************************************************************/
/* remove an argument from the list */
void del_arg(int *argc, char* argv[], int j)
{
  int jj;
  for (jj=j+1;jj<*argc;jj++)
    argv[jj-1]=argv[jj];
  (*argc)--;
  argv[*argc]=0;
}
