#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

  /** NOTE: I like the C-style i/o better, so you'll see scanf's,
            fprintf's, etc. throughout. On the other hand, I use
            C++-style memory management features. */ 
             
  /* The labels for these are self-explanatory. */
  FILE *ElePartFile, *EleInFile, *OutFile;
  int  numelems, numproc, i, j;
  char geomfile[120], FileName[120];

  strcpy(geomfile,argv[1]); /*Base geometry file name. You should have the
                              Metis-readable file geomfile, the vertex
                              file geomfile.node, the element file
                              geomfile.ele, and the partition files
                              geomfile.npart.numproc and
                              geomfile.epart.numproc to run this code.*/
  numproc  = atoi(argv[2]); /*Number of processors, into which the initial
                              geometry was split.*/
  i = 0, j = 0;             /*Multipurpose counters.*/
  int *ElementPartition;
  ElementPartition = NULL;  /*Array of element partition labels.*/

  /***Read in number of elements.***/
  sprintf(FileName, "%s.ele", geomfile);
  printf("Reading number of elements from: %s ...", FileName);
  EleInFile = fopen(FileName, "r");
  fscanf(EleInFile, "%d", &numelems);
  fclose(EleInFile);
  printf(" Done.\n");

  /***Read in element partition file.***/
  sprintf(FileName, "%s.epart.%d", geomfile, numproc);
  printf("Reading element partition from: %s ...", FileName);
  ElePartFile = fopen(FileName, "r");
  ElementPartition = new int[numelems];
  for (i=0; i<numelems; i++) {
    fscanf(ElePartFile, "%d", &ElementPartition[i]);
  } 
  fclose(ElePartFile);
  printf(" Done.\n");


  // Write triangle-format partition file.
  sprintf(FileName, "%s.part", geomfile);
  printf("Writing Triangle partitioning info to: %s ...", FileName);
  OutFile = fopen(FileName, "w");
  fprintf(OutFile, "%d %d\n", numelems, numproc);
  for (j=0; j<numelems; j++) {
    fprintf(OutFile, "   %d %d\n", j+1, ElementPartition[j]+1);
  }
  fclose(OutFile);
  printf(" Done.\n");


  delete [] ElementPartition;

  return(0);
}
  
