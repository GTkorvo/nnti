#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

using namespace std;

typedef struct {
  int Gindex;
  double Pos[2];
  int NumOwners;
} NodeType;

typedef struct {
  int Vertex[3];
} ElementType;

typedef struct {
  int Vertex[2];
  int Type;
} EdgeType;


main(int argc, char** argv) {

  /** NOTE: I like the C-style i/o better, so you'll see scanf's,
            fprintf's, etc. throughout. On the other hand, I use
            C++-style memory management features. */ 
             
  /* The labels for these are self-explanatory. */
  FILE *NodePartFile, *ElePartFile, *NodeInFile, *EleInFile,
       *EdgeInFile, *OutFile;
  int  numnodes, numelems, numedges, numproc,
       i, j, k,
       uniquenode_sd, sharednode_sd, interiornode_sd, commonnode_sd,  
       ele_sd, edge_sd, *endnodes_sd, 
       *NodePartition, *ElementPartition,
       *TmpIndx, *TmpIndxEdge, *SortIndx,
       *NumUniqueNode_sd, *NumSharedNode_sd, *NumInteriorNode_sd,
       *NumCommonNode_sd, *NumElem_sd, *NumEdge_sd, *NumBdryNode_sd,
       **AllSortIndx_sd, **InteriorIndx_sd, **CommonIndx_sd, **BdryNode_sd;
  NodeType    *Node, **Node_sd, **DDNode_sd;
  ElementType *Element, *TmpElement, **Element_sd;
  EdgeType *Edge, **Edge_sd;
  char geomfile[120], FileName[120];

  strcpy(geomfile,argv[1]); /*Base geometry file name. You should have the
                              Metis-readable file geomfile, the vertex
                              file geomfile.node, the element file
                              geomfile.ele, and the partition files
                              geomfile.npart.numproc and
                              geomfile.epart.numproc to run this code.*/
  numnodes = atoi(argv[2]); /*Number of nodes (triangle vertices). Must be
                              read in, since not available if quadratic 
                              elements are used.*/ 
  numproc  = atoi(argv[3]); /*Number of processors, into which the initial
                              geometry was split.*/
  i = 0, j = 0, k = 0;      /*Multipurpose counters.*/
  uniquenode_sd = 0;        /*Number of unique nodes per subdomain.*/
  ele_sd = 0;               /*Number of elements per subdomain.*/
  sharednode_sd = 0;        /*Number of shared nodes per subdomain.*/
  edge_sd = 0;              /*Number of edges per subdomain.*/
  numelems = 0;             /*TOTAL number of elements.*/
  numedges = 0;             /*TOTAL number of edges.*/
  NodePartition = NULL;     /*Array of node partition labels.*/
  ElementPartition = NULL;  /*Array of element partition labels.*/
  endnodes_sd = NULL;       /*Points to the end of the unique node list,
                              see below.*/
  TmpIndx = NULL;           /*Array of temporary node indices, see below.*/ 
  TmpIndxEdge = NULL;       /*Array of temporary edge indices, see below.*/ 
  SortIndx = NULL;          /*Array of sorted node indices, see below.*/ 
  NumUniqueNode_sd = NULL;  /*Array of total numbers of unique nodes per
                              subdomain.*/
  NumSharedNode_sd = NULL;  /*Array of total numbers of shared nodes per
                              subdomain.*/
  NumInteriorNode_sd=NULL;  /*Array of total numbers of unique nodes per
                              subdomain.*/
  NumCommonNode_sd = NULL;  /*Array of total numbers of shared nodes per
                              subdomain.*/
  NumBdryNode_sd = NULL;    /*Array of total numbers of boundary nodes per
                              subdomain.*/
  NumElem_sd = NULL;        /*Array of total numbers of elements per
                              subdomain.*/
  NumEdge_sd = NULL;        /*Array of total numbers of edges per
                              subdomain.*/
  Node = NULL;              /*Array containing all nodes.*/
  Element = NULL;           /*Array containing all elements.*/
  TmpElement = NULL;        /*Temporary array of elements, see below.*/
  Edge = NULL;              /*Array containing all edges.*/
  Node_sd = NULL;           /*Array of numproc arrays of nodes in a sd.*/
  Element_sd = NULL;        /*Array of numproc arrays of elems in a sd.*/
  Edge_sd = NULL;           /*Array of numproc arrays of edges in a sd.*/
  BdryNode_sd = NULL;       /*Array of numproc arrays of boundary nodes in a sd.*/

  Node_sd            = new NodeType* [numproc];
  DDNode_sd          = new NodeType* [numproc];
  Element_sd         = new ElementType* [numproc];
  Edge_sd            = new EdgeType* [numproc];
  NumUniqueNode_sd   = new int [numproc];
  NumSharedNode_sd   = new int [numproc];
  NumInteriorNode_sd = new int [numproc];
  NumCommonNode_sd   = new int [numproc];
  NumBdryNode_sd     = new int [numproc];
  NumElem_sd         = new int [numproc];
  NumEdge_sd         = new int [numproc];

  /***Read in input node file.***/
  sprintf(FileName, "%s.node", geomfile);
  printf("Reading nodes from: %s ...", FileName);
  NodeInFile = fopen(FileName, "r");
  //fscanf(NodeInFile, "%*d%*d%*d%*d");
  fscanf(NodeInFile, "%*[^\n]");   /* Skip to the End of the Line */
  fscanf(NodeInFile, "%*1[\n]");   /* Skip One Newline */
  Node = new NodeType[numnodes];
  for (i=0; i<numnodes; i++) {
    //fscanf(NodeInFile, "%d%lf%lf%*d", &Node[i].Gindex, &Node[i].Pos[0],
    //       &Node[i].Pos[1]);
    fscanf(NodeInFile, "%d%lf%lf", &Node[i].Gindex, &Node[i].Pos[0],
           &Node[i].Pos[1]);
    fscanf(NodeInFile, "%*[^\n]");   /* Skip to the End of the Line */
    fscanf(NodeInFile, "%*1[\n]");   /* Skip One Newline */
    Node[i].NumOwners = 0;
  }
  fclose(NodeInFile); 
  printf(" Done.\n");


  /***Read in input element file.***/
  sprintf(FileName, "%s.ele", geomfile);
  printf("Reading elements from: %s ...", FileName);
  EleInFile = fopen(FileName, "r");
  //fscanf(EleInFile, "%d%*d%*d", &numelems);
  fscanf(EleInFile, "%d", &numelems);
  fscanf(EleInFile, "%*[^\n]");   /* Skip to the End of the Line */
  fscanf(EleInFile, "%*1[\n]");   /* Skip One Newline */
  Element = new ElementType[numelems];
  for (i=0; i<numelems; i++) {
    //fscanf(EleInFile, "%*d%d%d%d%*d%*d%*d", &Element[i].Vertex[0],
    //       &Element[i].Vertex[1], &Element[i].Vertex[2]);
    fscanf(EleInFile, "%*d%d%d%d", &Element[i].Vertex[0],
           &Element[i].Vertex[1], &Element[i].Vertex[2]);
    fscanf(EleInFile, "%*[^\n]");   /* Skip to the End of the Line */
    fscanf(EleInFile, "%*1[\n]");   /* Skip One Newline */
  } 
  fclose(EleInFile);
  printf(" Done.\n");


  /***Read in input edge file.***/
  sprintf(FileName, "%s.edge", geomfile);
  printf("Reading edges from: %s ...", FileName);
  EdgeInFile = fopen(FileName, "r");
  //fscanf(EdgeInFile, "%d%*d%", &numedges);
  fscanf(EdgeInFile, "%d", &numedges);
  fscanf(EdgeInFile, "%*[^\n]");   /* Skip to the End of the Line */
  fscanf(EdgeInFile, "%*1[\n]");   /* Skip One Newline */
  Edge = new EdgeType[numedges];
  for (i=0; i<numedges; i++) {
    fscanf(EdgeInFile, "%*d%d%d%d%", &Edge[i].Vertex[0],
           &Edge[i].Vertex[1], &Edge[i].Type);
  } 
  fclose(EdgeInFile);
  printf(" Done.\n");


  /***Read in node partition file.***/
  sprintf(FileName, "%s.npart.%d", geomfile, numproc);
  printf("Reading node partition from: %s ...", FileName);
  NodePartFile = fopen(FileName, "r");
  NodePartition = new int[numnodes];
  for (i=0; i<numnodes; i++) {
    fscanf(NodePartFile, "%d", &NodePartition[i]);
  } 
  fclose(NodePartFile);
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


  /*** Chop up the global info into local subdomain files (which
       still contain global indices) ***/

  // This is an upper bound for the total number of nodes, elements,
  // and edges per single subdomain.
  TmpIndx = new int[numnodes];
  SortIndx = new int[numnodes];
  TmpElement = new ElementType[numelems];
  TmpIndxEdge = new int[numedges];

  for (i=0; i<numproc; i++) {

    /*First, extract global indices for subdomain array of nodes.*/
    uniquenode_sd = 0; 
    for (j=0; j<numnodes; j++) {
      if (NodePartition[j] == i) {
        TmpIndx[uniquenode_sd]        = Node[j].Gindex;
        uniquenode_sd++;
      }
    }

    /*Second, extract the subdomain array of elements.*/
    ele_sd = 0; 
    for (j=0; j<numelems; j++) {
      if (ElementPartition[j] == i) {
        TmpElement[ele_sd].Vertex[0] = Element[j].Vertex[0];
        TmpElement[ele_sd].Vertex[1] = Element[j].Vertex[1];
        TmpElement[ele_sd].Vertex[2] = Element[j].Vertex[2];
        ele_sd++;
      }
    }

    /*Determine shared nodes in the current subdomain.*/
    sharednode_sd = 0;
    for (j=0; j<ele_sd; j++) {
      for (k=0; k<3; k++) {
        if (! binary_search(&(TmpIndx[0]), &(TmpIndx[uniquenode_sd]),
            TmpElement[j].Vertex[k])) {
          TmpIndx[uniquenode_sd+sharednode_sd] = TmpElement[j].Vertex[k];
          sharednode_sd++;
        }
      }
    }
    sort(&(TmpIndx[uniquenode_sd]),
         &(TmpIndx[uniquenode_sd+sharednode_sd]));
    endnodes_sd = unique(&(TmpIndx[uniquenode_sd]),
                         &(TmpIndx[uniquenode_sd+sharednode_sd]));
    /* Do some pointer arithmetic to determine the number of additional
       (shared) nodes. Note: the previous figure included duplicates.*/
    sharednode_sd = endnodes_sd-&(TmpIndx[uniquenode_sd]);

    /*Determine edges in the current subdomain.*/
    // The following only works if there are no elements that belong to
    // a domain, but are completely (on all sides) surrounded by other
    // domains.
    //for (j=0; j<uniquenode_sd+sharednode_sd; j++)
    //  SortIndx[j] = TmpIndx[j];
    //sort(&(SortIndx[0]), &(SortIndx[uniquenode_sd+sharednode_sd]));
    //edge_sd = 0;
    //for (j=0; j<numedges; j++) {
    //  if (binary_search(&(SortIndx[0]),
    //        &(SortIndx[uniquenode_sd+sharednode_sd]), Edge[j].Vertex[0])
    //      &&
    //      binary_search(&(SortIndx[0]),
    //        &(SortIndx[uniquenode_sd+sharednode_sd]), Edge[j].Vertex[1]))
    //  {
    //      TmpIndxEdge[edge_sd] = j;
    //      edge_sd++;
    //  }
    //}
    edge_sd = 0;
    for (j=0; j<numedges; j++) {
      // only need boundary edges
      if (Edge[j].Type > 0) {
        for (k=0; k<ele_sd; k++) {
          if ( ((Edge[j].Vertex[0] == TmpElement[k].Vertex[0]) ||
                (Edge[j].Vertex[0] == TmpElement[k].Vertex[1]) ||
                (Edge[j].Vertex[0] == TmpElement[k].Vertex[2])) &&
               ((Edge[j].Vertex[1] == TmpElement[k].Vertex[0]) ||
                (Edge[j].Vertex[1] == TmpElement[k].Vertex[1]) ||
                (Edge[j].Vertex[1] == TmpElement[k].Vertex[2])) ) {
 
               TmpIndxEdge[edge_sd] = j;
               edge_sd++;
          }
        }
      }
    }

    /*Record node info per subdomain.*/
    Node_sd[i] = new NodeType[uniquenode_sd+sharednode_sd];
    NumUniqueNode_sd[i] = uniquenode_sd;
    NumSharedNode_sd[i] = sharednode_sd;
    j = 0;
    while (&TmpIndx[j] != endnodes_sd) {
      Node_sd[i][j].Gindex = TmpIndx[j];
      Node_sd[i][j].Pos[0] = Node[TmpIndx[j]-1].Pos[0];
      Node_sd[i][j].Pos[1] = Node[TmpIndx[j]-1].Pos[1];
      Node_sd[i][j].NumOwners = 0;
      j++;
    }
    /*Record element info per subdomain.*/
    Element_sd[i] = new ElementType[ele_sd];
    NumElem_sd[i] = ele_sd;
    for (j=0; j<ele_sd; j++) {
      Element_sd[i][j].Vertex[0] = TmpElement[j].Vertex[0];
      Element_sd[i][j].Vertex[1] = TmpElement[j].Vertex[1];
      Element_sd[i][j].Vertex[2] = TmpElement[j].Vertex[2];
    }
    /*Record edge info per subdomain.*/
    Edge_sd[i] = new EdgeType[edge_sd];
    NumEdge_sd[i] = edge_sd;
    for (j=0; j<edge_sd; j++) {
      Edge_sd[i][j].Vertex[0] = Edge[TmpIndxEdge[j]].Vertex[0];
      Edge_sd[i][j].Vertex[1] = Edge[TmpIndxEdge[j]].Vertex[1];
      Edge_sd[i][j].Type      = Edge[TmpIndxEdge[j]].Type;
    }

  }


  AllSortIndx_sd = new int* [numproc];
  InteriorIndx_sd = new int* [numproc];
  CommonIndx_sd = new int* [numproc];
  BdryNode_sd = new int* [numproc];
  
  for (i=0; i<numproc; i++) {
    int num = NumUniqueNode_sd[i]+NumSharedNode_sd[i];     
    AllSortIndx_sd[i] = new int[num];
    for (j=0; j<num; j++) {
      AllSortIndx_sd[i][j] = Node_sd[i][j].Gindex;
    }
    sort(AllSortIndx_sd[i], AllSortIndx_sd[i]+num);
  }

  for (i=0; i<numproc; i++) {
    int num = NumUniqueNode_sd[i]+NumSharedNode_sd[i];     


    /*Determine common nodes, i.e. the ones on the boundary between
      the current and any other subdomain.*/

    // WARNING: An unverified upper bound. Might need adjusting !!!
    CommonIndx_sd[i] = new int[4*num];

    int totalcommon = 0;
    int * jcommon = CommonIndx_sd[i];
    for (j=0; j<i; j++) {
      int jnum = NumUniqueNode_sd[j]+NumSharedNode_sd[j];
      jcommon = set_intersection(
                  AllSortIndx_sd[i], AllSortIndx_sd[i]+num,
                  AllSortIndx_sd[j], AllSortIndx_sd[j]+jnum,
                  jcommon);
    }
    for (j=i+1; j<numproc; j++) {
      int jnum = NumUniqueNode_sd[j]+NumSharedNode_sd[j];
      jcommon = set_intersection(
                  AllSortIndx_sd[i], AllSortIndx_sd[i]+num,
                  AllSortIndx_sd[j], AllSortIndx_sd[j]+jnum,
                  jcommon);
    }
    totalcommon = jcommon - CommonIndx_sd[i]; 
    sort(CommonIndx_sd[i], CommonIndx_sd[i]+totalcommon);
    jcommon = unique(CommonIndx_sd[i], CommonIndx_sd[i]+totalcommon);
    totalcommon = jcommon - CommonIndx_sd[i];
    NumCommonNode_sd[i] = totalcommon;


    /*Determine nodes interior to a subdomain (including unshared
      boundaries).*/

    InteriorIndx_sd[i] = new int[num];

    int totalinterior = 0;
    set_difference(AllSortIndx_sd[i], AllSortIndx_sd[i] + num,
                   CommonIndx_sd[i], CommonIndx_sd[i] + totalcommon,
                   InteriorIndx_sd[i]);
    totalinterior = num - totalcommon;
    NumInteriorNode_sd[i] = totalinterior;


    /* Determine all boundary vertices in a subdomain/ */
    int * lastindx = 0;
    NumBdryNode_sd[i] = NumEdge_sd[i];
    BdryNode_sd[i] = new int[2*NumEdge_sd[i]];
    for (j=0; j<NumEdge_sd[i]; j++) {
      BdryNode_sd[i][j] = Edge_sd[i][j].Vertex[0];
      BdryNode_sd[i][j+NumEdge_sd[i]] = Edge_sd[i][j].Vertex[1];
    }
    sort(BdryNode_sd[i], BdryNode_sd[i]+2*NumEdge_sd[i]);
    lastindx  = unique(BdryNode_sd[i], BdryNode_sd[i]+2*NumEdge_sd[i]);
    NumBdryNode_sd[i] = lastindx - BdryNode_sd[i];


    /* Record the interior/common node info. */
    DDNode_sd[i] = new NodeType[num];
    for (j=0; j<totalinterior; j++) {
      DDNode_sd[i][j].Gindex = InteriorIndx_sd[i][j];
      DDNode_sd[i][j].Pos[0] = Node[InteriorIndx_sd[i][j]-1].Pos[0];
      DDNode_sd[i][j].Pos[1] = Node[InteriorIndx_sd[i][j]-1].Pos[1];   
      DDNode_sd[i][j].NumOwners = 0;
    }
    for (j=totalinterior; j<num; j++) {
      DDNode_sd[i][j].Gindex = CommonIndx_sd[i][j-totalinterior];
      DDNode_sd[i][j].Pos[0] =
        Node[CommonIndx_sd[i][j-totalinterior]-1].Pos[0];
      DDNode_sd[i][j].Pos[1] =
        Node[CommonIndx_sd[i][j-totalinterior]-1].Pos[1];   
      DDNode_sd[i][j].NumOwners = 0;
    }
    
  }



  /*Update number of owners for each node.*/
  for (i=0; i<numproc; i++) {
    for(j=0; j<NumUniqueNode_sd[i]+NumSharedNode_sd[i]; j++) {
      Node[Node_sd[i][j].Gindex-1].NumOwners++;
    }
  }

  for (i=0; i<numproc; i++) {

    // Write subdomain file.
    sprintf(FileName, "%s.%03d", geomfile, i);
    printf("Writing subdomain partition to: %s ...", FileName);
    OutFile = fopen(FileName, "w");

    /*Write node info for the nonoverlapping case.*/
    fprintf(OutFile, "%d %d\n", NumUniqueNode_sd[i], NumSharedNode_sd[i]);
    for (j=0; j<NumUniqueNode_sd[i]+NumSharedNode_sd[i]; j++) {
      fprintf(OutFile, "%d %16.15e %16.15e %d\n", Node_sd[i][j].Gindex,
              Node_sd[i][j].Pos[0], Node_sd[i][j].Pos[1],
              Node[Node_sd[i][j].Gindex-1].NumOwners);
    }

    /*Write node info for the overlapping case.*/
    fprintf(OutFile, "%d %d\n", NumInteriorNode_sd[i], NumCommonNode_sd[i]);
    for (j=0; j<NumInteriorNode_sd[i]+NumCommonNode_sd[i]; j++) {
      fprintf(OutFile, "%d %16.15e %16.15e %d\n", DDNode_sd[i][j].Gindex,
              DDNode_sd[i][j].Pos[0], DDNode_sd[i][j].Pos[1],
              Node[DDNode_sd[i][j].Gindex-1].NumOwners);
    }       

    /*Write element info.*/
    fprintf(OutFile, "%d\n", NumElem_sd[i]);
    for (j=0; j<NumElem_sd[i]; j++) {
      fprintf(OutFile, "%d %d %d\n", Element_sd[i][j].Vertex[0],
              Element_sd[i][j].Vertex[1], Element_sd[i][j].Vertex[2]); 
    }

    /*Write edge info.*/
    fprintf(OutFile, "%d\n", NumEdge_sd[i]);
    for (j=0; j<NumEdge_sd[i]; j++) {
      fprintf(OutFile, "%d %d %d\n", Edge_sd[i][j].Vertex[0],
              Edge_sd[i][j].Vertex[1], Edge_sd[i][j].Type); 
    }
    printf(" Done.\n");

    /*Write boundary node info.*/
    // Include this if you need it in the preprocessing step.
    // fprintf(OutFile, "%d\n", NumBdryNode_sd[i]);
    // for (j=0; j<NumBdryNode_sd[i]; j++) {
    //   fprintf(OutFile, "%d\n", BdryNode_sd[i][j]); 
    // }

    printf(" Done.\n");
    fclose(OutFile);

  }


  delete [] Node;
  delete [] Element;
  delete [] Edge;
  delete [] NodePartition;
  delete [] ElementPartition;
  delete [] TmpIndx;
  delete [] SortIndx;
  delete [] TmpElement;
  delete [] NumUniqueNode_sd;
  delete [] NumSharedNode_sd;
  delete [] NumInteriorNode_sd;
  delete [] NumCommonNode_sd;
  delete [] NumBdryNode_sd;
  delete [] NumElem_sd;
  delete [] NumEdge_sd;
  delete [] Node_sd;
  delete [] AllSortIndx_sd;
  delete [] BdryNode_sd;
  delete [] Element_sd;
  delete [] Edge_sd;



  return(0);
}
  
