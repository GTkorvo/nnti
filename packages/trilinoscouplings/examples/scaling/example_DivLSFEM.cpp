// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_01.cpp
    \brief  Example creation of mass and stiffness matrices for div-curl system on a hexadedral mesh using div-conforming elements.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
    \remark Sample command line
    \code   ./example_DivLSFEM.exe 10 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII.h"
#include "../mesh_spec_lt/im_ne_nemesisI.h"

using namespace std;
using namespace Intrepid;

// Functions to evaluate exact solution and its derivatives
int evalu(double & uExact0, 
          double & uExact1, 
          double & uExact2, 
          double & x, 
          double & y, 
          double & z);

double evalDivu(double & x, 
                double & y, 
                double & z);

int evalCurlu(double & curlu0, 
              double & curlu1, 
              double & curlu2, 
              double & x, 
              double & y, 
              double & z, 
              double & mu);

int evalCurlCurlu(double & curlCurlu0, 
                  double & curlCurlu1, 
                  double & curlCurlu2, 
                  double & x, 
                  double & y, 
                  double & z, 
                  double & mu);

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 1 ),
                      std::invalid_argument,
                      ">>> ERROR (example_01): Invalid number of arguments. See code listing for requirements.");
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 12)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Div-Curl System on Hexahedral Mesh                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


// ************************************ GET INPUTS **************************************

  // Input file
    std::string   xmlInFileName = "input.xml";

  // Read xml file into parameter list
    Teuchos::ParameterList inputList;

   if(xmlInFileName.length()) {
      std::cout << "\nReading a parameter list from the XML file \""<<xmlInFileName<<"\" ...\n";
      Teuchos::updateParametersFromXmlFile(xmlInFileName,&inputList);
      std::cout << "\nParameter list read from the XML file \""<<xmlInFileName<<"\":\n\n";
      inputList.print(std::cout,2,true,true);
    }

  // Get pamgen mesh definition
    std::string meshInput = Teuchos::getParameter<std::string>(inputList,"meshInput");
    std::cout << meshInput << "\n";
 

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int numEdgesPerElem = hex_8.getEdgeCount();
    int numFacesPerElem = hex_8.getSideCount();
    int numNodesPerEdge = 2;
    int numNodesPerFace = 4;
    int numEdgesPerFace = 4;
    int spaceDim = hex_8.getDimension();

   // Build reference element edge to node map
    FieldContainer<int> refEdgeToNode(numEdgesPerElem,numNodesPerEdge);
    for (int i=0; i<numEdgesPerElem; i++){
        refEdgeToNode(i,0)=hex_8.getNodeMap(1, i, 0);
        refEdgeToNode(i,1)=hex_8.getNodeMap(1, i, 1);
    }

   // Build reference element face to node map
    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
        refFaceToNode(i,0)=hex_8.getNodeMap(2, i, 0);
        refFaceToNode(i,1)=hex_8.getNodeMap(2, i, 1);
        refFaceToNode(i,2)=hex_8.getNodeMap(2, i, 2);
        refFaceToNode(i,3)=hex_8.getNodeMap(2, i, 3);
    }


// *********************************** GENERATE MESH ************************************

    std::cout << "Generating mesh ... \n\n";

    int NX = 20;
    int NY = 20;
    int NZ = 20;

  
   // Generate mesh with Pamgen
    int dim=3;
    int rank=0;
    int numProcs=1;
    long int maxInt = 100000000;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
    
   // Get mesh size info
    char title[100];
    int numDim;
    int numNodes;
    int numElems;
    int numElemBlk;
    int numNodeSets;
    int numSideSets;
    int id = 0;

    im_ex_get_init(id, title, &numDim, &numNodes, 
                                &numElems, &numElemBlk, &numNodeSets,
                                &numSideSets);

   // Read node coordinates and place in field container
    FieldContainer<double> nodeCoord(numNodes,dim);
    double * nodeCoordx = new double [numNodes];
    double * nodeCoordy = new double [numNodes];
    double * nodeCoordz = new double [numNodes];
    im_ex_get_coord(id,nodeCoordx,nodeCoordy,nodeCoordz);
    for (int i=0; i<numNodes; i++) {          
      nodeCoord(i,0)=nodeCoordx[i];
      nodeCoord(i,1)=nodeCoordy[i];
      nodeCoord(i,2)=nodeCoordz[i];
    }
    delete [] nodeCoordx;
    delete [] nodeCoordy;
    delete [] nodeCoordz;

   // Get element block information
    int * blockIds = new int [numElemBlk];
    im_ex_get_elem_blk_ids(id, blockIds);

    int  *nodesPerElement    = new int [numElemBlk];
    int  *elementAttributes  = new int [numElemBlk];
    int  *elements           = new int [numElemBlk];
    char **elementTypes      = new char * [numElemBlk];
    int  **elmtNodeLinkage   = new int * [numElemBlk];

    for(int b = 0; b < numElemBlk; b ++){
      elementTypes[b] = new char [MAX_STR_LENGTH + 1];
      im_ex_get_elem_block(id,
                              blockIds[b],
                              elementTypes[b],
                              (int*)&(elements[b]),
                              (int*)&(nodesPerElement[b]),
                              (int*)&(elementAttributes[b]));
    }

    for(int b = 0; b < numElemBlk; b++){
      elmtNodeLinkage[b] =  new int [nodesPerElement[b]*elements[b]];
      im_ex_get_elem_conn(id,blockIds[b],elmtNodeLinkage[b]);
    }

   // Get mu value for each block of elements from parameter list
    double  *mu = new double [numElemBlk];
    for(int b = 0; b < numElemBlk; b++){
       stringstream muBlock;
       muBlock.clear();
       muBlock << "mu" << b;
       std::cout << muBlock.str() << "\n";
       mu[b] = inputList.get(muBlock.str(),1.0);
    }

   // Get node-element connectivity and set mu value for element
    int ielem = 0;
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    FieldContainer<double> muVal(numElems);
    for(int b = 0; b < numElemBlk; b++){
      for(int el = 0; el < elements[b]; el++){
        for (int j=0; j<numNodesPerElem; j++) {
          elemToNode(ielem,j) = elmtNodeLinkage[b][el*numNodesPerElem + j]-1;
        }
        muVal(ielem) = mu[b];     
        ielem ++;
      }
    }

   // Print mesh information  
    int numEdges = (NX)*(NY + 1)*(NZ + 1) + (NX + 1)*(NY)*(NZ + 1) + (NX + 1)*(NY + 1)*(NZ);
    int numFaces = (NX)*(NY)*(NZ + 1) + (NX)*(NY + 1)*(NZ) + (NX + 1)*(NY)*(NZ);
    std::cout << " Number of Elements: " << numElems << " \n";
    std::cout << "    Number of Nodes: " << numNodes << " \n";
    std::cout << "    Number of Edges: " << numEdges << " \n";
    std::cout << "    Number of Faces: " << numFaces << " \n\n";
  
  // Get edge connectivity
    FieldContainer<int> edgeToNode(numEdges, numNodesPerEdge);
    FieldContainer<int> elemToEdge(numElems, numEdgesPerElem);
    int iedge = 0;
    int inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
           if (i < NX){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + 1;
               if (j < NY && k < NZ){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,0) = iedge;
                  if (j > 0)
                     elemToEdge(ielem-NX,2) = iedge; 
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,4) = iedge; 
                  if (j > 0 && k > 0)
                     elemToEdge(ielem-NX*NY-NX,6) = iedge; 
                }
               else if (j == NY && k == NZ){
                  ielem=i+(NY-1)*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,6) = iedge;
                }
               else if (k == NZ && j < NY){
                  ielem=i+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,4) = iedge;
                  if (j > 0)
                    elemToEdge(ielem-NX,6) = iedge;
                }
               else if (k < NZ && j == NY){
                  ielem=i+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,2) = iedge;
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,6) = iedge;
                }
               iedge++;
            }
           if (j < NY){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + NX+1;
               if (i < NX && k < NZ){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,3) = iedge;
                  if (i > 0)
                     elemToEdge(ielem-1,1) = iedge; 
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,7) = iedge; 
                  if (i > 0 && k > 0)
                     elemToEdge(ielem-NX*NY-1,5) = iedge; 
                }
               else if (i == NX && k == NZ){
                  ielem=NX-1+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,5) = iedge;
                }
               else if (k == NZ && i < NX){
                  ielem=i+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,7) = iedge;
                  if (i > 0)
                    elemToEdge(ielem-1,5) = iedge;
                }
               else if (k < NZ && i == NX){
                  ielem=NX-1+j*NX+k*NX*NY;
                  elemToEdge(ielem,1) = iedge;
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,5) = iedge;
                }
               iedge++;
            }
           if (k < NZ){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + (NX+1)*(NY+1);
               if (i < NX && j < NY){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,8) = iedge;
                  if (i > 0)
                     elemToEdge(ielem-1,9) = iedge; 
                  if (j > 0)
                     elemToEdge(ielem-NX,11) = iedge; 
                  if (i > 0 && j > 0)
                     elemToEdge(ielem-NX-1,10) = iedge; 
                }
               else if (i == NX && j == NY){
                  ielem=NX-1+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,10) = iedge;
                }
               else if (j == NY && i < NX){
                  ielem=i+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,11) = iedge;
                  if (i > 0)
                    elemToEdge(ielem-1,10) = iedge;
                }
               else if (j < NY && i == NX){
                  ielem=NX-1+j*NX+k*NX*NY;
                  elemToEdge(ielem,9) = iedge;
                  if (j > 0)
                     elemToEdge(ielem-NX,10) = iedge;
                }
               iedge++;
            }
            inode++;
         }
      }
   }

   // Get face connectivity
    FieldContainer<int> faceToNode(numFaces, numNodesPerFace);
    FieldContainer<int> elemToFace(numElems, numFacesPerElem);
    FieldContainer<int> faceToEdge(numFaces, numEdgesPerFace);
    int iface = 0;
    inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
           if (i < NX && k < NZ) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + 1;
              faceToNode(iface,2)=inode + (NX+1)*(NY+1)+1;
              faceToNode(iface,3)=inode + (NX+1)*(NY+1);
              if (j < NY) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,0);
                 faceToEdge(iface,1)=elemToEdge(ielem,9);
                 faceToEdge(iface,2)=elemToEdge(ielem,4);
                 faceToEdge(iface,3)=elemToEdge(ielem,8);
                 elemToFace(ielem,0)=iface;
                 if (j > 0) {
                    elemToFace(ielem-NX,2)=iface;
                 }
              }
              else if (j == NY) {
                 ielem=i+(NY-1)*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,2);
                 faceToEdge(iface,1)=elemToEdge(ielem,10);
                 faceToEdge(iface,2)=elemToEdge(ielem,6);
                 faceToEdge(iface,3)=elemToEdge(ielem,11);
                 elemToFace(ielem,2)=iface;
              }
              iface++;
           }
           if (j < NY && k < NZ) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + NX+1;
              faceToNode(iface,2)=inode + (NX+1)*(NY+1) + NX+1;
              faceToNode(iface,3)=inode + (NX+1)*(NY+1);
              if (i < NX) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,3);
                 faceToEdge(iface,1)=elemToEdge(ielem,11);
                 faceToEdge(iface,2)=elemToEdge(ielem,7);
                 faceToEdge(iface,3)=elemToEdge(ielem,8);
                 elemToFace(ielem,3)=iface;
                 if (i > 0) {
                    elemToFace(ielem-1,1)=iface;
                 }
              }
              else if (i == NX) {
                 ielem=NX-1+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,1);
                 faceToEdge(iface,1)=elemToEdge(ielem,10);
                 faceToEdge(iface,2)=elemToEdge(ielem,5);
                 faceToEdge(iface,3)=elemToEdge(ielem,9);
                 elemToFace(ielem,1)=iface;
              }
              iface++;
           }
           if (i < NX && j < NY) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + 1;
              faceToNode(iface,2)=inode + NX+2;
              faceToNode(iface,3)=inode + NX+1;
              if (k < NZ) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,0);
                 faceToEdge(iface,1)=elemToEdge(ielem,1);
                 faceToEdge(iface,2)=elemToEdge(ielem,2);
                 faceToEdge(iface,3)=elemToEdge(ielem,3);
                 elemToFace(ielem,4)=iface;
                 if (k > 0) {
                    elemToFace(ielem-NX*NY,5)=iface;
                 }
              }
              else if (k == NZ) {
                 ielem=i+j*NX+(NZ-1)*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,4);
                 faceToEdge(iface,1)=elemToEdge(ielem,5);
                 faceToEdge(iface,2)=elemToEdge(ielem,6);
                 faceToEdge(iface,3)=elemToEdge(ielem,7);
                 elemToFace(ielem,5)=iface;
              }
              iface++;
           }
          inode++;
         }
      }
   }
 
   // Output element to face connectivity
    ofstream el2fout("elem2face.dat");
    ofstream el2nout("elem2node.dat");
    for (int i=0; i<numElems; i++) {
      for (int l=0; l<numFacesPerElem; l++) {
         el2fout << elemToFace(i,l) << "  ";
      } 
      el2fout << "\n";
      for (int m=0; m<numNodesPerElem; m++) {
        el2nout << elemToNode(i,m) << "  ";
      } 
      el2nout << "\n";
    }
    el2fout.close();
    el2nout.close();

   // Output face to edge and face to node connectivity
    ofstream f2edout("face2edge.dat");
    ofstream f2nout("face2node.dat");
    for (int k=0; k<numFaces; k++) {
       for (int i=0; i<numEdgesPerFace; i++) {
           f2edout << faceToEdge(k,i) << "  ";
       } 
       for (int j=0; j<numNodesPerFace; j++) {
           f2nout << faceToNode(k,j) << "  ";
       } 
       f2edout << "\n";
       f2nout << "\n";
    }
    f2edout.close();
    f2nout.close();

   // Container indicating whether a face is on the boundary (1-yes 0-no)
    FieldContainer<int> faceOnBoundary(numFaces);
    FieldContainer<int> edgeOnBoundary(numEdges);

 // Get boundary (side set) information
    int * sideSetIds = new int [numSideSets];
    int numSidesInSet;
    int numDFinSet;
    im_ex_get_side_set_ids(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
        if (numSidesInSet > 0){
          int * sideSetElemList = new int [numSidesInSet];
          int * sideSetSideList = new int [numSidesInSet];
          im_ex_get_side_set(id,sideSetIds[i],sideSetElemList,sideSetSideList);
          for (int j=0; j<numSidesInSet; j++) {
             int iface = sideSetSideList[j]-1;
             faceOnBoundary(elemToFace(sideSetElemList[j]-1,iface))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),0))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),1))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),2))=1;
             edgeOnBoundary(faceToEdge(elemToFace(sideSetElemList[j]-1,iface),3))=1;
          }
          delete [] sideSetElemList;
          delete [] sideSetSideList;
       }
    }

    delete [] sideSetIds;

   //TEMP
    ofstream fFaceout("faceOnBndy.dat");
    for (int i=0; i<numFaces; i++){
       fFaceout << faceOnBoundary(i) <<"\n";
    }
    fFaceout.close();
    ofstream fEdgeout("edgeOnBndy.dat");
    for (int i=0; i<numEdges; i++){
       fEdgeout << edgeOnBoundary(i) <<"\n";
    }
    fEdgeout.close();

    // Print coords
    FILE *f=fopen("coords.dat","w");
    inode=0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {          
          fprintf(f,"%22.16e %22.16e %22.16e\n",nodeCoord(inode,0),nodeCoord(inode,1),nodeCoord(inode,2));
          inode++;
        }
      }
    }
    fclose(f);


// **************************** INCIDENCE MATRIX **************************************

   // Edge to face incidence matrix
    std::cout << "Building incidence matrix ... \n\n";

    Epetra_SerialComm Comm;
    Epetra_Map globalMapD(numFaces, 0, Comm);
    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_Map globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix DCurl(Copy, globalMapD, globalMapC, 2);

    double vals[4];
    vals[0]=0.5; vals[1]=0.5; vals[2]=-0.5; vals[3]=-0.5;
    for (int j=0; j<numFaces; j++){
        int rowNum = j;
        int colNum[4];
        colNum[0] = faceToEdge(j,0);
        colNum[1] = faceToEdge(j,1);
        colNum[2] = faceToEdge(j,2);
        colNum[3] = faceToEdge(j,3);
        DCurl.InsertGlobalValues(1, &rowNum, 4, colNum, vals);
    }


// ************************************ CUBATURE **************************************

   // Get numerical integration points and weights
    std::cout << "Getting cubature ... \n\n";

    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);

   // Get numerical integration points and weights for hexahedron face
    //             (needed for rhs boundary term)

    // Define topology of the face parametrization domain as [-1,1]x[-1,1]
    CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Define cubature
    DefaultCubatureFactory<double>  cubFactoryFace;
    Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactoryFace.create(paramQuadFace, 3);
    int cubFaceDim    = hexFaceCubature -> getDimension();
    int numFacePoints = hexFaceCubature -> getNumPoints();

    // Define storage for cubature points and weights on [-1,1]x[-1,1]
    FieldContainer<double> paramGaussWeights(numFacePoints);
    FieldContainer<double> paramGaussPoints(numFacePoints,cubFaceDim);

    // Define storage for cubature points on workset faces
    hexFaceCubature -> getCubature(paramGaussPoints, paramGaussWeights);


// ************************************** BASIS ***************************************

   // Define basis 
    std::cout << "Getting basis ... \n\n";
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexHDivBasis;
    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsD = hexHDivBasis.getCardinality();
    int numFieldsG = hexHGradBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexDVals(numFieldsD, numCubPoints, spaceDim); 
     FieldContainer<double> hexDivs(numFieldsD, numCubPoints); 
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> worksetDVals(numFieldsD, numFacePoints, spaceDim); 

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDivs, cubPoints, OPERATOR_DIV);
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


    std::cout << "Building mass and stiffness matrices ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Containers for nodes, edge and face signs 
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
    FieldContainer<double> hexEdgeSigns(numCells, numFieldsC);
    FieldContainer<double> hexFaceSigns(numCells, numFieldsD);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> weightedMeasureMu(numCells, numCubPoints);
    FieldContainer<double> hexCValsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCValsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for element HDIV mass matrix
    FieldContainer<double> massMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> hexDValsTransformed(numCells, numFieldsD, numCubPoints, spaceDim);
    FieldContainer<double> hexDValsTransformedWeighted(numCells, numFieldsD, numCubPoints, spaceDim);
   // Containers for element HDIV stiffness matrix
    FieldContainer<double> stiffMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> hexDivsTransformed(numCells, numFieldsD, numCubPoints);
    FieldContainer<double> hexDivsTransformedWeighted(numCells, numFieldsD, numCubPoints);
   // Containers for element HGRAD mass matrix
    FieldContainer<double> massMatrixG(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah(numCells, numCubPoints);
    FieldContainer<double> gD(numCells, numFieldsD);
    FieldContainer<double> hD(numCells, numFieldsD);
    FieldContainer<double> gDBoundary(numCells, numFieldsD);
    FieldContainer<double> refGaussPoints(numFacePoints,spaceDim);
    FieldContainer<double> worksetGaussPoints(numCells,numFacePoints,spaceDim);
    FieldContainer<double> worksetJacobians(numCells, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet(numCells, numFacePoints);
    FieldContainer<double> worksetFaceTu(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceTv(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceN(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetVFieldVals(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDValsTransformed(numCells, numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> curluFace(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDataCrossField(numCells, numFieldsD, numFacePoints, spaceDim);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix MassD(Copy, globalMapD, numFieldsD);
    Epetra_FECrsMatrix MassG(Copy, globalMapG, numFieldsG);
    Epetra_FECrsMatrix StiffD(Copy, globalMapD, numFieldsD);
    Epetra_FEVector rhsD(globalMapD);

    ofstream fSignsout("faceSigns.dat");

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

     // Face signs
      for (int j=0; j<numFacesPerElem; j++) {
         hexFaceSigns(0,j) = -1.0;
         for (int i=0; i<numNodesPerFace; i++) {
           int indf=i+1;
           if (indf > numNodesPerFace) indf=0;
           if (elemToNode(k,refFaceToNode(j,0))==faceToNode(elemToFace(k,j),i) &&
               elemToNode(k,refFaceToNode(j,1))==faceToNode(elemToFace(k,j),indf))
                hexFaceSigns(0,j) = 1.0;
          }
         fSignsout << hexFaceSigns(0,j) << "  ";
       }
       fSignsout << "\n";

     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;
       }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );


// ************************** Compute element HCurl mass matrices *******************************

     // transform to physical coordinates 
      fst::HCURLtransformVALUE<double>(hexCValsTransformed, hexJacobInv, 
                                   hexCVals);

     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMu(nC,nPt) = weightedMeasure(nC,nPt) * muVal(k);
        }
      }

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCValsTransformedWeighted,
                                   weightedMeasureMu, hexCValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixC,
                             hexCValsTransformed, hexCValsTransformedWeighted,
                             COMP_BLAS);
     // apply edge signs
      fst::applyLeftFieldSigns<double>(massMatrixC, hexEdgeSigns);
      fst::applyRightFieldSigns<double>(massMatrixC, hexEdgeSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************** Compute element HDiv mass matrices *******************************

     // transform to physical coordinates 
      fst::HDIVtransformVALUE<double>(hexDValsTransformed, hexJacobian, hexJacobDet,
                                   hexDVals);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDValsTransformedWeighted,
                                   weightedMeasure, hexDValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixD,
                             hexDValsTransformed, hexDValsTransformedWeighted,
                             COMP_BLAS);
     // apply face signs
      fst::applyLeftFieldSigns<double>(massMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(massMatrixD, hexFaceSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = massMatrixD(0,row,col);
            MassD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************ Compute element HDiv stiffness matrices *****************************

      // transform to physical coordinates 
      fst::HDIVtransformDIV<double>(hexDivsTransformed, hexJacobDet,
                                    hexDivs);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDivsTransformedWeighted,
                                   weightedMeasure, hexDivsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(stiffMatrixD,
                             hexDivsTransformed, hexDivsTransformedWeighted,
                             COMP_BLAS);

     // apply face signs
      fst::applyLeftFieldSigns<double>(stiffMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(stiffMatrixD, hexFaceSigns);

     // assemble into global matrix
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = stiffMatrixD(0,row,col);
            StiffD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }
// ************************** Compute element HGrad mass matrices *******************************

     // transform to physical coordinates
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixG,
                             hexGValsTransformed, hexGValsTransformedWeighted, COMP_BLAS);

      // assemble into global matrix
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = massMatrixG(0,row,col);
            MassG.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }


// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
       FieldContainer<double> rhsDatah(numCells, numCubPoints);
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          double du1, du2, du3;

          evalCurlCurlu(du1, du2, du3, x, y, z, muVal(k));
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          rhsDatah(0,nPt) = evalDivu(x, y, z);
       }

     // integrate (g,curl w) term
      fst::integrate<double>(gD, rhsDatag, hexDValsTransformedWeighted,
                             COMP_BLAS);

     // integrate (h,div w) term
      fst::integrate<double>(hD, rhsDatah, hexDivsTransformedWeighted,
                             COMP_BLAS);

     // apply signs
      fst::applyFieldSigns<double>(gD, hexFaceSigns);
      fst::applyFieldSigns<double>(hD, hexFaceSigns);

     // calculate boundary term
      for (int i = 0; i < numFacesPerElem; i++){
        if (faceOnBoundary(elemToFace(k,i))){

         // map Gauss points on quad to reference face: paramGaussPoints -> refGaussPoints
            CellTools::mapToReferenceSubcell(refGaussPoints,
                                   paramGaussPoints,
                                   2, i, hex_8);

         // get basis values at points on reference cell
            hexHDivBasis.getValues(worksetDVals, refGaussPoints, OPERATOR_VALUE);

         // compute Jacobians at Gauss pts. on reference face for all parent cells
            CellTools::setJacobian(worksetJacobians, refGaussPoints,
                         hexNodes, hex_8);
            CellTools::setJacobianDet(worksetJacobDet, worksetJacobians);

         // transform to physical coordinates
            fst::HDIVtransformVALUE<double>(worksetDValsTransformed, worksetJacobians,
                                   worksetJacobDet, worksetDVals);

         // map Gauss points on quad from ref. face to face workset: refGaussPoints -> worksetGaussPoints
            CellTools::mapToPhysicalFrame(worksetGaussPoints,
                                refGaussPoints,
                                hexNodes, hex_8);

         // compute face tangents
            CellTools::getPhysicalFaceTangents(worksetFaceTu,
                                     worksetFaceTv,
                                     paramGaussPoints,
                                     worksetJacobians,
                                     i, hex_8);

         // face outer normals (relative to parent cell) are uTan x vTan
            RealSpaceTools<double>::vecprod(worksetFaceN, worksetFaceTu, worksetFaceTv);


         // evaluate curl u at face points
           for(int nPt = 0; nPt < numFacePoints; nPt++){

             double x = worksetGaussPoints(0, nPt, 0);
             double y = worksetGaussPoints(0, nPt, 1);
             double z = worksetGaussPoints(0, nPt, 2);

             evalCurlu(curluFace(0,nPt,0), curluFace(0,nPt,1), curluFace(0,nPt,2), x, y, z, muVal(k));
           }

         // compute the cross product of curluFace with basis and multiply by weights
           for (int nF = 0; nF < numFieldsD; nF++){
              for(int nPt = 0; nPt < numFacePoints; nPt++){
                  worksetDataCrossField(0,nF,nPt,0) = (curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,2)
                                 - curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,1))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,1) = (curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,0)
                                 - curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,2))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,2) = (curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,1)
                                 - curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,0))
                                  *paramGaussWeights(nPt);
              } //nPt
           } //nF

          // integrate
           fst::integrate<double>(gDBoundary, worksetFaceN, worksetDataCrossField,
                             COMP_CPP);

          // apply signs
           fst::applyFieldSigns<double>(gDBoundary, hexFaceSigns);

          // add into hC term
            for (int nF = 0; nF < numFieldsD; nF++){
                gD(0,nF) = gD(0,nF) - gDBoundary(0,nF);
            }

        } // if faceOnBoundary
      } // numFaces

    // assemble into global vector
     for (int row = 0; row < numFieldsD; row++){
           int rowIndex = elemToFace(k,row);
           double val = hD(0,row)+gD(0,row);
           rhsD.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
 } // *** end element loop ***

  // Assemble over multiple processors, if necessary
   DCurl.GlobalAssemble(); DCurl.FillComplete(MassC.RowMap(),MassD.RowMap()); 
   MassC.GlobalAssemble();  MassC.FillComplete();
   MassD.GlobalAssemble();  MassD.FillComplete();
   MassG.GlobalAssemble();  MassG.FillComplete();
   StiffD.GlobalAssemble(); StiffD.FillComplete();
   rhsD.GlobalAssemble();


   // Build the inverse diagonal for MassC
   Epetra_CrsMatrix MassCinv(Copy,MassC.RowMap(),MassC.RowMap(),1);
   Epetra_Vector DiagC(MassC.RowMap());

   DiagC.PutScalar(1.0);
   MassC.Multiply(false,DiagC,DiagC);
   for(int i=0; i<DiagC.MyLength(); i++) {
     DiagC[i]=1.0/DiagC[i];
   }
   for(int i=0; i<DiagC.MyLength(); i++) {
     int CID=MassC.GCID(i);
     MassCinv.InsertGlobalValues(MassC.GRID(i),1,&(DiagC[i]),&CID);
   }
   MassCinv.FillComplete();

  // Set value to zero on diagonal that corresponds to boundary edge
   for(int i=0;i<numEdges;i++) {
     if (edgeOnBoundary(i)){
      double val=0.0;
      MassCinv.ReplaceGlobalValues(i,1,&val,&i);
     }
   }


    int numEntries;
    double *values;
    int *cols;

  // Adjust matrices and rhs due to boundary conditions
   for (int row = 0; row<numFaces; row++){
      MassD.ExtractMyRowView(row,numEntries,values,cols);
        for (int i=0; i<numEntries; i++){
           if (faceOnBoundary(cols[i])) {
             values[i]=0;
          }
       }
      StiffD.ExtractMyRowView(row,numEntries,values,cols);
        for (int i=0; i<numEntries; i++){
           if (faceOnBoundary(cols[i])) {
             values[i]=0;
          }
       }
    }
   for (int row = 0; row<numFaces; row++){
       if (faceOnBoundary(row)) {
          int rowindex = row;
          StiffD.ExtractMyRowView(row,numEntries,values,cols);
          for (int i=0; i<numEntries; i++){
             values[i]=0;
          }
          MassD.ExtractMyRowView(row,numEntries,values,cols);
          for (int i=0; i<numEntries; i++){
             values[i]=0;
          }
         rhsD[0][row]=0;
         double val = 1.0;
         StiffD.ReplaceGlobalValues(1, &rowindex, 1, &rowindex, &val);
       }
    }

/*
// Get boundary faces and apply zeros and ones to rhs and StiffD
     int numBCFaces=0;
     for (int i=0; i<numFaces; i++){
         if (faceOnBoundary(i)){
             numBCFaces++;
         }
      }
      int * BCFaces = new int [numBCFaces];
      int indbc=0;
      for (int i=0; i<numFaces; i++){
         if (faceOnBoundary(i)){
            BCFaces[indbc]=i;
            indbc++;
            rhsD[0][i]=0;
         }
      }
   //   ML_Epetra::Apply_OAZToMatrix(BCEdges, numBCEdges, StiffD);
      delete [] BCFaces;
*/

   
  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m0_matrix.dat",MassG);
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",MassC);
   EpetraExt::RowMatrixToMatlabFile("mag_m1inv_matrix.dat",MassCinv);
   EpetraExt::RowMatrixToMatlabFile("mag_m2_matrix.dat",MassD);
   EpetraExt::RowMatrixToMatlabFile("mag_k2_matrix.dat",StiffD);
   EpetraExt::RowMatrixToMatlabFile("mag_t1_matrix.dat",DCurl);
   EpetraExt::MultiVectorToMatrixMarketFile("rhs2_vector.dat",rhsD,0,0,false);

   fSignsout.close();

  // Clean up
    delete [] nodesPerElement;
    delete [] elementAttributes;
    delete [] elements;

    for (int b = 0; b < numElemBlk; b++){
       delete [] elmtNodeLinkage[b];
       delete [] elementTypes[b];
    }

 // delete mesh
 Delete_Pamgen_Mesh();

 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);
 
 return 0;
}

// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {

   // function 1
    uExact0 = exp(y+z)*(x+1.0)*(x-1.0);
    uExact1 = exp(x+z)*(y+1.0)*(y-1.0);
    uExact2 = exp(x+y)*(z+1.0)*(z-1.0);
  
 /*
   // function 2
    uExact0 = cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    uExact1 = cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    uExact2 = cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
 
 */  
 /*
   // function 3
    uExact0 = x*x-1.0;
    uExact1 = y*y-1.0;
    uExact2 = z*z-1.0;

   // function 4
    uExact0 = sin(M_PI*x);
    uExact1 = sin(M_PI*y);
    uExact2 = sin(M_PI*z);
 */  

   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
   
   // function 1
    double divu = 2.0*x*exp(y+z)+2.0*y*exp(x+z)+2.0*z*exp(x+y);

   // function 2
   //double divu = 2.0*x*cos(M_PI*y)*cos(M_PI*z) + 2.0*y*cos(M_PI*x)*cos(M_PI*z)
   //               + 2.0*z*cos(M_PI*x)*cos(M_PI*y);
   
   // function 3
   // double divu = 2.0*(x + y + z);
   
   // function 4
   // double divu = M_PI*(cos(M_PI*x)+cos(M_PI*y)+cos(M_PI*z));

   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, 
                double & x, double & y, double & z, double & mu)
 {
  
   // function 1
    double duxdy = exp(y+z)*(x+1.0)*(x-1.0);
    double duxdz = exp(y+z)*(x+1.0)*(x-1.0);
    double duydx = exp(x+z)*(y+1.0)*(y-1.0);
    double duydz = exp(x+z)*(y+1.0)*(y-1.0);
    double duzdx = exp(x+y)*(z+1.0)*(z-1.0);
    double duzdy = exp(x+y)*(z+1.0)*(z-1.0);
 

  /*
   // function 2
    double duxdy = -M_PI*sin(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    double duxdz = -M_PI*sin(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0);
    double duydx = -M_PI*sin(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    double duydz = -M_PI*sin(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double duzdx = -M_PI*sin(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double duzdy = -M_PI*sin(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0);
  */

    curlu0 = (duzdy - duydz)/mu;
    curlu1 = (duxdz - duzdx)/mu;
    curlu2 = (duydx - duxdy)/mu;
 
  /*
   // function 3 and 4
    curlu0 = 0;
    curlu1 = 0;
    curlu2 = 0;
  */
  
   return 0;
 }

// Calculates curl of the curl of exact solution u
 int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, 
                    double & x, double & y, double & z, double & mu)
{
   
   // function 1
    double dcurlu0dy = exp(x+y)*(z+1.0)*(z-1.0) - 2.0*y*exp(x+z);
    double dcurlu0dz = 2.0*z*exp(x+y) - exp(x+z)*(y+1.0)*(y-1.0); 
    double dcurlu1dx = 2.0*x*exp(y+z) - exp(x+y)*(z+1.0)*(z-1.0); 
    double dcurlu1dz = exp(y+z)*(x+1.0)*(x-1.0) - 2.0*z*exp(x+y);
    double dcurlu2dx = exp(x+z)*(y+1.0)*(y-1.0) - 2.0*x*exp(y+z);
    double dcurlu2dy = 2.0*y*exp(x+z) - exp(y+z)*(x+1.0)*(x-1.0);
                       

 /*
   // function 2
    double dcurlu0dy = -M_PI*M_PI*cos(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0)
                           + 2.0*y*M_PI*sin(M_PI*z)*cos(M_PI*x);
    double dcurlu0dz = -2.0*z*M_PI*sin(M_PI*y)*cos(M_PI*x)
                          + M_PI*M_PI*cos(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double dcurlu1dx = -2.0*x*M_PI*sin(M_PI*z)*cos(M_PI*y)
                          + M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double dcurlu1dz = -M_PI*M_PI*cos(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0)
                           + 2.0*z*M_PI*sin(M_PI*x)*cos(M_PI*y);
    double dcurlu2dx = -M_PI*M_PI*cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0)
                           + 2.0*x*M_PI*sin(M_PI*y)*cos(M_PI*z);
    double dcurlu2dy = -2.0*y*M_PI*sin(M_PI*x)*cos(M_PI*z)
                          + M_PI*M_PI*cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
 */
                       
    curlCurlu0 = (dcurlu2dy - dcurlu1dz)/mu;
    curlCurlu1 = (dcurlu0dz - dcurlu2dx)/mu;
    curlCurlu2 = (dcurlu1dx - dcurlu0dy)/mu;
 
 /*
   // function 3 and 4
    curlCurlu0 = 0.0;
    curlCurlu1 = 0.0;
    curlCurlu2 = 0.0;
 */

    return 0;
}
