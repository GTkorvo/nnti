
#ifndef EDT_CRSGRAPH_MAPCOLORING_H
#define EDT_CRSGRAPH_MAPCOLORING_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_MapColoring;

namespace EpetraExt {

class CrsGraph_MapColoring : public StructuralTransform<Epetra_CrsGraph,Epetra_MapColoring> {

  bool verbose_;

 public:

  ~CrsGraph_MapColoring() {}

  CrsGraph_MapColoring( bool verbose = false )
  : verbose_(verbose)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_MAPCOLORING_H
