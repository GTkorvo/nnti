#ifndef PHX_GRID_QUAD_H
#define PHX_GRID_QUAD_H

#include "Teuchos_TestForException.hpp"

#include "phx_grid_Element.h"
#include "phx_grid_Segment.h"

namespace phx {
namespace grid {

class Quad : public Element
{
  public:
    Quad(const int numDimensions)
    {
      TEST_FOR_EXCEPTION(numDimensions <= 1, std::out_of_range,
                         "numDimensions = " << numDimensions << ", should be > 1");

      setLabel("phx::grid::Quad");
      setNumVertices(4);
      setNumDimensions(numDimensions);
      setNumComponents(4);
      Segment component(numDimensions);
      for (int i = 0; i < 4; ++i)
        setComponent(i, component);
    }
}; 

} // namespace grid
} // namespace phx
#endif
