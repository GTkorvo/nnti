#ifndef RBGEN_ISVD_SINGLEUDV_H
#define RBGEN_ISVD_SINGLEUDV_H

#include "RBGen_ISVDSingle.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  class ISVD_SingleUDV : public virtual ISVDUDV, public virtual ISVDSingle {
    public:
      ISVD_SingleUDV();
      virtual ~ISVD_SingleUDV() {};
  };

} // end of RBGen namespace

#endif // RBGEN_ISVD_SINGLEUDV_H
