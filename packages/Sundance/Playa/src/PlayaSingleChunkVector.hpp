/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SINGLECHUNKVECTOR_HPP
#define PLAYA_SINGLECHUNKVECTOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorBaseDecl.hpp"

extern "C"
{
void daxpy_(int*, double*, double*, int*, double*, int*);
}


namespace Playa
{
/**
 * Base class for vector types that have all on-processor data in a single
 * contiguous chunk
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class SingleChunkVector : public VectorBase<Scalar>
{
public:
  /** */
  SingleChunkVector() : rewound_(true) {}
  /** virtual dtor */
  virtual ~SingleChunkVector() {}

  /** */
  virtual ConstDataChunk<Scalar> nextConstChunk() const 
    {rewound_ = false; return ConstDataChunk<Scalar>(chunkSize(), dataPtr());}

  /** */
  virtual NonConstDataChunk<Scalar> nextChunk() 
    {rewound_ = false; return NonConstDataChunk<Scalar>(chunkSize(), dataPtr());}

  /** */
  virtual bool hasMoreChunks() const 
    {return rewound_;}

  /** */
  virtual void rewind() const 
    {rewound_=true;}

  /** \name Access to local elements */
  //@{
  /** read the element at the given local index */
  virtual const double& operator[](int localIndex) const 
    {return dataPtr()[localIndex];}

  /** writable access to the element at the given local index */
  virtual double& operator[](int localIndex) 
    {return dataPtr()[localIndex];}
  //@}

  /** */
  virtual int chunkSize() const = 0 ;
  /** */
  virtual const Scalar* dataPtr() const = 0 ;
  /** */
  virtual Scalar* dataPtr() = 0 ;

  virtual void update(const Scalar& alpha, const VectorBase<Scalar>* other,
    const Scalar& gamma)
    {
      const SingleChunkVector<Scalar>* sco 
        = dynamic_cast<const SingleChunkVector<Scalar>* >(other);
      TEST_FOR_EXCEPT(sco==0);

      Scalar* const myVals = this->dataPtr();
      const Scalar* const yourVals = sco->dataPtr();
      int n = chunkSize();
      int stride = 1;
      if (gamma==1.0)
      {
        daxpy_(&n, (double*) &alpha, (double*) yourVals, &stride, 
          myVals, &stride);
      }
      else
      {
        for (int i=0; i<n; i++)
        {
          myVals[i] = gamma*myVals[i] + alpha*yourVals[i];
        }
      }
    }

  /** */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma) 
    {
      const SingleChunkVector<Scalar>* scx 
        = dynamic_cast<const SingleChunkVector<Scalar>* >(x);
      TEST_FOR_EXCEPT(scx==0);
      const SingleChunkVector<Scalar>* scy 
        = dynamic_cast<const SingleChunkVector<Scalar>* >(y);
      TEST_FOR_EXCEPT(scy==0);

      Scalar* const myVals = this->dataPtr();
      const Scalar* const xVals = scx->dataPtr();
      const Scalar* const yVals = scy->dataPtr();
      int n = chunkSize();
      if (gamma==1.0)
      {
        for (int i=0; i<n; i++)
        {
          myVals[i] += alpha*xVals[i] + beta*yVals[i];
        }
      }
      else
      {
        for (int i=0; i<n; i++)
        {
          myVals[i] = gamma*myVals[i] + alpha*xVals[i] + beta*yVals[i];
        }
      }
    }

  /** */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma, const VectorBase<Scalar>* z,
    const Scalar& delta) 
    {
      const SingleChunkVector<Scalar>* scx 
        = dynamic_cast<const SingleChunkVector<Scalar>* >(x);
      TEST_FOR_EXCEPT(scx==0);
      const SingleChunkVector<Scalar>* scy 
        = dynamic_cast<const SingleChunkVector<Scalar>* >(y);
      TEST_FOR_EXCEPT(scy==0);
      const SingleChunkVector<Scalar>* scz
        = dynamic_cast<const SingleChunkVector<Scalar>* >(z);
      TEST_FOR_EXCEPT(scz==0);

      Scalar* const myVals = this->dataPtr();
      const Scalar* const xVals = scx->dataPtr();
      const Scalar* const yVals = scy->dataPtr();
      const Scalar* const zVals = scz->dataPtr();

      int n = chunkSize();
      if (delta==1.0)
      {
        for (int i=0; i<n; i++)
        {
          myVals[i] += alpha*xVals[i] + beta*yVals[i] + gamma*zVals[i];
        }
      }
      else
      {
        for (int i=0; i<n; i++)
        {
          myVals[i] = delta*myVals[i] + alpha*xVals[i] + beta*yVals[i] 
            + gamma*zVals[i];
        }
      }
    }

private:
  mutable bool rewound_;
};
}


#endif
