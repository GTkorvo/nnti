/*Paul
04-Aug-2002 Status: Templated for class T. All Epetra methods except Distributor. Fully documented. Switched to images.
03-Sept-2002 Took out Directory and ImageID methods. Templated for PacketType, OrdinalType.
12-Oct-2002 Added some consts (still some left). Updated for Common->Compiler_Directives renaming.
30-Oct-2002 Updated for Compiler_Directives -> ConfigDefs renaming.
12-Nov-2002 Changed remaining template<class...> to template<typename...>
19-Nov-2002 myImageID and numImages moved back from Platform.
06-Feb-2003 Updated const syntax.
*/

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_COMM_HPP_
#define _TEUCHOS_COMM_HPP_

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

//! Teuchos::Comm:  The Teuchos Communication Abstract Base Class.
/*! The Teuchos Comm class is an interface that encapsulates the general
  information and services needed for other Teuchos classes to run on a parallel computer.
  
  Comm currently has one default implementation, via SerialComm, for serial execution.
  (A second default implementation is planned. MpiComm will be for MPI
  distributed memory execution.)  It is meant to insulate the user from
  the specifics of communication that are not required for normal
  manipulation of linear algebra objects.  Most Comm interfaces are similar to MPI
  interfaces, except that the type of data is not required as an argument since C++ can bind
  to the appropriate interface based on argument typing.
*/

template<typename PacketType, typename OrdinalType>
class Comm {
  public:

  //@{ \name Constructor/Destructor Methods
  //! Destructor
  virtual ~Comm() {};
  //@}

	//@{ \name Image Info Methods
	//! getMyImageID
	virtual int getMyImageID() const = 0;
	//! getNumImages
	virtual int getNumImages() const = 0;
	//@}
  
  //@{ \name Barrier Methods
  //! Barrier. Each image must stop until all images have reached the barrier.
  virtual void barrier() const = 0;
  //@}

  //@{ \name Broadcast Methods
  //! Broadcast
  /*!Take list of input values from the root image and sends to all other images.
    \param myVals InOut
           On entry, the root image contains the list of values.  On exit,
	   all images will have the same list of values.  Note that values must be
	   aalocated on all images before the broadcast.
    \param count In
           On entry, contains the length of the list of myVals.
    \param root In
           On entry, contains the imageID from which all images will receive a copy of myVals.
  */
  virtual void broadcast(PacketType* myVals, OrdinalType const count, int const root) const = 0;
  //@}

  //@{ \name Gather Methods
  //! Gather All function.
  /*! Take list of input values from all images in the communicator and creates an ordered contiguous list of
    those values on each image.
    \param myVals In
           On entry, contains the list of values to be sent to all images.
    \param allVals Out
           On exit, contains the list of values from all images. Must be of size numImages*count.
    \param count In
           On entry, contains the length of the list of values.
  */
  virtual void gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const = 0;
  //@}

  //@{ \name Sum Methods
  //! Global Sum function.
  /*!Take list of input values from all images in the communicator, computes the sum and returns the
    sum to all images.
    \param partialSums In
           On entry, contains the list of values, usually partial sums computed locally,
	   to be summed across all images.
    \param globalSums Out
           On exit, contains the list of values summed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  virtual void sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const = 0;
  //@}
	
  //@{ \name Max/Min Methods
  //! Global Max function.
  /*! Take list of input values from all images in the communicator, computes the max and returns the
    max to all images.
    \param partialMaxs In
           On entry, contains the list of values, usually partial maxs computed locally;
	   using these Partial Maxs, the max across all images will be computed.
    \param globalMaxs Out
           On exit, contains the list of maxs computed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  virtual void maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const = 0;
  //! Global Min function.
  /*! Take list of input values from all images in the communicator, computes the min and returns the
    min to all images.
    \param partialMins In
           On entry, contains the list of values, usually partial mins computed locally;
	   using these Partial Mins, the min across all images will be computed.
    \param globalMins Out
           On exit, contains the list of mins computed across all images.
    \param count In
           On entry, contains the length of the list of values.
  */
  virtual void minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const = 0;
  //@}

  //@{ \name Parallel Prefix Methods
  //! Scan Sum function.
  /*! Take list of input values from all images in the communicator, computes the scan sum and returns it 
    to all images such that image i contains the sum of values from image 0 up to and including
   image i.
    \param myVals In
           On entry, contains the list of values to be summed across all images.
    \param scanSums Out
           On exit, contains the list of values summed across images 0 through i.
    \param count In
           On entry, contains the length of the list of values.
  */
  virtual void scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const = 0;
  //@}

	//@{ \name I/O Methods
	//! printInfo
	virtual void printInfo(ostream& os) const = 0;
	//@}
	
}; // class Comm

} // namespace Teuchos

#endif // _TEUCHOS_COMM_HPP_
