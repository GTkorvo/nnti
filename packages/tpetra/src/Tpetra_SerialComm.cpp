namespace Tpetra {

  // default constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm() 
    : Object("Tpetra::Comm[Serial]") {}
  
  // copy constructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::SerialComm(SerialComm<PacketType, OrdinalType> const& comm) 
    : Object(comm.label()) {}
  
  // destructor
  template<typename PacketType, typename OrdinalType>
  SerialComm<PacketType, OrdinalType>::~SerialComm() {}
  
  // gather all
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::gatherAll(PacketType* myVals, PacketType* allVals, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      allVals[i] = myVals[i];
  }
  
  // sum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::sumAll(PacketType* partialSums, PacketType* globalSums, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalSums[i] = partialSums[i];
  }
  
  // min/max
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::maxAll(PacketType* partialMaxs, PacketType* globalMaxs, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalMaxs[i] = partialMaxs[i];
  }
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::minAll(PacketType* partialMins, PacketType* globalMins, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      globalMins[i] = partialMins[i];
  }
  
  // scansum
  template<typename PacketType, typename OrdinalType>
  void SerialComm<PacketType, OrdinalType>::scanSum(PacketType* myVals, PacketType* scanSums, OrdinalType const count) const {
    for(OrdinalType i = 0; i < count; i++)
      scanSums[i] = myVals[i];
  }

	// print

  
} // namespace Tpetra
