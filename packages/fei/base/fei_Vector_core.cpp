/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fstream.hpp"

#include "fei_Vector_core.hpp"
#include "fei_VectorSpace.hpp"
#include "fei_Reducer.hpp"
#include "fei_SSVec.hpp"
#include "fei_SSMat.hpp"
#include "snl_fei_RecordCollection.hpp"
#include "snl_fei_VectorTraits_SSVec.hpp"
#include "snl_fei_CommUtils.hpp"
#include "fei_TemplateUtils.hpp"

#undef fei_file
#define fei_file "fei_Vector_core.cpp"

#include "fei_ErrMacros.hpp"

fei::Vector_core::Vector_core(fei::SharedPtr<fei::VectorSpace> vecSpace,
                              int numLocalEqns)
  : eqnComm_(),
    vecSpace_(vecSpace),
    intCommUtils_(),
    firstLocalOffset_(0),
    lastLocalOffset_(0),
    numLocal_(0),
    work_indices_(),
    work_indices2_(),
    haveFEVector_(false),
    remotelyOwned_(),
    overlapAlreadySet_(false),
    dbgprefix_("Vcore: ")
{
  intCommUtils_ = vecSpace->getCommUtils();

  eqnComm_.reset(new fei::EqnComm(intCommUtils_->getCommunicator(),numLocalEqns));
  remotelyOwned_.resize(intCommUtils_->numProcs());
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    remotelyOwned_[i] = new SSVec;
  }

  const std::vector<int>& offsets = eqnComm_->getGlobalOffsets();
  firstLocalOffset_ = offsets[intCommUtils_->localProc()];
  lastLocalOffset_ = offsets[intCommUtils_->localProc()+1] - 1;
  numLocal_ = lastLocalOffset_ - firstLocalOffset_ + 1;

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os<<dbgprefix_<<" ctor firstLocal="<<firstLocalOffset_<<", lastLocal="
     <<lastLocalOffset_<<FEI_ENDL;
  }
}

fei::Vector_core::~Vector_core()
{
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    delete remotelyOwned_[i];
  }
}

void fei::Vector_core::setOverlap(int numRemoteEqns,
                                  const int* remoteEqns)
{
  if (numRemoteEqns == 0 && remoteEqns == NULL) {
    if (overlapAlreadySet_) return;
  }

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"setOverlap"<<FEI_ENDL;
  }

  int localProc = intCommUtils_->localProc();

  if (numRemoteEqns != 0 && remoteEqns != NULL) {
    for(int i=0; i<numRemoteEqns; ++i) {
      int proc = eqnComm_->getOwnerProc(remoteEqns[i]);
      if (proc == localProc) continue;

      remotelyOwned_[proc]->addEntry(remoteEqns[i], 0.0);
    }
  }
  else {
    int numEqns = vecSpace_->getNumIndices_SharedAndOwned();
    std::vector<int> eqns(numEqns);
    vecSpace_->getIndices_SharedAndOwned(numEqns, &eqns[0], numEqns);

    for(int i=0; i<numEqns; ++i) {
      int proc = eqnComm_->getOwnerProc(eqns[i]);
      if (proc == localProc) continue;

      remotelyOwned_[proc]->addEntry(eqns[i], 0.0);
    }
  }

  overlapAlreadySet_ = true;
}

int fei::Vector_core::scatterToOverlap()
{
  if (intCommUtils_->numProcs() == 1 || haveFEVector()) {
    return(0);
  }

#ifndef FEI_SER
  if (!overlapAlreadySet_) {
    setOverlap();
  }

  //...and now the overlap is whatever is in our remotelyOwned_ vectors.

  //first find out which procs we'll be receiving from.
  std::vector<int> recvProcs;
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    if ((int)i == intCommUtils_->localProc()) continue;
    if (remotelyOwned_[i]->length() == 0) continue;

    recvProcs.push_back((int)i);
  }

  //find out the send-procs.
  std::vector<int> sendProcs;
  intCommUtils_->mirrorProcs(recvProcs, sendProcs);

  //declare arrays to send from, and corresponding sizes
  std::vector<std::vector<int> > send_ints(sendProcs.size());
  std::vector<std::vector<double> > send_doubles(sendProcs.size());
  std::vector<int> send_sizes(sendProcs.size());

  std::vector<MPI_Request> mpiReqs(sendProcs.size()+recvProcs.size());
  std::vector<MPI_Status> mpiStatuses(sendProcs.size()+recvProcs.size());
  int tag1 = 20070430;
  int tag2 = 20070460;

  //first, the procs we're going to send to, have to let us know
  //how much data we're supposed to send. So we have to receive
  //sizes and then indices from the "send"-procs.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    MPI_Irecv(&send_sizes[i], 1, MPI_INT, sendProcs[i],
              tag1, intCommUtils_->getCommunicator(), &mpiReqs[i]);
  }

  //now we'll send the sizes of our remotely-owned data to the
  //procs that we will be receiving the data from, and also the
  //indices that we want to receive.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];

    int size = remotelyOwned_[proc]->length();
    MPI_Send(&size, 1, MPI_INT, proc, tag1, intCommUtils_->getCommunicator());
  }
 
  MPI_Waitall(sendProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

  //now resize our send_ints and send_doubles arrays, and post the recvs
  //for indices that we're supposed to pack.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];
    int size = send_sizes[i];
    send_ints[i].resize(size);
    MPI_Irecv(&(send_ints[i][0]), size, MPI_INT, proc, tag1,
              intCommUtils_->getCommunicator(), &mpiReqs[i]);
    send_doubles[i].resize(size);
  }

  //now send the indices that we want to receive data for.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    int size = remotelyOwned_[proc]->indices().length();
    MPI_Send(remotelyOwned_[proc]->indices().dataPtr(), size, MPI_INT,
             proc, tag1, intCommUtils_->getCommunicator());
  }

  MPI_Waitall(sendProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

  //now post our recvs.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    int size = remotelyOwned_[proc]->indices().length();
    MPI_Irecv(remotelyOwned_[proc]->coefs().dataPtr(), size, MPI_DOUBLE,
              proc, tag2, intCommUtils_->getCommunicator(), &mpiReqs[i]);
  }

  //now pack and send the coefs that the other procs need from us.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    int num = send_sizes[i];
    int err = copyOutOfUnderlyingVector(num, &(send_ints[i][0]),
                                        &(send_doubles[i][0]), 0);
    if (err != 0) {
      FEI_COUT << "fei::Vector_core::scatterToOverlap ERROR getting data to send."<<FEI_ENDL;
      return(err);
    }

    MPI_Send(&(send_doubles[i][0]), num, MPI_DOUBLE, proc,
             tag2, intCommUtils_->getCommunicator());
  }

  MPI_Waitall(recvProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

#endif  //#ifndef FEI_SER

  return(0);
}

int fei::Vector_core::copyOut(int numValues,
				  const int* indices,
				  double* values,
				  int vectorIndex) const
{
  const std::vector<SSVec*>& remote = remotelyOwned();

  for(int i=0; i<numValues; ++i) {
    int ind = indices[i];

    int local = ind - firstLocalOffset_;
    if (local < 0 || local >= numLocal_) {
      if (ind < 0) {
	continue;
      }

      int proc = eqnComm_->getOwnerProc(ind);

      int insertPoint = -1;
      int idx = snl_fei::binarySearch(ind, remote[proc]->indices(), insertPoint);
      if (idx < 0) {
	FEI_CERR << "fei::Vector_core::copyOut: proc " << intCommUtils_->localProc()
	     << ", index " << ind << " not in remotelyOwned_ vec object."<<FEI_ENDL;
	ERReturn(-1);
      }
      else {
	values[i] = remote[proc]->coefs()[idx];
      }
    }
    else {
      CHK_ERR( copyOutOfUnderlyingVector(1, &ind, &(values[i]), vectorIndex) );
    }
  }

  return(0);
}

int fei::Vector_core::giveToVector(int numValues,
				       const int* indices,
				       const double* values,
				       bool sumInto,
				       int vectorIndex)
{
  std::vector<SSVec*>& remote = remotelyOwned();

  for(int i=0; i<numValues; ++i) {
    int ind = indices[i];
    double val = values[i];

    if (ind < 0) {
//      throw std::runtime_error("negative index not allowed");
      //preservation of existing behavior: certain Sierra scenarios
      //involve passing negative indices for positions that should be
      //ignored... so we'll continue rather than throwing.
      continue;
    }
    int proc = eqnComm_->getOwnerProc(ind);
    int local = ind - firstLocalOffset_;
    if (local < 0 || local >= numLocal_) {
      if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
        FEI_OSTREAM& os = *output_stream_;
        os << dbgprefix_<<"giveToVector remote["<<proc<<"]("
         <<ind<<","<<val<<")"<<FEI_ENDL;
      }

      if (sumInto) {
	CHK_ERR( remote[proc]->addEntry(ind, val) );
      }
      else {
	CHK_ERR( remote[proc]->putEntry(ind, val) );
      }
    }
    else {
      int err = giveToUnderlyingVector(1, &ind, &val, sumInto, vectorIndex);
      if (err != 0) {
	FEI_CERR << "giveToVector sumIn ERROR, ind: " << ind
	     << ", val: " << val << FEI_ENDL;
	ERReturn(-1);
      }
    }
  }

  return(0);
}

int fei::Vector_core::assembleFieldData(int fieldID,
					    int idType,
					    int numIDs,
					    const int* IDs,
					    const double* data,
					    bool sumInto,
					    int vectorIndex)
{
  if (vecSpace_.get() == NULL) ERReturn(-1);

  int fieldSize = vecSpace_->getFieldSize(fieldID);

  work_indices_.resize(numIDs*fieldSize);
  int* indicesPtr = &work_indices_[0];

  CHK_ERR( vecSpace_->getGlobalIndices(numIDs, IDs, idType, fieldID,
					indicesPtr) );

  CHK_ERR( giveToVector(numIDs*fieldSize, indicesPtr, data, sumInto, vectorIndex) );

  return(0);
}

int fei::Vector_core::gatherFromOverlap(bool accumulate)
{
  if (intCommUtils_->numProcs() == 1 || haveFEVector()) {
    return(0);
  }

#ifndef FEI_SER
  //first create the list of procs we'll be sending to.
  std::vector<int> sendProcs;
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    if ((int)i == intCommUtils_->localProc()) continue;
    if (remotelyOwned_[i]->length() == 0) continue;

    sendProcs.push_back(i);
  }

  std::vector<int> recvProcs;
  intCommUtils_->mirrorProcs(sendProcs, recvProcs);

  //declare arrays to hold the indices and coefs we'll be receiving.
  std::vector<std::vector<int> > recv_ints(recvProcs.size());
  std::vector<std::vector<double> > recv_doubles(recvProcs.size());
  std::vector<int> recv_sizes(recvProcs.size());

  std::vector<MPI_Request> mpiReqs(recvProcs.size()*2);
  std::vector<MPI_Status> mpiStatuses(recvProcs.size()*2);
  int tag1 = 20070430;
  int tag2 = 20070460;

  //post the recvs for the sizes.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    MPI_Irecv(&recv_sizes[i], 1, MPI_INT, proc,
              tag1, intCommUtils_->getCommunicator(), &mpiReqs[i]);
  }

  //send the sizes of data we'll be sending.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];
    int size = remotelyOwned_[proc]->length();
    MPI_Send(&size, 1, MPI_INT, proc, tag1, intCommUtils_->getCommunicator());
  }

  MPI_Waitall(recvProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

  //now post the recvs for the data.
  unsigned offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    int size = recv_sizes[i];
    std::vector<int>& recv_ints_i = recv_ints[i];
    std::vector<double>& recv_doubles_i = recv_doubles[i];
    recv_ints_i.resize(size);
    recv_doubles_i.resize(size);
    MPI_Irecv(&(recv_ints_i[0]), size, MPI_INT, proc,
              tag1, intCommUtils_->getCommunicator(), &mpiReqs[offset++]);
    MPI_Irecv(&(recv_doubles_i[0]), size, MPI_DOUBLE, proc,
              tag2, intCommUtils_->getCommunicator(), &mpiReqs[offset++]);
  }

  //now send the outgoing data.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];
    int size = remotelyOwned_[proc]->length();
    MPI_Send(remotelyOwned_[proc]->indices().dataPtr(), size, MPI_INT,
             proc, tag1, intCommUtils_->getCommunicator());
    MPI_Send(remotelyOwned_[proc]->coefs().dataPtr(), size, MPI_DOUBLE,
             proc, tag2, intCommUtils_->getCommunicator());

    remotelyOwned_[proc]->coefs() = 0.0;
  }

  MPI_Waitall(recvProcs.size()*2, &mpiReqs[0], &mpiStatuses[0]);

  //now store the data we've received.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int num = recv_sizes[i];
    std::vector<int>& recv_ints_i = recv_ints[i];
    std::vector<double>& recv_doubles_i = recv_doubles[i];
    int err = giveToUnderlyingVector(num, &(recv_ints_i[0]),
                                     &(recv_doubles_i[0]), accumulate, 0);
    if (err != 0) {
      FEI_COUT << "fei::Vector_core::gatherFromOverlap ERROR storing recvd data" << FEI_ENDL;
      return(err);
    }
  }

#endif  //#ifndef FEI_SER

  return(0);
}

int fei::Vector_core::copyOutFieldData(int fieldID,
					   int idType,
					   int numIDs,
					   const int* IDs,
					   double* data,
					   int vectorIndex)
{
  if (vecSpace_.get() == NULL) ERReturn(-1);

  int fieldSize = vecSpace_->getFieldSize(fieldID);

  if (haveFEVector_) {
    snl_fei::RecordCollection* collection = NULL;
    CHK_ERR( vecSpace_->getRecordCollection(idType, collection) );
    int nodeNumber;
    int dofOffset;
    int numInstances;
    int foffset;
    std::vector<int>& eqnNums = vecSpace_->getEqnNumbers();
    int* vspcEqnPtr = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

    int offset = 0;
    for(int i=0; i<numIDs; ++i) {
      fei::Record* node = collection->getRecordWithID(IDs[i]);
      nodeNumber = node->getNumber();
      int* eqnNumbers = vspcEqnPtr+node->getOffsetIntoEqnNumbers();
      node->getFieldMask()->getFieldEqnOffset(fieldID, foffset, numInstances);
      dofOffset = eqnNumbers[foffset] - eqnNumbers[0];
      for(int j=0; j<fieldSize; ++j) {
	CHK_ERR( copyOut_FE(nodeNumber, dofOffset+j, data[offset++]));
      }
    }
  }
  else {
    work_indices_.resize(numIDs*fieldSize*2);
    int* indicesPtr = &work_indices_[0];

    CHK_ERR( vecSpace_->getGlobalIndices(numIDs, IDs, idType,
					 fieldID, indicesPtr) );

    CHK_ERR( copyOut(numIDs*fieldSize, indicesPtr, data) );
  }

  return(0);
}

int fei::Vector_core::writeToFile(const char* filename,
                                    bool matrixMarketFormat)
{
  int numProcs = intCommUtils_->numProcs();
  int localProc =intCommUtils_->localProc();

  double coef;

  static char mmbanner[] = "%%MatrixMarket matrix array real general";

  for(int p=0; p<numProcs; ++p) {
    intCommUtils_->Barrier();
    if (p != localProc) continue;

    FEI_OFSTREAM* outFile = NULL;
    if (p==0) {
      outFile = new FEI_OFSTREAM(filename, IOS_OUT);
      FEI_OFSTREAM& ofref = *outFile;
      if (matrixMarketFormat) {
        ofref << mmbanner << FEI_ENDL;
        ofref << eqnComm_->getGlobalOffsets()[numProcs] << " 1" << FEI_ENDL;
      }
      else {
        ofref << eqnComm_->getGlobalOffsets()[numProcs] << FEI_ENDL;
      }
    }
    else outFile = new FEI_OFSTREAM(filename, IOS_APP);
    FEI_OFSTREAM& ofref = *outFile;
    ofref.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    ofref.precision(13);

    for(int i=firstLocalOffset_; i<=lastLocalOffset_; ++i) {
      CHK_ERR( copyOut(1, &i, &coef) );
      if (matrixMarketFormat) {
        ofref << " " << coef << FEI_ENDL;
      }
      else {
        ofref << i << " " << coef << FEI_ENDL;
      }
    }
    
    delete outFile;
  } 
    
  return(0);
}

int fei::Vector_core::writeToStream(FEI_OSTREAM& ostrm,
					bool matrixMarketFormat)
{
  int numProcs = intCommUtils_->numProcs();
  int localProc =intCommUtils_->localProc();

  double coef;

  static char mmbanner[] = "%%MatrixMarket matrix array real general";

  IOS_FMTFLAGS oldf = ostrm.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
  ostrm.precision(13);

  for(int proc=0; proc<numProcs; ++proc) {
    intCommUtils_->Barrier();
    if (proc != localProc) continue;

    if (proc==0) {
      if (matrixMarketFormat) {
	ostrm << mmbanner << FEI_ENDL;
	ostrm << eqnComm_->getGlobalOffsets()[numProcs] << " 1" << FEI_ENDL;
      }
      else {
	ostrm << eqnComm_->getGlobalOffsets()[numProcs] << FEI_ENDL;
      }
    }

    for(int p=0; p<localProc; ++p) {
      for(int ii=0; ii<remotelyOwned_[p]->length(); ++ii) {
	if (matrixMarketFormat) {
	  ostrm << " " << remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
	else {
	  ostrm << " " << remotelyOwned_[p]->indices()[ii] << " "
		<< remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
      }
    }

    for(int i=firstLocalOffset_; i<=lastLocalOffset_; ++i) {
      CHK_ERR( copyOut(1, &i, &coef) );
      if (matrixMarketFormat) {
	ostrm << " " << coef << FEI_ENDL;
      }
      else {
	ostrm << " " << i << " " << coef << FEI_ENDL;
      }
    }

    for(int p=localProc+1; p<numProcs; ++p) {
      for(int ii=0; ii<remotelyOwned_[p]->length(); ++ii) {
	if (matrixMarketFormat) {
	  ostrm << " " << remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
	else {
	  ostrm << " " << remotelyOwned_[p]->indices()[ii] << " "
		<< remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
      }
    }
  }

  ostrm.setf(oldf, IOS_FLOATFIELD);

  return(0);
}

