#include "MesquiteError.hpp"

using namespace Mesquite;

/* Special handler for the bad_alloc exception thrown by new */ 
void Mesquite::out_of_store()
{
  std::cerr << "MESQUITE ERROR: Operator new failed: out of store.\n";
  throw std::bad_alloc();
}

// Constructors
MsqError::MsqError(bool e, enum Error_Codes ec, std::string m)
{
  errorOn = e;
  printError=true;
  errorCode = ec;
  if (m!="")
    msg = m;
}

MsqError::MsqError(Error_Codes ec, std::string m)
{
  errorOn = true;
  printError=true;
  errorCode = ec;
  if (m != std::string(""))
    msg = m;
}

/*! 
The copy constructor has protected access since users should always pass
errors by reference.
The copy constructor sends a warning if used, since MsqError object
should always be passed by reference.
*/
MsqError::MsqError(const MsqError &A)
{
  std::cerr << "WARNING: the MsqError() copy constructor has been used\n"
            << "         MsqError objects should always be passed by reference";
  errorOn = A.errorOn;
  errorCode = A.errorCode;
  if (A.msg!="")
    msg = A.msg;
  printError=A.printError;
}


void MsqError::set_msg(const std::string &m)
{
  errorOn = true;
  msg = m;
}

void MsqError::set_error_code(enum Error_Codes ec)
{
  errorOn = true;
  errorCode = ec;
}

void MsqError::reset()
{
  errorOn = false;
  errorCode = MSQ_NO_ERROR;
  msg = "";
  printError=true;
}

/*! \fn MsqError::handler(int line, const char *func,const char* file, const char *dir)
  This functions provides by default the following error handling mechanism:
  - prints the contextual error message, if any was set with set_msg(...) .
  - handles a standard error code, if any was set with set_error_code(...) .
  - prints stack tracing information.
  */
#undef __FUNC__
#define __FUNC__ "MsqError::handler"
void MsqError::handler(int line, const char *func,
                       const char* file, const char *dir)
{
  
  if(printError){
    if (msg!="" && errorCode!=MSQ_PRINT_STACK) {
      std::cerr << "MESQUITE ERROR: " << msg << std::endl;
    }
    switch(errorCode){
      case MSQ_MEM_ERR:
        std::cerr << "MESQUITE ERROR:  Out of memory. \n";
        break;
      case MSQ_NULL_ERR:
        std::cerr << "MESQUITE ERROR:  Null pointer. \n";
        break;
      case MSQ_INIT_ERR:
        std::cerr << "MESQUITE ERROR:  Data Structure Not Initialized. \n";
        break;
      case MSQ_INPUT_ERR:
        std::cerr << "MESQUITE ERROR:  Incorrect Input \n";
        break;
      case MSQ_FILE_OPEN_ERR:
        std::cerr << "MESQUITE ERROR:  File open error \n";
        break;
      case MSQ_FREE_ERR:
        std::cerr << "MESQUITE ERROR:  Error freeing memory \n";
        break;
      case MSQ_INVALID_MESH_ERR:
        std::cerr << "MESQUITE ERROR:  Invalid Mesh; use SMuntangle to create a valid mesh prior to smoothing \n";
        break;
      case MSQ_DIVIDE_BY_ZERO_ERR:
        std::cerr << "MESQUITE ERROR:  Division by zero \n";
        break;
      case MSQ_DATA_ERR:
        std::cerr << "MESQUITE ERROR:  Incorrect data \n";
        break;
      case MSQ_NO_PD_STORAGE_MODE:
        std::cerr << "MESQUITE ERROR: no storage mode set in PatchData object.\n";
        break;
      case MSQ_NO_ERROR:
        break;
      default:
        break;
        
    }
    std::cerr << "MESQUITE ERROR: " << func << "()  line " << line
              << " in " << dir <<"/"<< file << std::endl;
  
#ifdef MESQUITE_PRINT_ERROR_STACK
    errorCode = MSQ_PRINT_STACK; /* set it to print the stack */
#else
    errorCode = MSQ_DO_NOT_PRINT;
    printError=false;
#endif
  }
  
  return ;
}

