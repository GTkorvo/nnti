#ifndef TEUCHOS_ARRAY_H
#define TEUCHOS_ARRAY_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos
{

  /**
   * Array is a templated array class derived from the STL vector, but with
   * index boundschecking and an extended interface.
   */
  template<class T>
  class Array : public std::vector<T>
  {
  public:
    /** Empty ctor */
    Array();

    /** Allocate an array with n elements */
    Array(int n);

    /** Allocate n elements, and fill with value */
    Array(int n, const T& t);

    /** Stick a new entry at the end of the array. Resize to allow
     * space for the new entry
     */
    inline Array<T>& append(const T& entry) {push_back(entry); return *this;}

    /**
     * Return number of elements in the array. Equivalent to size(), but
     * included for backwards compatibiity.
     */
    int length() const {return size();}

    /**
     * Read/Write access to a the i-th element, with optional bounsdchecking.
     */
    inline T& operator[](int i);
    /**
     * Read-only access to a the i-th element, with optional boundschecking.
     */
    inline const T& operator[](int i) const;

    /** Write as a string */
    std::string toString() const ;

    /** 
     * indicate whether Array has been compiled with boundschecking on 
     */
    static bool hasBoundsChecking();

  private:

    /** check for a bounds violation if HAVE_ARRAY_BOUNDSCHECK has been
     * defined as 1. */
    void indexCheckCrash(int i) const;
  };

  /** \relates Array */
  template<class T> std::ostream& operator<<(std::ostream& os, 
                                             const Array<T>& array);

  /** \relates Array */
  template<class T> int hashCode(const Array<T>& array);

  /** \relates Array */
  template<class T> std::string toString(const Array<T>& array);


  template<class T> inline Array<T>::Array()
    : std::vector<T>()
  {}

  template<class T> inline Array<T>::Array(int n)
    : std::vector<T>(n)
  {}

  template<class T> inline Array<T>::Array(int n, const T& t)
    : std::vector<T>(n, t)
  {}

  template<class T> inline
  T& Array<T>::operator[](int i) {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    indexCheckCrash(i);
#endif
    return std::vector<T>::operator[](i);
  }

  template<class T> inline
  const T& Array<T>::operator[](int i) const {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    indexCheckCrash(i);
#endif
    return std::vector<T>::operator[](i);
  }

  template<class T> inline
  bool Array<T>::hasBoundsChecking()
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK  
    return true;
#else
    return false;
#endif
  }

  template<class T> inline
  void Array<T>::indexCheckCrash(int i) const
  {
    TEST_FOR_EXCEPTION(i<0 || i>=size(), std::range_error,
                       "Array<T>::indexCheckCrash: "
                       "index " << i << "out of range [0, "<< size() << ")");
      
  }

  // print in form (), (1), or (1,2)
  template<class T> inline ostream& operator<<(ostream& os, const Array<T>& array)
  {
    return os << toString(array);
  }

  template<class T> inline int hashCode(const Array<T>& array)
  {
    int rtn = hashCode(array.length());
    for (int i=0; i<array.length(); i++)
      {
        rtn += hashCode(array[i]);
      }
    return rtn;
  }

  template<class T> inline std::string Array<T>::toString() const
  {
    std::string rtn = "{";

    for (int i=0; i<length(); i++)
      {
        rtn += Teuchos::toString(operator[](i));
        if (i<length()-1) rtn += ", ";
      }
    rtn += "}";

    return rtn;
  }

  template<class T> inline std::string toString(const Array<T>& array)
  {
    return array.toString();
  }


  /** \relates Array create an array with one entry */
  template<class T> inline
  Array<T> tuple(const T& a)
  {
    Array<T> rtn(1, a);
    return rtn;
  }

  /** \relates Array create an array with two entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b)
  {
    Array<T> rtn(2);
    rtn[0] = a;
    rtn[1] = b;
    return rtn;
  }

  /** \relates Array create an array with three entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c)
  {
    Array<T> rtn(3);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    return rtn;
  }

  /** \relates Array create an array with four entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d)
  {
    Array<T> rtn(4);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    return rtn;
  }

  /** \relates Array create an array with five entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
  {
    Array<T> rtn(5);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    return rtn;
  }


  /** \relates Array create an array with six entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f)
  {
    Array<T> rtn(6);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    return rtn;
  }

  /** \relates Array create an array with seven entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g)
  {
    Array<T> rtn(7);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    return rtn;
  }

  /** \relates Array create an array with eight entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h)
  {
    Array<T> rtn(8);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    return rtn;
  }

  /** \relates Array create an array with nine entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h, const T& i)
  {
    Array<T> rtn(9);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    rtn[8] = i;
    return rtn;
  }


  /** \relates Array create an array with ten entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h, const T& i, const T& j)
  {
    Array<T> rtn(10);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    rtn[8] = i;
    rtn[9] = j;
    return rtn;
  }
}

#endif // TEUCHOS_ARRAY_H

