#ifndef _fei_ArrayUtils_hpp_
#define _fei_ArrayUtils_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fwd.hpp"

#include <algorithm>

namespace fei {


  /** Lower bound finds the first entry in list that is not less than item.
      A binary search is used, and list is assumed to be sorted.
   */
  template<typename T>
  inline int lowerBound(const T& item,
		 const T* list,
		 int len)
  {
    //The correctness of this function is tested in
    // fei/test_utils/test_misc.cpp, in test_misc::serialtest2().

    if (len < 1) return 0;

    unsigned start = 0;
    unsigned end = len - 1;

    while(end - start > 23) {
      unsigned mid = (start + end) >> 1;
      if (list[mid] < item) start = mid;
      else end = mid;
    }

    while (start+7<=end ) {
      if ( item <= list[start+0] ) return start+0;
      if ( item <= list[start+1] ) return start+1;
      if ( item <= list[start+2] ) return start+2;
      if ( item <= list[start+3] ) return start+3;
      if ( item <= list[start+4] ) return start+4;
      if ( item <= list[start+5] ) return start+5;
      if ( item <= list[start+6] ) return start+6;
      if ( item <= list[start+7] ) return start+7;
      start+=8;
    }

    while (start<=end && item > list[start] ) 
      start++;

    return(start);
  }

  /** Binary search of a list that's assumed to be sorted.

      @param item to be searched for.
      @param list List to be searched.
      @param len Length of list.

      @return offset Offset at which item was found, or -1 if not found.
   */
  template<typename T>
  inline int binarySearch(const T& item,
        const T* list,
        int len)
    {
      int result = lowerBound(item,list,len);
      if ( result < len && item == list[result] )
	return result;
      return -1;
    }

  /** sort the specified array, and move the contents
    of the specified companions array to match the new order.
    This is an implementation of the insertion sort algorithm. */
  template<typename T>
  inline void insertion_sort_with_companions(int len, int* array, T* companions)
    {
      int i, j, index;
      T companion;

      for (i=1; i < len; i++) {
        index = array[i];
        companion = companions[i];
        j = i;
        while ((j > 0) && (array[j-1] > index))
        {
          array[j] = array[j-1];
          companions[j] = companions[j-1];
          j = j - 1;
        }
        array[j] = index;
        companions[j] = companion;
      }
    }


  /** Binary search of a list that's assumed to be sorted.

      @param item to be searched for.
      @param list List to be searched.
      @param len Length of list.
      @param insertPoint If item is not found, this is the offset into list at
      which item could be inserted while maintaining sortedness. Not referenced
      if item is found.

      @return offset Offset at which item was found, or -1 if not found.
   */
  template<typename T>
  inline int binarySearch(const T& item,
		     const T* list,
		     int len,
		     int& insertPoint)
    {
      //The correctness of this function is tested in src/utils/test_Set.C,
      //in the function test_Set::test2.
      insertPoint = lowerBound(item,list,len);
      if ( insertPoint < len && item == list[insertPoint] )
	return insertPoint;
      return -1;
    }

  /** Binary search of an std::vector that's assumed to be sorted.
   */
  template<typename T>
  inline int binarySearch(const T& item, const std::vector<T>& list, int& insertPoint)
    {
      if (list.size() == 0) {
        insertPoint = 0;
        return(-1);
      }
      return( binarySearch(item, &list[0], list.size(), insertPoint) );
    }

  /** Binary search of an std::vector that's assumed to be sorted.
   */
  template<typename T>
  inline int binarySearch(const T& item, const std::vector<T>& list)
    {
      if (list.size() == 0) return(-1);
      return( binarySearch(item, &list[0], list.size()) );
    }

  /** Perform a binary search but limit the search to a given range.
      @param item Value to be searched for.
      @param list
      @param listLength
      @param start Starting offset of search 'window'.
      @param end Ending offset of search 'window'. end should be less than 
      listLength.
      @param insertPoint
      @return offset position at which item was found. If not found, returns -1.
      (Since 0-based indexing is used, 'end' can't be greater than listLength-1.)
  */
  template<typename T>
  inline int binarySearch(const T& item, const T* list, int /*listLength*/,
		     int start, int end, int& insertPoint)
    {
      int result = binarySearch(item,list+start,end-start+1,insertPoint);
      if ( result >= 0 )
	return result+start;
      insertPoint+= start;
      return -1;
    }

   /** Perform a binary search for each item in a sorted input list.

       \param numItems   number of items to be searched for
       \param items      list of items (length numItems) to be searched for
       \param offsets    list (length numItems) allocated by caller. On exit,
                         offsets[i] contains the offset of item in 'list', or
                         -1 if item is not present in 'list'.
       \param list       array (length 'listLength') to be searched
       \param listLength length of input array 'list'
    */
  template<typename T>
  inline int binarySearch(int numItems, const T* items, int* offsets,
		     const T* list, int listLength)
    {
      int i;
      if (numItems < 1 || listLength < 1) {
        if (listLength < 1) {
          for(i=0; i<numItems; ++i) offsets[i] = -1;
        }
      }

      int tmp, start = 0;
      int end = listLength -1;
      int insertPoint = -1;
      for(i=0; i<numItems; ++i) {
        tmp = binarySearch(items[i], list, listLength, start, end, insertPoint);
        start = tmp > -1 ? tmp : insertPoint;
        offsets[i] = tmp;
      }

      return(0);
    }

  /** Insert an item into a sorted list, maintaining sortedness.
      If the item is inserted, return the offset at which it was inserted.
      If the item was already present, return -1.
   */
  template<class T>
  inline int sortedListInsert(const T& item, std::vector<T>& list)
    {
      typename std::vector<T>::iterator iter =
        std::lower_bound(list.begin(), list.end(), item);

      if (iter == list.end() || *iter != item) {
        iter = list.insert(iter, item);
        return( iter - list.begin() );
      }

      return(-1);
    }

  /** Insert an item into a sorted list, maintaining sortedness.
   */
  template<class T>
  inline int sortedListInsert(const T& item, T*& list, int& len, int& allocLen)
    {
      int i, insertPoint = -1;
      int index = binarySearch(item, list, len, insertPoint);
      if (index < 0) {
        if (len >= allocLen) {
          allocLen = len+2;
          T* newlist = new T[allocLen];
          for(i=0; i<insertPoint; ++i) newlist[i] = list[i];
          newlist[insertPoint] = item;
          for(i=insertPoint; i<len; ++i) newlist[i+1] = list[i];
          delete [] list;
          list = newlist;
        }
        else {
          for(i=len; i>insertPoint; --i) {
            list[i] = list[i-1];
          }
          list[insertPoint] = item;
        }
        ++len;
        return(insertPoint);
      }

      return(-1);
    }

  /** Insert an item into a list at a specified position.
   */
  template<class T>
  inline int listInsert(const T& item, int offset, T*& list,
		   int& usedLength, int& allocatedLength,
		   int allocChunkSize=200)
    {
      if (offset < 0 || offset > usedLength) {
        return(-1);
      }

      if (usedLength < allocatedLength) {
        for(int i=usedLength; i>offset; --i) {
          list[i] = list[i-1];
        }
        list[offset] = item;
        ++usedLength;
        return(0);
      }

      T* newlist = new T[allocatedLength+allocChunkSize];

      allocatedLength += allocChunkSize;
      int i;
      for(i=0; i<offset; ++i) {
        newlist[i] = list[i];
      }

      newlist[offset] = item;

      for(i=offset+1; i<=usedLength; ++i) {
        newlist[i] = list[i-1];
      }

      ++usedLength;
      delete [] list;
      list = newlist;
      return(0);
    }

  /** Simple exhaustive search of a list.
      @return offset at which item is found, or -1 if not found.
  */
  template<class T>
  inline int searchList(const T& item, const T* list, int len)
    {
      for(int i=0; i<len; ++i) {
        if (list[i] == item) {
          return(i);
        }
      }
      return(-1);
    }

} //namespace fei

#endif // _fei_ArrayUtils_hpp_

