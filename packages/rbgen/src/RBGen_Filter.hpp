#ifndef RBGEN_FILTER_H
#define RBGEN_FILTER_H

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace RBGen {

  enum SortType {
    LARGEST,
    SMALLEST
  };

  //! Class for selecting desired singular values
  template <class ScalarType>
  class Filter {
    public: 

    //@{ @name Constructor/Destructor.

    //! Default constructor.
    Filter(SortType which = LARGEST) {};

    //! Destructor.
    virtual ~Filter() {};
    //@}

    //@{ @name Filter Method
    virtual vector<int> filter(const vector<ScalarType> &svals) = 0;
    //@}
  };

  //! Range-based filter
  template <class ScalarType>
  class RangeFilter : public Filter<ScalarType> {
    public: 
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    RangeFilter(SortType which, int minRank = 1, int maxRank = 1) 
      : which_(which) {
      TEST_FOR_EXCEPTION(minRank < 1,logic_error,"RangeFilter: minRank must be > 1");
      TEST_FOR_EXCEPTION(maxRank < 1,logic_error,"RangeFilter: maxRank must be > 1");
      minRank_ = minRank;
      maxRank_ = maxRank;
    };

    //! Destructor.
    virtual ~RangeFilter() {};
    //@}

    //@{ @name Filter Method
    //! \brief Pass at most maxRank and at most minRank values through the filter.
    //  \note It is assumed that svals are sorted in decreasing order.
    vector<int> filter(const vector<ScalarType> &svals) {
      int n = (unsigned int)svals.size();
      if      (n > maxRank_) n = maxRank_;
      else if (minRank_ < n) n = minRank_;
      vector<int> ret(n);
      if (LARGEST == which_) {
        for (int i=0; i<n; ++i) {
          ret[i] = i;
        }
      }
      else if (SMALLEST == which_) {
        for (int i=svals.size()-1; i>svals.size()-n-1; --i) {
          ret[i] = i;
        }
      }
      return ret;
    }
    //@}

    private: 
      SortType which_;
      int minRank_;
      int maxRank_;
  };

  //! Threshold-based filter
  template <class ScalarType>
  class ThreshFilter : public Filter<ScalarType> {
    public: 
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    ThreshFilter(SortType which, bool absthresh, 
                 typename Teuchos::ScalarTraits<ScalarType>::magnitudeType thresh)
      : which_(which), absthresh_(absthresh) {
      TEST_FOR_EXCEPTION(thresh < Teuchos::ScalarTraits<ScalarType>::zero(),
                         logic_error,"ThreshFilter: minRank must be > 1");
      thresh_ = thresh;
    };

    //! Destructor.
    virtual ~ThreshFilter() {};
    //@}

    //@{ @name Filter Method
    //! \brief Return an index for the singular values falling within the
    //   threshold.
    //
    //   The threshold definition depends on \c which:
    //     - ::LARGEST: value \c cand passes the filter if \f$cand \geq tval\f$
    //     - ::SMALLEST: value \c cand passes the filter if \f$cand \leq tval\f$
    //   where \c tval depends on \c absthresh:
    //     - \c false: 
    //        - \f$tval \doteq thresh \cdot max(svals)\f$
    //        - \f$tval \doteq thresh \cdot min(svals)\f$
    //     - \c true: \f$tval \doteq thresh\f$
    //
    //  \note It is assumed that svals are sorted in decreasing order.
    vector<int> filter(const vector<ScalarType> &svals) {
      const int last = (unsigned int)svals.size() - 1;

      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tval;
      if (absthresh_) {
        tval = thresh_;
      }
      else {
        if (LARGEST == which_) {
          tval = thresh_ * svals[0];
        }
        else if (SMALLEST == which_) {
          tval = thresh_ * svals[last];
        }
      }

      vector<int> ret;
      if (LARGEST == which_) {
        int num = find(svals.begin(),svals.end(),bind2nd(less<ScalarType>,tval)) - svals.begin();
        ret.resize(num);
        for (int i=0; i<num; ++i) {
          ret[i] = i;
        }
      }
      else if (SMALLEST == which_) {
        int num = svals.end() - find(svals.begin(),svals.end(),bind2nd(less<ScalarType>,tval)) + 1;
        ret.resize(num);
        for (int i=last; i>last-n; --i) {
          ret[i] = i;
        }
      }
      return ret;
    }
    //@}

    private: 
      SortType which_;
      bool absthresh_;
      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType thresh_;
  };

  //! Composite filter
  template <class ScalarType>
  class CompFilter : public Filter<ScalarType> {
    public: 

    enum CompType {
      AND,
      OR
    };

    //@{ @name Constructor/Destructor.

    //! Default constructor.
    CompFilter(CompType andor, 
                 const Teuchos::RefCountPtr<Filter<ScalarType> > &f1,
                 const Teuchos::RefCountPtr<Filter<ScalarType> > &f2 ) 
      : andor_(andor), f1_(f1), f2_(f2) {
      TEST_FOR_EXCEPTION(f1_ == Teuchos::null || f2_ == Teuchos::null,
                         logic_error,"CompFilter: Component filters must be non-null.");
    };

    //! Destructor.
    virtual ~CompFilter() {};
    //@}

    //@{ @name Filter Method
    //! \brief 
    //  \note It is assumed that svals are sorted in decreasing order.
    vector<int> filter(const vector<ScalarType> &svals) {
      vector<int> ind1 = f1->filter(svals);
      vector<int> ind2 = f2->filter(svals);
      vector<int> ret;
      if (AND == andor_) {
        set_intersection(ind1.begin(),ind1.end(),ind2.begin(),ind2.end());
      }
      else if (OR == andor_) {
        set_union(ind1.begin(),ind1.end(),ind2.begin(),ind2.end());
      }
      return ret;
    }
    //@}

    private: 
      CompType andor_;
      Teuchos::RefCountPtr<Filter<ScalarType> > f1_, f2_;
  };

}
#endif
