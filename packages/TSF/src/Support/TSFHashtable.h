#ifndef TSFHASHTABLE_H
#define TSFHASHTABLE_H

#include "TSFDefs.h"
#include "TSFArray.h"
#include "TSFHashUtils.h"

namespace TSF
{
  using std::string;

  /** \ingroup Containers
   * helper class for hashtable, representing a single key, value pair.
   */
  template<class Key, class Value> class HashPair
    {
    public:
      /** empty ctor */
      inline HashPair() : key_(), value_() {;}
      /** construct with a key and value */
      inline HashPair(const Key& key, const Value& value)
        : key_(key), value_(value) {;}

      Key key_;
      Value value_;
    };

  /**
     \ingroup Containers
     Templated hashtable class.
     @author Kevin Long
  */


  template<class Key, class Value> class TSFHashtable
    {
    public:
      /**
       * Create an empty TSFHashtable
       */
      inline TSFHashtable(int capacity=101, double rehashDensity = 0.8);

      /**
       * Check for the presence of a key
       */
      inline bool containsKey(const Key& key) const ;
      /**
       * Get the value indexed by key
       */
      inline const Value& get(const Key& key) const ;
      /**
       * Put a new (key, value) pair in the table.
       */
      inline void put(const Key& key, const Value& value);

      /**
       * Remove from the table the element given by key.
       */
      inline void remove(const Key& key);

      /**
       * Get the number of elements in the table
       */
      inline int size() const {return count_;}

      /**
       * Get lists of keys and values in TSFArray form
       */
      inline void arrayify(TSFArray<Key>& keys, TSFArray<Value>& values) const ;

      /**
       * Return the average degeneracy (average number of entries per hash code).
       */
      inline double avgDegeneracy() const {return avgDegeneracy_;}

      /**
       * Return the density of the hashtable (num entries / capacity)
       */
      inline double density() const {return ((double)count_)/((double) capacity_);}


      /**
       * Set the density at which to do a rehash
       */
      inline void setRehashDensity(double rehashDensity);

    private:

      inline void rehash();
      inline int nextPrime(int newCap) const ;
      inline void accumulateAvgFill(int n) const ;

      TSFArray<TSFArray<HashPair<Key, Value> > > data_;
      int count_;
      int capacity_;
      mutable Value mostRecentValue_;
      mutable Key mostRecentKey_;

      mutable int nHits_;
      mutable double avgDegeneracy_;
      double rehashDensity_;
    };

  /*
    template<class Key, class Value>
    inline string toString(const TSFHashtable<Key, Value>& h);
  */

  template<class Key, class Value>
    ostream& operator<<(ostream& os, const TSFHashtable<Key, Value>& h);

  template<class Key, class Value> inline
    TSFHashtable<Key, Value>::TSFHashtable(int capacity, double rehashDensity)
    : data_(), count_(0), capacity_(TSFHashUtils::nextPrime(capacity)),
    nHits_(0), avgDegeneracy_(0), rehashDensity_(rehashDensity)
    {
      data_.resize(capacity_);
    }

  template<class Key, class Value> inline
    bool TSFHashtable<Key, Value>::containsKey(const Key& key) const
    {
      try
        {
          const TSFArray<HashPair<Key, Value> >& candidates
            = data_[hashCode(key) % capacity_];

          accumulateAvgFill(candidates.length());

          for (int i=0; i<candidates.length(); i++)
            {
              const HashPair<Key, Value>& c = candidates[i];
              if (c.key_ == key)
                {
                  //          (Key&) mostRecentKey_ = key;
                  //(Value&) mostRecentValue_ = c.value_;
                  return true;
                }
            }
        }
      catch(exception& e)
        {
          TSFError::trace(e, "in TSFHashtable::containsKey()");
        }
      return false;
    }

  template<class Key, class Value> inline
    void TSFHashtable<Key, Value>::put(const Key& key, const Value& value)
    {
      try
        {
          int index = hashCode(key) % capacity_;

          TSFArray<HashPair<Key, Value> >& local = data_[index];

          // check for duplicate key
          for (int i=0; i<local.length(); i++)
            {
              if (local[i].key_ == key)
                {
                  local[i].value_ = value;
                  return;
                }
            }

          // no duplicate key, so increment element count by one.
          count_++;

          // check for need to resize.
          if ((double) count_ > rehashDensity_ * (double) capacity_)
            {
              capacity_ = TSFHashUtils::nextPrime(capacity_+1);
              rehash();
              // recaluate index
              index = hashCode(key) % capacity_;
            }

          data_[index].append(HashPair<Key, Value>(key, value));
        }
      catch(exception& e)
        {
          TSFError::trace(e, "in TSFHashtable::put()");
        }
    }



  template<class Key, class Value> inline
    void TSFHashtable<Key, Value>::rehash()
    {
      try
        {
          TSFArray<TSFArray<HashPair<Key, Value> > > tmp(capacity_);

          for (int i=0; i<data_.length(); i++)
            {
              for (int j=0; j<data_[i].length(); j++)
                {
                  int newIndex = hashCode(data_[i][j].key_) % capacity_;
                  tmp[newIndex].append(data_[i][j]);
                }
            }

          data_ = tmp;
        }
      catch(exception& e)
        {
          TSFError::trace(e, "in TSFHashtable::rehash()");
        }
    }

  template<class Key, class Value> inline
    void TSFHashtable<Key, Value>::arrayify(TSFArray<Key>& keys, TSFArray<Value>& values) const
    {
      keys.reserve(size());
      values.reserve(size());

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              keys.append(data_[i][j].key_);
              values.append(data_[i][j].value_);
            }
        }
    }

  template<class Key, class Value>  inline
    string toString(const TSFHashtable<Key, Value>& h)
    {
      TSFArray<Key> keys;
      TSFArray<Value> values;
      h.arrayify(keys, values);

      string rtn = "[";
      for (int i=0; i<keys.length(); i++)
        {
          rtn += "{" + toString(keys[i]) + ", " + toString(values[i])
            + "}";
          if (i < keys.length()-1) rtn += ", ";
        }
      rtn += "]";

      return rtn;
    }

  template<class Key, class Value> inline
    const Value& TSFHashtable<Key, Value>::get(const Key& key) const
    {
      try
        {
          if (!containsKey(key)) TSFError::raise(toString(key)
                                                 + "key not found in TSFHashtable"
                                                 + toString(*this));

          const TSFArray<HashPair<Key, Value> >& candidates
            = data_[hashCode(key) % capacity_];

          accumulateAvgFill(candidates.length());

          for (int i=0; i<candidates.length(); i++)
            {
              const HashPair<Key, Value>& c = candidates[i];
              if (c.key_ == key)
                {
                  return c.value_;
                }
            }
        }
      catch(exception& e)
        {
          TSFError::trace(e, "in TSFHashtable::get()");
        }
      return mostRecentValue_;
    }

  template<class Key, class Value> inline
    void TSFHashtable<Key, Value>::remove(const Key& key)
    {
      if (!containsKey(key)) TSFError::raise(toString(key)
                                             + " not found in TSFHashtable" );

      count_--;
      int h = hashCode(key) % capacity_;
      const TSFArray<HashPair<Key, Value> >& candidates = data_[h];

      for (int i=0; i<candidates.length(); i++)
        {
          const HashPair<Key, Value>& c = candidates[i];
          if (c.key_ == key)
            {
              data_[h].remove(i);
              break;
            }
        }
    }

  template<class Key, class Value> inline
    void TSFHashtable<Key, Value>::accumulateAvgFill(int n) const
    {
      avgDegeneracy_ = ((double) nHits_)/(nHits_ + 1.0) * avgDegeneracy_ + ((double) n)/(nHits_ + 1.0);
      nHits_++;
    }

  template<class Key, class Value>  inline
    ostream& operator<<(ostream& os, const TSFHashtable<Key, Value>& h)
    {
      return os << toString(h);
    }


}

#endif
