// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_FANCY_O_STREAM_HPP
#define TEUCHOS_FANCY_O_STREAM_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"

namespace Teuchos {


/** \brief Stream buffering class that performs the magic of indenting
 * data sent to an std::ostream object.
 *
 * \ingroup teuchos_outputting_grp
 *
 * Note, this is not a user-level class.  Users should use
 * <tt>basic_FancyOStream</tt>.
 */
template<typename CharT, typename Traits>
class basic_FancyOStream_buf : public std::basic_streambuf<CharT,Traits>
{
public:
  
  /** \brief . */
  typedef CharT                             char_type;
  /** \brief . */
  typedef Traits  					                traits_type;
  /** \brief . */
  typedef typename traits_type::int_type 		int_type;
  /** \brief . */
  typedef typename traits_type::pos_type 		pos_type;
  /** \brief . */
  typedef typename traits_type::off_type 		off_type;

  /** \brief . */
  basic_FancyOStream_buf(
    const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr
    ,const int                                                     startingTab
    ,const bool                                                    showLinePrefix
    ,const int                                                     maxLenLinePrefix
    ,const bool                                                    showTabCount
    ,const bool                                                    showProcRank
    );

  /** \brief . */
  void initialize(
    const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr
    ,const int                                                     startingTab
    ,const bool                                                    showLinePrefix
    ,const int                                                     maxLenLinePrefix
    ,const bool                                                    showTabCount
    ,const bool                                                    showProcRank
    );

  /** \brief . */
  RefCountPtr<std::basic_ostream<char_type,traits_type> > getOStream();

  /** \brief . */
  void setTabIndentStr(const std::basic_string<char_type,traits_type> &tabIndentStr);

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTabIndentStr() const;

  /** \brief .*/
  void setShowLinePrefix(const bool showLinePrefix);

  /** \brief .*/
  bool getShowLinePrefix() const;

  /** \brief .*/
  void setMaxLenLinePrefix(const int maxLenLinePrefix);

  /** \brief .*/
  int getMaxLenLinePrefix() const;

  /** \brief . */
  void setShowTabCount(const bool showTabCount);

  /** \brief . */
  bool getShowTabCount() const;

  /** \brief . */
  void setShowProcRank(const bool showProcRank);

  /** \brief . */
  bool getShowProcRank() const;

  /** \brief .*/
  void setProcRankAndSize( const int procRank, const int numProcs );

  /** \brief .*/
  int getProcRank() const;

  /** \brief .*/
  int getNumProcs() const;

  /** \brief . */
  void setOutputToRootOnly( const int rootRank );

  /** \brief . */
  int getOutputToRootOnly() const;

  /** \brief . */
  void pushTab(const int tabs);

  /** \brief . */
  void popTab();

  /** \brief . */
  void pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix);

  /** \brief . */
  void popLinePrefix();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTopLinePrefix() const;

  /** \brief . */
  void pushDisableTabbing();

  /** \brief . */
  void popDisableTabbing();
  
protected:
  
  /** \name Protected overridden functions from std::basic_streambuf */
  //@{

  /** \brief . */
  std::streamsize xsputn(const char_type* s, std::streamsize n);

  /** \brief . */
  int_type overflow(int_type c);

#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS

  void imbue(const locale& l) 
    {
      std::cerr << "\ncalled imbue()\n";
      std::basic_streambuf<CharT,Traits>::imbue(l);
    }
  
  basic_streambuf<char_type,Traits>* 
  setbuf(char_type* s, streamsize n)
    {
      std::cerr << "\ncalled setbuf()\n";
      return std::basic_streambuf<CharT,Traits>::setbuf(s,n);
    }
      
  pos_type 
  seekoff(off_type a, ios_base::seekdir b,ios_base::openmode c)
    {
      std::cerr << "\ncalled seekoff()\n";
      return std::basic_streambuf<CharT,Traits>::seekoff(a,b,c);
    }

  pos_type 
  seekpos(pos_type a, ios_base::openmode b)
    {
      std::cerr << "\ncalled seekpos()\n";
      return std::basic_streambuf<CharT,Traits>::seekpos(a,b);
    }
  
  int 
  sync()
    {
      std::cerr << "\ncalled sync()\n";
      return std::basic_streambuf<CharT,Traits>::sync();
    }
  
  streamsize 
  showmanyc()
    {
      std::cerr << "\ncalled showmanyc()\n";
      return std::basic_streambuf<CharT,Traits>::showmanyc();
    }
  
  streamsize 
  xsgetn(char_type* s, streamsize n)
    {
      std::cerr << "\ncalled xsgetn()\n";
      return std::basic_streambuf<CharT,Traits>::xsgetn(s,n);
    }
  
  int_type 
  underflow()
    {
      std::cerr << "\ncalled underflow()\n";
      return std::basic_streambuf<CharT,Traits>::underflow();
    }

  int_type 
  uflow() 
    {
      std::cerr << "\ncalled uflow()\n";
      return std::basic_streambuf<CharT,Traits>::uflow();
    }

  int_type 
  pbackfail(int_type c = traits_type::eof())
    {
      std::cerr << "\ncalled pbackfail()\n";
      return std::basic_streambuf<CharT,Traits>::pbackfail(c);
    }

#endif // TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS

  //@}

private:

  // ////////////////////////
  // Private types

  typedef std::basic_string<char_type,traits_type>   string_t;
  typedef std::deque<int>                            tabIndentStack_t;
  typedef std::deque<string_t>                       linePrefixStack_t;

  // ////////////////////////
  // Private data members

  RefCountPtr<std::basic_ostream<char_type,traits_type> >  oStreamSet_;
  RefCountPtr<std::basic_ostream<char_type,traits_type> >  oStream_;
  std::basic_string<char_type,traits_type>                 tabIndentStr_;
  bool                                                     showLinePrefix_;
  int                                                      maxLenLinePrefix_;
  bool                                                     showTabCount_;
  bool                                                     showProcRank_;
  int                                                      rootRank_;
  int                                                      procRank_;
  int                                                      numProcs_;
  int                                                      rankPrintWidth_;
  
  int                      tabIndent_;
  tabIndentStack_t         tabIndentStack_;
  linePrefixStack_t        linePrefixStack_;
  int                      enableTabbingStack_;
  
  bool                     wroteNewline_;

  // ////////////////////////
  // Private member functions

  void writeChars( const char_type s[], std::streamsize n );

  void writeFrontMatter();

  // Not defined and not to be called
  basic_FancyOStream_buf();
  basic_FancyOStream_buf(const basic_FancyOStream_buf<CharT,Traits>&);
  basic_FancyOStream_buf<CharT,Traits> operator=(const basic_FancyOStream_buf<CharT,Traits>&);

};

/** \brief std::ostream subclass that performs the magic of indenting data
 * sent to an std::ostream object among other things.
 *
 * Use the typedef <tt>FancyOStream</tt> for support for the <tt>char</tt>
 * character type.
 *
 * Indentation of the stream is accomplished through creating
 * <tt>basic_OSTab</tt> objects.
 *
 * In addition to indenting output, this stream object can also print various
 * types of information at the beginning of each line.  The type of information
 * supported is:
 * <ul>
 * <li> Processor rank: Set using <tt>setShowProcRank()</tt>.
 * <li> Line prefix name: Set using <tt>showLinePrefix()</tt> and <tt>OSTab::OSTab()</tt>.
 * <li> Tab counts (useful for debugging): Set using <tt>setShowTabCount()</tt>.
 * </ul>
 *
 * See <tt>FancyOutputting_test.cpp</tt> for examples of how this class is
 * used and the output it generates.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_FancyOStream : public basic_ostream<CharT, Traits>
{
public:

  /** \name Public types */
  //@{

  /** \brief . */
  typedef CharT 					                  char_type;
  /** \brief . */
  typedef Traits 					                  traits_type;
  /** \brief . */
  typedef typename traits_type::int_type 		int_type;
  /** \brief . */
  typedef typename traits_type::pos_type 		pos_type;
  /** \brief . */
  typedef typename traits_type::off_type 		off_type;
  /** \brief . */

  /** \brief . */
  typedef basic_FancyOStream_buf<CharT,Traits> 	streambuf_t;
  /** \brief . */
  typedef basic_ostream<char_type, traits_type>	ostream_t;

  //@}

  /** \name Public client functions */
  //@{

  /** \brief . */
  explicit
  basic_FancyOStream(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "  "
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    ,const bool                                                    showProcRank      = false
    );

  /** \brief . */
  void initialize(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "  "
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    ,const bool                                                    showProcRank      = false
    );

  /** \brief. */
  RefCountPtr<std::basic_ostream<char_type,traits_type> > getOStream();

  /** \brief . */
  basic_FancyOStream& setTabIndentStr(const std::basic_string<char_type,traits_type> &tabIndentStr);

  /** \brief. */
  const std::basic_string<char_type,traits_type>& getTabIndentStr() const;

  /** \brief Set if processor rank, line prefixes, and tab counts are shown or not .*/
  basic_FancyOStream& setShowAllFrontMatter(const bool showAllFrontMatter);

  /** \brief .*/
  basic_FancyOStream& setShowLinePrefix(const bool showLinePrefix);

  /** \brief .*/
  basic_FancyOStream& setMaxLenLinePrefix(const int maxLenLinePrefix);

  /** \brief . */
  basic_FancyOStream& setShowTabCount(const bool showTabCount);

  /** \brief . */
  basic_FancyOStream& setShowProcRank(const bool showProcRank);

  /** \brief . */
  basic_FancyOStream& setProcRankAndSize( const int procRank, const int numProcs );

  /** \brief . */
  basic_FancyOStream& setOutputToRootOnly( const int rootRank );

  /** \brief . */
  void copyAllOutputOptions(const basic_FancyOStream<CharT,Traits> &oStream);

  //@}

  /** \name Functions designed to be used by basic_OSTab */
  //@{

  /** \brief. */
  void pushTab(const int tabs = 1);

  /** \brief. */
  void popTab();

  /** \brief . */
  void pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix);

  /** \brief . */
  void popLinePrefix();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTopLinePrefix() const;

  /** \brief . */
  void pushDisableTabbing();

  /** \brief . */
  void popDisableTabbing();

  //@}
  
private:

  streambuf_t	streambuf_;

  // Not defined and not to be called
  basic_FancyOStream();
  basic_FancyOStream(const basic_FancyOStream<CharT,Traits>&);
  basic_FancyOStream<CharT,Traits> operator=(const basic_FancyOStream<CharT,Traits>&);

};

/** \brief Get a FancyOStream from an std::ostream object.
 *
 * If the object already is a FancyOStream, then nothing has to be done.
 * Otherwise, a temp FancyOStream is created for this purpose.
 *
 * \relates basic_FancyOStream
 */
template <typename CharT, typename Traits>
RefCountPtr<basic_FancyOStream<CharT,Traits> >
getFancyOStream( const RefCountPtr<std::basic_ostream<CharT,Traits> > &out )
{
  RefCountPtr<basic_FancyOStream<CharT,Traits> >
    fancyOut = rcp_dynamic_cast<basic_FancyOStream<CharT,Traits> >(out);
  if(fancyOut.get())
    return fancyOut;
  return rcp(new basic_FancyOStream<CharT,Traits>(out));
}

/** \brief Tabbing class for helping to create formated, indented output for a
 * <tt>basic_FancyOStream</tt> object.
 *
 * Use the typedef <tt>OSStream</tt> for support for the <tt>char</tt>
 * character type.
 *
 * This class is used to create tab indents and set line prefix names for
 * output that is generated by a <tt>basic_FancyOStream</tt> object.
 *
 * \relates basic_FancyOStream
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_OSTab
{
public:

  /** \brief . */
  static const int DISABLE_TABBING = -99999; // This magic number should be just fine!
  /** \brief. */
  basic_OSTab(
    const RefCountPtr<basic_FancyOStream<CharT,Traits> >   &fancyOStream
    ,const int                                             tabs = 1
    ,const std::basic_string<CharT,Traits>                 linePrefix = ""
    )
    :fancyOStream_(fancyOStream)
    ,tabs_(tabs)
    ,linePrefix_(linePrefix)
    {
      updateState();
    }
  /** \brief. */
  basic_OSTab( const basic_OSTab &osTab )
    :fancyOStream_(osTab.fancyOStream_)
    ,tabs_(osTab.tabs_)
    {
      updateState();
    }
  /** \brief. */
  ~basic_OSTab()
    {
      if(fancyOStream_.get()) {
        if(tabs_==DISABLE_TABBING)
          fancyOStream_->popDisableTabbing();
        else
          fancyOStream_->popTab();
        if(linePrefix_.length()) fancyOStream_->popLinePrefix();
      }
    }
  /** \brief. */
  basic_OSTab<CharT,Traits>& operator=( const basic_OSTab &osTab )
    {
      fancyOStream_ = osTab.fancyOStream_;
      tabs_ = osTab.tabs_;
      updateState();
      return *this;
    }
  /** \brief. */
  basic_OSTab<CharT,Traits>& incrTab(const int tabs = 1)
    {
      tabs_ += tabs;
      if(fancyOStream_.get()) {
        fancyOStream_->popTab();
        fancyOStream_->pushTab(tabs_);
      }
      return *this;
    }
  /** \brief. */
  RefCountPtr<basic_FancyOStream<CharT,Traits> > getOStream() const
    {
      return fancyOStream_;
    }
  
private:
  
  RefCountPtr<basic_FancyOStream<CharT,Traits> >  fancyOStream_;
  int                                             tabs_;
  std::basic_string<CharT,Traits>                 linePrefix_;

  void updateState()
    {
      if(fancyOStream_.get()) {
        if(tabs_==DISABLE_TABBING)
          fancyOStream_->pushDisableTabbing();
        else
          fancyOStream_->pushTab(tabs_);
        if(linePrefix_.length()) fancyOStream_->pushLinePrefix(linePrefix_);
      }
    }
  
};

// ///////////////////////////////
// Typedefs

/** \brief .
 * \ingroup teuchos_outputting_grp
 */
typedef basic_FancyOStream<char> FancyOStream;

/** \brief .
 * \ingroup teuchos_outputting_grp
 */
typedef basic_OSTab<char> OSTab;

/** \brief .
 * \ingroup teuchos_outputting_grp
 */
#define TEUCHOS_OSTAB ::Teuchos::OSTab __localThisTab = this->getOSTab()

// ////////////////////////////////
// Defintions

//
// basic_FancyOStream_buf
//

template<typename CharT, typename Traits>
basic_FancyOStream_buf<CharT,Traits>::basic_FancyOStream_buf(
  const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  this->initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::initialize(
  const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  oStreamSet_ = oStream;
  oStream_ = oStream;
  tabIndentStr_ = tabIndentStr;
  showLinePrefix_ = showLinePrefix;
  maxLenLinePrefix_ = maxLenLinePrefix;
  showTabCount_ = showTabCount;
  showProcRank_ = showProcRank;
  rootRank_ = -1;
  procRank_ = GlobalMPISession::getRank();
  numProcs_ = GlobalMPISession::getNProc();
  rankPrintWidth_ = int(std::log10(float(numProcs_)))+1;
  tabIndent_ = startingTab;
  tabIndentStack_.resize(0);
  linePrefixStack_.resize(0);
  wroteNewline_ = true;
  enableTabbingStack_ = 0;
}

template<typename CharT, typename Traits>
RefCountPtr<std::basic_ostream<CharT,Traits> >
basic_FancyOStream_buf<CharT,Traits>::getOStream()
{
  return oStreamSet_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setTabIndentStr(const std::basic_string<char_type,traits_type> &tabIndentStr)
{
  tabIndentStr_ = tabIndentStr;
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream_buf<CharT,Traits>::getTabIndentStr() const
{
  return tabIndentStr_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  showLinePrefix_ = showLinePrefix;
}

template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowLinePrefix() const
{
  return showLinePrefix_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setMaxLenLinePrefix(const int maxLenLinePrefix)
{
  TEST_FOR_EXCEPT( !(maxLenLinePrefix>=5) );
  maxLenLinePrefix_ = maxLenLinePrefix;
}

template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getMaxLenLinePrefix() const
{
  return maxLenLinePrefix_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  showTabCount_ = showTabCount;
}

template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowTabCount() const
{
  return showTabCount_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowProcRank(const bool showProcRank)
{
  showProcRank_ = showProcRank;
}

template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowProcRank() const
{
  return showProcRank_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setProcRankAndSize( const int procRank, const int numProcs )
{
  procRank_ = procRank;
  numProcs_ = numProcs;
}

template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getProcRank() const
{
  return procRank_;
}

template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getNumProcs() const
{
  return numProcs_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setOutputToRootOnly( const int rootRank  )
{
  rootRank_ = rootRank;
  if(rootRank >= 0) {
    if(rootRank == procRank_)
      oStream_ = oStreamSet_;
    else
      oStream_ = Teuchos::rcp(new oblackholestream());
  }
  else {
    oStream_ = oStreamSet_;
  }
}

template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getOutputToRootOnly() const
{
  return rootRank_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushTab(const int tabs)
{
  tabIndentStack_.push_back(tabs);
  tabIndent_ += tabs;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popTab()
{
  tabIndent_ -= tabIndentStack_.back();
  tabIndentStack_.pop_back();
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix)
{
  TEST_FOR_EXCEPT(static_cast<int>(linePrefix.length()) > maxLenLinePrefix_);
  linePrefixStack_.push_back(linePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popLinePrefix()
{
  linePrefixStack_.pop_back();
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream_buf<CharT,Traits>::getTopLinePrefix() const
{
  return linePrefixStack_.back();
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushDisableTabbing()
{
  ++enableTabbingStack_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popDisableTabbing()
{
  --enableTabbingStack_;
}

// protected

template<typename CharT, typename Traits>
std::streamsize basic_FancyOStream_buf<CharT,Traits>::xsputn(const char_type* s, std::streamsize n)
{
#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS
  std::cerr << "\ncalled xsputn()\n";
#endif
  writeChars(s,n);
  return n;
}

template<typename CharT, typename Traits>
typename basic_FancyOStream_buf<CharT,Traits>::int_type 
basic_FancyOStream_buf<CharT,Traits>::overflow(int_type c)
{
#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS
  std::cerr << "\ncalled overflow()\n";
#endif
  if(c != traits_type::eof()) {
    const char_type cc[] = { traits_type::to_char_type(c) };
    this->writeChars(cc,1);
  }
  return traits_type::not_eof(c);
  //return std::basic_streambuf<CharT,Traits>::overflow(c);
}

// private

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeChars( const char_type s[], std::streamsize n )
{
  if(n == 0) return;
  std::streamsize p = 0, first_p = 0;
  bool done_outputting = false;
  const char_type newline = '\n';
  while( !done_outputting ) {
    // Find the next newline
    for( p = first_p; p < n; ++p ) {
      if(s[p] == newline) {
        break;
      }
    }
    if(p == n) {
      // We did not find a newline at the end!
      --p;
      done_outputting = true;
    }
    else if( p == n-1 && s[p] == newline ) {
      // The last character in the string is a newline
      done_outputting = true;
    }
    // Write the beginning of the line if we need to
    if(wroteNewline_) {
      writeFrontMatter();
      wroteNewline_ = false;
    }
    // Write up to the newline or the end of the string
    oStream_->write(s+first_p,p-first_p+1);
    if(s[p] == newline) {
      wroteNewline_ = true;
      if(rootRank_ < 0)
        oStream_->flush();
    }
    // Update for next search
    if(!done_outputting)
      first_p = p+1; 
  }
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeFrontMatter()
{
  bool didOutput = false;
  if(showProcRank_) {
    *oStream_ << "p=" << std::right << std::setw(rankPrintWidth_) << procRank_;
    didOutput = true;
  }
  if(showLinePrefix_) {
    if(didOutput)
      *oStream_ << ", ";
    *oStream_ << std::left << std::setw(maxLenLinePrefix_);
    if(linePrefixStack_.size()) {
      *oStream_ << this->getTopLinePrefix();
    }
    else {
      *oStream_ << "";
    }
    didOutput = true;
  }
  if(showTabCount_) {
    if(didOutput)
      *oStream_ << ", ";
    *oStream_ << "tabs=" << std::right << std::setw(2) << tabIndent_;
    didOutput = true;
  }
  // ToDo: Add the Prefix name if asked
  // ToDo: Add the processor number if asked
  // ToDo: Add the number of indents if asked
  if(didOutput) {
    *oStream_ << " |" << tabIndentStr_;
  }
  if(enableTabbingStack_==0) {
    for( int i = 0; i < tabIndent_; ++i )
      *oStream_ << tabIndentStr_;
  }
}

//
// basic_FancyOStream
//

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>::basic_FancyOStream(
  const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
 )
  : ostream_t(NULL), streambuf_(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank)
{
  this->init(&streambuf_);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::initialize(
  const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  streambuf_.initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank);
  this->init(&streambuf_);
}

template<typename CharT, typename Traits>
RefCountPtr<std::basic_ostream<CharT,Traits> >
basic_FancyOStream<CharT,Traits>::getOStream()
{
  return streambuf_.getOStream();
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setTabIndentStr(const std::basic_string<char_type,traits_type> &tabIndentStr)
{
  streambuf_.setTabIndentStr(tabIndentStr);
  return *this;
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::getTabIndentStr() const
{
  return streambuf_.getTabIndentStr();
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowAllFrontMatter(const bool showAllFrontMatter)
{
  streambuf_.setShowLinePrefix(showAllFrontMatter);
  streambuf_.setShowTabCount(showAllFrontMatter);
  streambuf_.setShowProcRank(showAllFrontMatter);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  streambuf_.setShowLinePrefix(showLinePrefix);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
 basic_FancyOStream<CharT,Traits>::setMaxLenLinePrefix(const int maxLenLinePrefix)
{
  streambuf_.setMaxLenLinePrefix(maxLenLinePrefix);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  streambuf_.setShowTabCount(showTabCount);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowProcRank(const bool showProcRank)
{
  streambuf_.setShowProcRank(showProcRank);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setProcRankAndSize( const int procRank, const int numProcs )
{
  streambuf_.setProcRankAndSize(procRank,numProcs);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setOutputToRootOnly( const int rootRank )
{
  streambuf_.setOutputToRootOnly(rootRank);
  return *this;
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::copyAllOutputOptions( const basic_FancyOStream<CharT,Traits> &oStream )
{
  //streambuf_.setTabIndentStr(oStream.streambuf_.getTabIndentStr());
  streambuf_.setShowLinePrefix(oStream.streambuf_.getShowLinePrefix());
  streambuf_.setMaxLenLinePrefix(oStream.streambuf_.getMaxLenLinePrefix());
  streambuf_.setShowTabCount(oStream.streambuf_.getShowTabCount());
  streambuf_.setShowProcRank(oStream.streambuf_.getShowProcRank());
  streambuf_.setProcRankAndSize(oStream.streambuf_.getProcRank(),oStream.streambuf_.getNumProcs());
  streambuf_.setOutputToRootOnly(oStream.streambuf_.getOutputToRootOnly());
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushTab(const int tabs)
{
  streambuf_.pushTab(tabs);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popTab()
{
  streambuf_.popTab();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix)
{
  streambuf_.pushLinePrefix(linePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popLinePrefix()
{
  streambuf_.popLinePrefix();
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::getTopLinePrefix() const
{
  return streambuf_.getTopLinePrefix();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushDisableTabbing()
{
  streambuf_.pushDisableTabbing();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popDisableTabbing()
{
  return streambuf_.popDisableTabbing();
}

} // namespace Teuchos

#endif // TEUCHOS_FANCY_O_STREAM_HPP
