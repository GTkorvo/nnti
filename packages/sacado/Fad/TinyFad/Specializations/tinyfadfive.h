#ifndef _tinyfadfive_h_
#define _tinyfadfive_h_

//*********************************************************
// This file is generated by generate.cc.
// Use this file for any modification
//*********************************************************

template <class T> class TinyFad<5,T> {
public:
  typedef T value_type;
  typedef T grad_type;
protected:

  int n;
  T val_;

  T dx0_;
  T dx1_;
  T dx2_;
  T dx3_;
  T dx4_;


public:
  void diff(const size_t ith, const size_t sz){
     n = ith+1;
     dx0_ = T(0.);
     dx1_ = T(0.);
     dx2_ = T(0.);
     dx3_ = T(0.);
     dx4_ = T(0.);

     switch(ith){
     case 0 : dx0_ = T(1.);break;
     case 1 : dx1_ = T(1.);break;
     case 2 : dx2_ = T(1.);break;
     case 3 : dx3_ = T(1.);break;
     case 4 : dx4_ = T(1.);break;
     default : std::cout << "ith = " << ith << "  out of definition set" << std::endl;exit(1);
     }
  }

  TinyFad(const T& ind, const int ini) : n(ini+1), val_(ind) {
     dx0_ = T(0.);
     dx1_ = T(0.);
     dx2_ = T(0.);
     dx3_ = T(0.);
     dx4_ = T(0.);

     switch(ini){
     case 0 : dx0_ = T(1.);break;
     case 1 : dx1_ = T(1.);break;
     case 2 : dx2_ = T(1.);break;
     case 3 : dx3_ = T(1.);break;
     case 4 : dx4_ = T(1.);break;
     default : std::cout << "ini = " << ini << "  out of definition set" << std::endl;exit(1);
     }
  }
  TinyFad() : n(0), val_(0.) {
     dx0_ = T(0.);
     dx1_ = T(0.);
     dx2_ = T(0.);
     dx3_ = T(0.);
     dx4_ = T(0.);
  }
  TinyFad(const No_Initialization &): n(0) {}
  TinyFad(const T& in) : n(0), val_(in) {
     dx0_ = T(0.);
     dx1_ = T(0.);
     dx2_ = T(0.);
     dx3_ = T(0.);
     dx4_ = T(0.);
  }
  TinyFad(const TinyFad<5,T> & in) : n(0), val_(in.val_) {
     dx0_ = in.dx0_;
     dx1_ = in.dx1_;
     dx2_ = in.dx2_;
     dx3_ = in.dx3_;
     dx4_ = in.dx4_;
  }

  ~TinyFad() {}

  int    N()       const {return n-1;}

  const T& val()     const { return val_;}
  T& val()                 { return val_;}

  const T& d0() const { return dx0_;}
  T& d0() { return dx0_;}

  const T& d1() const { return dx1_;}
  T& d1() { return dx1_;}

  const T& d2() const { return dx2_;}
  T& d2() { return dx2_;}

  const T& d3() const { return dx3_;}
  T& d3() { return dx3_;}

  const T& d4() const { return dx4_;}
  T& d4() { return dx4_;}

  T& dx(int i){
     switch(i){
     case 0 : return dx0_;
     case 1 : return dx1_;
     case 2 : return dx2_;
     case 3 : return dx3_;
     case 4 : return dx4_;
     default : std::cout << "i out of bounds" << std::endl;exit(1);
     }
  }
  const T& dx(int i) const {
     switch(i){
     case 0 : return dx0_;
     case 1 : return dx1_;
     case 2 : return dx2_;
     case 3 : return dx3_;
     case 4 : return dx4_;
     default : std::cout << "i out of bounds" << std::endl;exit(1);
     }
  }
  T& d(int i){
     switch(i){
     case 0 : return dx0_;
     case 1 : return dx1_;
     case 2 : return dx2_;
     case 3 : return dx3_;
     case 4 : return dx4_;
     default : std::cout << "i out of bounds" << std::endl;exit(1);
     }
  }
  const T& d(int i) const {
     switch(i){
     case 0 : return dx0_;
     case 1 : return dx1_;
     case 2 : return dx2_;
     case 3 : return dx3_;
     case 4 : return dx4_;
     default : std::cout << "i out of bounds" << std::endl;exit(1);
     }
  }

  TinyFad<5,T> & operator = (const TinyFad<5,T> & in){
     val_ = in.val_;

     dx0_ = in.dx0_;
     dx1_ = in.dx1_;
     dx2_ = in.dx2_;
     dx3_ = in.dx3_;
     dx4_ = in.dx4_;

     return *this;
  }

  TinyFad<5,T> & operator = (const T & in){
     val_ = in;

     dx0_ = T(0.);
     dx1_ = T(0.);
     dx2_ = T(0.);
     dx3_ = T(0.);
     dx4_ = T(0.);

     return *this;
  }

  TinyFad<5,T> & operator += (const TinyFad<5,T> & in){
     dx0_ += in.dx0_;
     dx1_ += in.dx1_;
     dx2_ += in.dx2_;
     dx3_ += in.dx3_;
     dx4_ += in.dx4_;
     val_ += in.val_;


    return *this;
  }
  TinyFad<5,T> & operator -= (const TinyFad<5,T> & in){
     dx0_ -= in.dx0_;
     dx1_ -= in.dx1_;
     dx2_ -= in.dx2_;
     dx3_ -= in.dx3_;
     dx4_ -= in.dx4_;
     val_ -= in.val_;

     return *this;
  }
  TinyFad<5,T> & operator *= (const TinyFad<5,T> & in){
     dx0_ = dx0_ * in.val_ + val_ * in.dx0_;
     dx1_ = dx1_ * in.val_ + val_ * in.dx1_;
     dx2_ = dx2_ * in.val_ + val_ * in.dx2_;
     dx3_ = dx3_ * in.val_ + val_ * in.dx3_;
     dx4_ = dx4_ * in.val_ + val_ * in.dx4_;
     val_ *= in.val_;

     return *this;
  }
  TinyFad<5,T> & operator /= (const TinyFad<5,T> & in){
     if (in.val_ == 0.) error("TinyFad & TinyFad::operator /= (const TinyFad & in), dividing by 0");
     dx0_ = ( dx0_ * in.val_ - val_ * in.dx0_ ) / in.val_ / in.val_ ;
     dx1_ = ( dx1_ * in.val_ - val_ * in.dx1_ ) / in.val_ / in.val_ ;
     dx2_ = ( dx2_ * in.val_ - val_ * in.dx2_ ) / in.val_ / in.val_ ;
     dx3_ = ( dx3_ * in.val_ - val_ * in.dx3_ ) / in.val_ / in.val_ ;
     dx4_ = ( dx4_ * in.val_ - val_ * in.dx4_ ) / in.val_ / in.val_ ;
     val_ /= in.val_;

     return *this;
  }

  TinyFad<5,T> & operator += (const T & in){     val_ += in;

     return *this;
  }
  TinyFad<5,T> & operator -= (const T & in){     val_ -= in;

     return *this;
  }
  TinyFad<5,T> & operator *= (const T & in){
     val_ *= in;

     dx0_ *= in;
     dx1_ *= in;
     dx2_ *= in;
     dx3_ *= in;
     dx4_ *= in;

     return *this;
  }
  TinyFad<5,T> & operator /= (const T & in){
     if ( in == T(0.) ) error("TinyFad & TinyFad::operator /= (const T & in), dividing by 0");
     val_ /= in;

     dx0_ /= in;
     dx1_ /= in;
     dx2_ /= in;
     dx3_ /= in;
     dx4_ /= in;

     return *this;
  }

  TinyFad<5,T> operator++(int){
     TinyFad<5,T> tmp(*this);
     tmp.val_++;
     return tmp;
  };
  TinyFad<5,T> operator--(int){
     TinyFad<5,T> tmp(*this);
     tmp.val_--;
     return tmp;
  };
  TinyFad<5,T> & operator++(){     ++val_;
     return *this;
  }
  TinyFad<5,T> & operator--(){     --val_;
     return *this;
  }
};


template <class T> inline TinyFad<5,T> operator + (const TinyFad<5,T>& in)
{
  return TinyFad<5,T>(in);
}

template <class T> inline TinyFad<5,T> operator - (const TinyFad<5,T>& in)
{
  TinyFad<5,T> tmp;
  tmp -= in;
  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator +(const TinyFad<5,L>& un, const TinyFad<5,R>& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0() + deux.d0();

  tmp.d1() = un.d1() + deux.d1();

  tmp.d2() = un.d2() + deux.d2();

  tmp.d3() = un.d3() + deux.d3();

  tmp.d4() = un.d4() + deux.d4();

  tmp.val() = un.val() + deux.val();

  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote>
operator +(const TinyFad<5,L>& un, const R& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0();

  tmp.d1() = un.d1();

  tmp.d2() = un.d2();

  tmp.d3() = un.d3();

  tmp.d4() = un.d4();

  tmp.val() = un.val() + deux;

  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator +(const L& un, const TinyFad<5,R>& deux) {
  return operator +(deux,un);
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator *(const TinyFad<5,L>& un, const TinyFad<5,R>& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0()*deux.val() + un.val() * deux.d0();

  tmp.d1() = un.d1()*deux.val() + un.val() * deux.d1();

  tmp.d2() = un.d2()*deux.val() + un.val() * deux.d2();

  tmp.d3() = un.d3()*deux.val() + un.val() * deux.d3();

  tmp.d4() = un.d4()*deux.val() + un.val() * deux.d4();

  tmp.val() = un.val() * deux.val();

  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator *(const TinyFad<5,L>& un, const R& deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0()*deux;

  tmp.d1() = un.d1()*deux;

  tmp.d2() = un.d2()*deux;

  tmp.d3() = un.d3()*deux;

  tmp.d4() = un.d4()*deux;

  tmp.val() = un.val() * deux;

  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator *(const L& un, const TinyFad<5,R>& deux) {

  return operator *(deux,un);
}


template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote >
operator -(const TinyFad<5,L> & un, const TinyFad<5,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0() - deux.d0();

  tmp.d1() = un.d1() - deux.d1();

  tmp.d2() = un.d2() - deux.d2();

  tmp.d3() = un.d3() - deux.d3();

  tmp.d4() = un.d4() - deux.d4();

  tmp.val() = un.val() - deux.val();

  return tmp;
}

template <class L, class R> inline
TinyFad<5,typename NumericalTraits<L,R>::promote>
operator -(const L & un, const TinyFad<5,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() -= deux.d0();

  tmp.d1() -= deux.d1();

  tmp.d2() -= deux.d2();

  tmp.d3() -= deux.d3();

  tmp.d4() -= deux.d4();

  tmp.val() = un - deux.val();

  return tmp;
}

template <class L, class R> inline
TinyFad<5, typename NumericalTraits<L,R>::promote >
operator -(const TinyFad<5,L> & un, const R & deux) {
  typedef typename NumericalTraits<L,R>::promote value_type;

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0();

  tmp.d1() = un.d1();

  tmp.d2() = un.d2();

  tmp.d3() = un.d3();

  tmp.d4() = un.d4();

  tmp.val() = un.val() - deux;

  return tmp;
}

template <class L, class R> inline
TinyFad<5, typename NumericalTraits<L,R>::promote >
operator /(const TinyFad<5,L> & un, const TinyFad<5,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux.val() == 0.) error("TinyFad & TinyFad::operator /(const TinyFad<5,L> & un, const TinyFad<5,R> & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );
  value_type dval = deux.val();

  tmp.d0() = ( un.d0()* deux.val() - un.val() * deux.d0() ) / dval / dval ;

  tmp.d1() = ( un.d1()* deux.val() - un.val() * deux.d1() ) / dval / dval ;

  tmp.d2() = ( un.d2()* deux.val() - un.val() * deux.d2() ) / dval / dval ;

  tmp.d3() = ( un.d3()* deux.val() - un.val() * deux.d3() ) / dval / dval ;

  tmp.d4() = ( un.d4()* deux.val() - un.val() * deux.d4() ) / dval / dval ;

  tmp.val() = un.val() / dval;

  return tmp;
}

template <class L, class R> inline
TinyFad<5, typename NumericalTraits<L,R>::promote >
operator /(const L & un, const TinyFad<5,R> & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux.val() == 0.) error("TinyFad & TinyFad::operator /(const L & un, const TinyFad<5,R> & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );
  value_type dval = deux.val();

  tmp.d0() = - un * deux.d0()  / dval / dval ;

  tmp.d1() = - un * deux.d1()  / dval / dval ;

  tmp.d2() = - un * deux.d2()  / dval / dval ;

  tmp.d3() = - un * deux.d3()  / dval / dval ;

  tmp.d4() = - un * deux.d4()  / dval / dval ;

  tmp.val() = un / dval;

  return tmp;
}

template <class  L, class R> inline
TinyFad<5, typename NumericalTraits<L,R>::promote >
operator /(const TinyFad<5,L> & un, const R & deux) {

  typedef typename NumericalTraits<L,R>::promote value_type;

  if (deux == 0.) error("TinyFad & TinyFad::operator /(const TinyFad<5,L> & un, const R & deux), dividing by 0");

  No_Initialization nothing;
  TinyFad<5,value_type> tmp( nothing );

  tmp.d0() = un.d0()  / deux;

  tmp.d1() = un.d1()  / deux;

  tmp.d2() = un.d2()  / deux;

  tmp.d3() = un.d3()  / deux;

  tmp.d4() = un.d4()  / deux;

  tmp.val() = un.val() / deux;

   return tmp;
}

template <class T> TinyFad<5,T> exp (const TinyFad<5,T>& in)
{
  TinyFad<5,T> tmp(exp(in.val()));

  tmp.d0() = in.d0()*exp(in.val());
  tmp.d1() = in.d1()*exp(in.val());
  tmp.d2() = in.d2()*exp(in.val());
  tmp.d3() = in.d3()*exp(in.val());
  tmp.d4() = in.d4()*exp(in.val());

  return tmp;
}

template <class T> TinyFad<5,T> log (const TinyFad<5,T>& in)
{
  if ( in.val() <= 0) error("TinyFad log (const TinyFad& in) : zero or negative value");
  TinyFad<5,T> tmp(log(in.val()));

  tmp.d0() = in.d0() / in.val();
  tmp.d1() = in.d1() / in.val();
  tmp.d2() = in.d2() / in.val();
  tmp.d3() = in.d3() / in.val();
  tmp.d4() = in.d4() / in.val();

  return tmp;
}

template <class T> TinyFad<5,T> sqrt (const TinyFad<5,T>& in)
{
  if ( in.val() < 0. ) error("TinyFad<5,T> sqrt (const TinyFad& in) : negative value");
  TinyFad<5,T> tmp(sqrt(in.val()));

  bool test=true;
  if ( in.val() == T(0.) ){
    if ( in.d0() != T(0.) ) test = false;

    if ( in.d1() != T(0.) ) test = false;

    if ( in.d2() != T(0.) ) test = false;

    if ( in.d3() != T(0.) ) test = false;

    if ( in.d4() != T(0.) ) test = false;

    if ( !test )
      error("TinyFad<5,T> sqrt (const TinyFad& in) : null value");
  }
  else {
    tmp.d0() = in.d0() / sqrt(in.val()) / 2.;
    tmp.d1() = in.d1() / sqrt(in.val()) / 2.;
    tmp.d2() = in.d2() / sqrt(in.val()) / 2.;
    tmp.d3() = in.d3() / sqrt(in.val()) / 2.;
    tmp.d4() = in.d4() / sqrt(in.val()) / 2.;
  }

  return tmp;
}

template <class T> TinyFad<5,T> sin (const TinyFad<5,T>& in)
{
  TinyFad<5,T> tmp( sin(in.val()) );

  tmp.d0() = in.d0() * cos( in.val() );
  tmp.d1() = in.d1() * cos( in.val() );
  tmp.d2() = in.d2() * cos( in.val() );
  tmp.d3() = in.d3() * cos( in.val() );
  tmp.d4() = in.d4() * cos( in.val() );

  return tmp;
}

template <class T> TinyFad<5,T> cos (const TinyFad<5,T>& in)
{
  TinyFad<5,T> tmp(cos(in.val()));

  tmp.d0() = - in.d0() * sin( in.val() );
  tmp.d1() = - in.d1() * sin( in.val() );
  tmp.d2() = - in.d2() * sin( in.val() );
  tmp.d3() = - in.d3() * sin( in.val() );
  tmp.d4() = - in.d4() * sin( in.val() );

  return tmp;
}

template <class T> TinyFad<5,T> tan (const TinyFad<5,T>& in)
{
  if ( in.val() == 0) error("TinyFad tan (const TinyFad& in) undiefined in 0.");
  TinyFad<5,T> tmp(tan(in.val()));

  T cosinus = cos(in.val());
  tmp.d0() = in.d0() / cosinus / cosinus;
  tmp.d1() = in.d1() / cosinus / cosinus;
  tmp.d2() = in.d2() / cosinus / cosinus;
  tmp.d3() = in.d3() / cosinus / cosinus;
  tmp.d4() = in.d4() / cosinus / cosinus;

  return tmp;
}

template <class T> TinyFad<5,T> pow (const TinyFad<5,T>& in, double e)
{
  TinyFad<5,T> tmp(pow(in.val(), e));

  tmp.d0() = e*in.d0()*pow(in.val(), e-1);
  tmp.d1() = e*in.d1()*pow(in.val(), e-1);
  tmp.d2() = e*in.d2()*pow(in.val(), e-1);
  tmp.d3() = e*in.d3()*pow(in.val(), e-1);
  tmp.d4() = e*in.d4()*pow(in.val(), e-1);

  return tmp;
}

template <class T> TinyFad<5,T> pow (const TinyFad<5,T>& un, const TinyFad<5,T>& deux)
{
  if (un.val() == 0) error("TinyFad pow (const TinyFad& un, const TinyFad& deux), un = 0. ");
  TinyFad<5,T> tmp(pow(un.val(), deux.val()));

  tmp.d0() = deux.d0() * log(un.val()) * pow(un.val(), deux.val())
	    + deux.val() * un.d0() * pow(un.val(), deux.val()-1);
  tmp.d1() = deux.d1() * log(un.val()) * pow(un.val(), deux.val())
	    + deux.val() * un.d1() * pow(un.val(), deux.val()-1);
  tmp.d2() = deux.d2() * log(un.val()) * pow(un.val(), deux.val())
	    + deux.val() * un.d2() * pow(un.val(), deux.val()-1);
  tmp.d3() = deux.d3() * log(un.val()) * pow(un.val(), deux.val())
	    + deux.val() * un.d3() * pow(un.val(), deux.val()-1);
  tmp.d4() = deux.d4() * log(un.val()) * pow(un.val(), deux.val())
	    + deux.val() * un.d4() * pow(un.val(), deux.val()-1);

  return tmp;
}

template <class T> TinyFad<5,T> pow (const TinyFad<5,T>& in, const int e)
{
  TinyFad<5,T> tmp( std::pow((double)in.val(), (double)e) );

  tmp.d0() = e*in.d0()*std::pow((double)in.val(), (double)e-1);
  tmp.d1() = e*in.d1()*std::pow((double)in.val(), (double)e-1);
  tmp.d2() = e*in.d2()*std::pow((double)in.val(), (double)e-1);
  tmp.d3() = e*in.d3()*std::pow((double)in.val(), (double)e-1);
  tmp.d4() = e*in.d4()*std::pow((double)in.val(), (double)e-1);

  return tmp;
}

template <class T> TinyFad<5,T> abs (const TinyFad<5,T>& in)
{
  int sign = in.val() > 0? 1:0;

  if (sign) return in;
  else return (-in);
}

template <class T> std::ostream& operator << (std::ostream& os, const TinyFad<5,T>& a)
{
  os.setf(std::ios::fixed,std::ios::scientific);
  os.width(12);
  os << a.val() << "  [";

  os.width(12);
  os << a.d0();
  os.width(12);
  os << a.d1();
  os.width(12);
  os << a.d2();
  os.width(12);
  os << a.d3();
  os.width(12);
  os << a.d4();

  os << "]";

  return os;
}


#endif
