/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Timer.hpp>
#include <iostream>
#include <cstdlib>

namespace TestTeamVector {
#define VECTORLENGTH 16
typedef Kokkos::TeamVectorPolicy<VECTORLENGTH> Policy;

struct my_complex {
  double re,im;
  int dummy;
  KOKKOS_INLINE_FUNCTION
  my_complex() {
    re = 0.0;
    im = 0.0;
    dummy = 0;
  }
  KOKKOS_INLINE_FUNCTION
  my_complex(const my_complex& src) {
    re = src.re;
    im = src.im;
    dummy = src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const double& val) {
    re = val;
    im = 0.0;
    dummy = 0;
  }
  KOKKOS_INLINE_FUNCTION
  my_complex& operator += (const my_complex& src) {
    re += src.re;
    im += src.im;
    dummy += src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile my_complex& src) volatile {
    re += src.re;
    im += src.im;
    dummy += src.dummy;
  }
  KOKKOS_INLINE_FUNCTION
  my_complex& operator *= (const my_complex& src) {
    double re_tmp = re*src.re - im*src.im;
    double im_tmp = re * src.im + im * src.re;
    re = re_tmp;
    im = im_tmp;
    dummy *= src.dummy;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION
  void operator *= (const volatile my_complex& src) volatile {
    double re_tmp = re*src.re - im*src.im;
    double im_tmp = re * src.im + im * src.re;
    re = re_tmp;
    im = im_tmp;
    dummy *= src.dummy;
  }
  KOKKOS_INLINE_FUNCTION
  bool operator == (const my_complex& src) {
    return (re == src.re) && (im == src.im) && ( dummy == src.dummy );
  }
  KOKKOS_INLINE_FUNCTION
  bool operator != (const my_complex& src) {
      return (re != src.re) || (im != src.im) || ( dummy != src.dummy );
  }
  KOKKOS_INLINE_FUNCTION
  bool operator != (const double& val) {
    return (re != val) ||
           (im != 0) || (dummy != 0);
  }
  KOKKOS_INLINE_FUNCTION
  my_complex& operator= (const int& val) {
    re = val;
    im = 0.0;
    dummy = 0;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION
  my_complex& operator= (const double& val) {
    re = val;
    im = 0.0;
    dummy = 0;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION
  operator double() {
    return re;
  }
};

#if defined (KOKKOS_HAVE_CXX11)

template<typename Scalar, class ExecutionSpace>
struct functor_team_for {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_team_for(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  unsigned team_shmem_size(int team_size) const {return team_size*13*sizeof(Scalar)+8;}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {

    typedef typename ExecutionSpace::scratch_memory_space shmem_space ;
    typedef Kokkos::View<Scalar*,shmem_space,Kokkos::MemoryUnmanaged> shared_int;
    shared_int values = shared_int(team.team_shmem(),team.team_size()*13);

    team.vector_single([&] () {
      values(team.team_rank()) = 0;
    });

    team.team_par_for(131,[&] (int i) {
      team.vector_single([&] () {
        values(team.team_rank()) += i - team.league_rank() + team.league_size() + team.team_size();
      });
    });

    team.team_barrier();
    if(team.team_rank()==0) {
      team.vector_single([&] () {
        Scalar test = 0;
        Scalar value = 0;
        for(int i = 0; i < 131; i++) {
          test += i - team.league_rank() + team.league_size() + team.team_size();
        }
        for(int i=0; i<team.team_size(); i++) {
          value += values(i);
        }
        if(test!=value) {
          printf("FAILED vector_par_for %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) test,(double) value);
          flag()=1;
        }
      });
    }
  }
};

template<typename Scalar, class ExecutionSpace>
struct functor_vec_single {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_vec_single(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {
    Scalar value = 0;

    team.vector_par_for(13,[&] (int i) {
      value = i;
    });

    team.vector_single([&] () {
    },value);

    Scalar value2 = 0;
    team.vector_par_reduce(13, [&] (int i, Scalar& val) {
      val += value;
    },value2);

    if(value2!=(value*13)) {
      printf("FAILED vector_single broadcast %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) value2,(double) value);
      flag()=1;
    }
  }
};

template<typename Scalar, class ExecutionSpace>
struct functor_vec_for {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_vec_for(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  unsigned team_shmem_size(int team_size) const {return team_size*13*sizeof(Scalar)+8;}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {

    typedef typename ExecutionSpace::scratch_memory_space shmem_space ;
    typedef Kokkos::View<Scalar*,shmem_space,Kokkos::MemoryUnmanaged> shared_int;
    shared_int values = shared_int(team.team_shmem(),team.team_size()*13);

    team.vector_par_for(13,[&] (int i) {
      values(13*team.team_rank() + i) = i - team.team_rank() - team.league_rank() + team.league_size() + team.team_size();
    });

    team.vector_single([&] () {
      Scalar test = 0;
      Scalar value = 0;
      for(int i = 0; i < 13; i++) {
        test += i - team.team_rank() - team.league_rank() + team.league_size() + team.team_size();
        value += values(13*team.team_rank() + i);
      }
      if(test!=value) {
        printf("FAILED vector_par_for %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) test,(double) value);
        flag()=1;
      }
    });
  }
};

template<typename Scalar, class ExecutionSpace>
struct functor_vec_red {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_vec_red(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {
    Scalar value = 0;

    team.vector_par_reduce(13,[&] (int i, Scalar& val) {
      val += i;
    }, value);

    team.vector_single([&] () {
      Scalar test = 0;
      for(int i = 0; i < 13; i++) {
        test+=i;
      }
      if(test!=value) {
        printf("FAILED vector_par_reduce %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) test,(double) value);
        flag()=1;
      }
    });
  }
};

template<typename Scalar, class ExecutionSpace>
struct functor_vec_red_join {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_vec_red_join(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {
    Scalar value = 1;

    team.vector_par_reduce(13,[&] (int i, Scalar& val) {
                                val = i+1; }
      , value ,
      [&] (Scalar& val, const Scalar& src) {val*=src;}
    );

    team.vector_single([&] () {
      Scalar test = 1;
      for(int i = 0; i < 13; i++) {
        test*=i+1;
      }
      if(test!=value) {
        printf("FAILED vector_par_reduce_join %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) test,(double) value);
        flag()=1;
      }
    });
  }
};

template<typename Scalar, class ExecutionSpace>
struct functor_vec_scan {
  typedef Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace> policy_type;
  typedef ExecutionSpace execution_space;

  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag;
  functor_vec_scan(Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> flag_):flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (typename policy_type::member_type team) const {
    Scalar reduce_value;
    team.vector_par_scan(13,[&] (int i, Scalar& val, bool final) {
      val += i;
      if(final) {
        Scalar test = 0;
        for(int k = 0; k <= i; k++) {
          test+=k;
        }
        if(test!=val) {
          printf("FAILED vector_par_scan %i %i %lf %lf\n",team.league_rank(),team.team_rank(),(double) test,(double) val);
          flag()=1;
        }
      }
    },reduce_value);
  }
};

#endif

template<typename Scalar,class ExecutionSpace>
bool test_scalar(int nteams, int team_size, int test) {
  Kokkos::View<int,Kokkos::LayoutLeft,ExecutionSpace> d_flag("flag");
  Kokkos::View<int,Kokkos::LayoutLeft,Kokkos::HostSpace> h_flag("flag");
  h_flag() = 0;
  Kokkos::deep_copy(d_flag,h_flag);
  #ifdef KOKKOS_HAVE_CXX11
  if(test==0)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_vec_red<Scalar, ExecutionSpace>(d_flag));
  if(test==1)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_vec_red_join<Scalar, ExecutionSpace>(d_flag));
  if(test==2)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_vec_scan<Scalar, ExecutionSpace>(d_flag));
  if(test==3)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_vec_for<Scalar, ExecutionSpace>(d_flag));
  if(test==4)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_vec_single<Scalar, ExecutionSpace>(d_flag));
  if(test==5)
  Kokkos::parallel_for( Kokkos::TeamVectorPolicy<VECTORLENGTH,ExecutionSpace>(nteams,team_size),
      functor_team_for<Scalar, ExecutionSpace>(d_flag));
  #endif
  Kokkos::deep_copy(h_flag,d_flag);

  return (h_flag() == 0);
}

template<class ExecutionSpace>
bool Test(int test) {
  bool passed = true;
  passed = passed && test_scalar<int, ExecutionSpace>(317,33,test);
  passed = passed && test_scalar<long long int, ExecutionSpace>(317,33,test);
  passed = passed && test_scalar<float, ExecutionSpace>(317,33,test);
  passed = passed && test_scalar<double, ExecutionSpace>(317,33,test);
  passed = passed && test_scalar<my_complex, ExecutionSpace>(317,33,test);
  return passed;
}

}

