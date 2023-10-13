/**
 * Copyright 2014-2016 Richard Pausch
 *
 * This file is part of Clara 2.
 *
 * Clara 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Clara 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Clara 2.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VECTOR_H
#define VECTOR_H
//#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>

/**
 * \brief  Vector class by Richard
 *
 *  This class provides operations like: +, -, *, %, /, =, +=, [], mag() \n
 *  usage: Vector<datatype, size> \n
 *  short for Vector<double, 3> --> R_vec
 */

template< typename T, unsigned N>
class Vector{

public:
  // constructor, destructor:
  /// basic constructor setting everything to zero
  inline Vector(){
    for (unsigned i=0; i<N; ++i)
    data[i] = T(0.0);
  };

  //! \brief special constructor for 3D
  /*! @param x first coordinate
      @param y second coordinate
      @param z third coordinate */
  inline Vector(double x,
                double y=0.,
                double z=0.){  
    assert( N==3 ); // better via template? one more dummy function
    data[0]=x;
    data[1]=y;
    data[2]=z;
  };

  /// copy constructor --> not necessary

  /// destructor
  inline ~Vector() {}  // not necessary


  // Getters and Setters:
  /// access data
  inline T & operator[](unsigned i){
    assert(i<N);
    return data[i];
  };

  /// access data
  const inline T & operator[](unsigned i) const{
    assert(i<N);
    return data[i];
  };



  // Calculations:

  /// addition
  inline Vector operator+ (const  Vector& v) const{
    Vector back;
    for (unsigned i=0; i<N; ++i)
    {
      back[i] = data[i] + v.data[i];
    }
    return back;
  };

  /// subtraction
  inline Vector operator- (const Vector& v) const{
    Vector back;
    for (unsigned i=0; i<N; ++i)
    {
      back[i] = data[i] - v.data[i];
    }
    return back;
  };

  /// magnitude
  inline T mag() const{
    T sqr_magnitude=0;
    for (unsigned i=0; i<N; ++i)
    {
      sqr_magnitude += data[i] * data[i];
    }
    return std::sqrt(sqr_magnitude);
  };

  /// scalar product
  inline T dot(const Vector& v) const{
    T result=0;
    for (unsigned i=0; i<N; ++i)
    {
      result += data[i] * v[i];
    }
    return result;
  };

  ///  multiplication with scalar
  inline Vector dot(const T a) const{
    Vector result;
    for (unsigned i=0; i<N; ++i)
    {
      result[i] = a* data[i];
    }
    return result;
  };

  /// cross-product-warning
  inline Vector cross(const Vector& v) const{
    std::cout << std::endl <<
    "---------------------------" << std::endl <<
    " cross product not defined for N != 3 " << std::endl <<
    "---------------------------" << std::endl;
    assert(false);
    Vector dummy;
    return dummy;
  };

  /// assign addition
  inline Vector & operator += (const Vector& v){
    for (unsigned i=0; i<N; ++i)
    {
      data[i] += v[i];
    }
    return *this;
  };

  inline Vector unit_vec() const{
    return *this / mag() ;
  };


  ////////////

  /// make real vector complex
  Vector<std::complex<T>, N> make_complex() const{
    Vector<std::complex<T>, N> result;
    for(unsigned i=0; i<N; ++i)
      result[i] = std::complex<T>(data[i]);

    return result;
  };

  /// make complex vector real (absolute value of complex number)
  Vector<double, N> abs() const{
    Vector<double, N> result;
    for(unsigned i = 0; i< N; ++i)
      result[i] = std::abs(data[i]);

    return result;
  };



private:
  /// stored data
  T data[N];

};


// Member function specialization:


/// 3D cross product
template<>
inline Vector<double, 3> Vector<double, 3>::cross(const Vector<double, 3>& v) const{
  Vector<double, 3> result;
  result[0] = data[1]*v[2]-v[1]*data[2];
  result[1] = data[2]*v[0]-v[2]*data[0];
  result[2] = data[0]*v[1]-v[0]*data[1];

  return result;
};



// global methods --> symbols for calculations and printing

template< typename T, unsigned N>  /// Vector * Vector --> scalar product
inline double operator * (const Vector<T, N>& a ,
                          const Vector<double, 3>& b )
{
  return a.dot(b);
};

template< typename T, unsigned N>  /// Vector * scalar
inline Vector<T,N> operator * (const Vector<T, N> & v,
                               const double a)
{
  return v.dot(a);
};

template< typename T, unsigned N>  /// scalar * Vector
inline Vector<T,N> operator * (const double a,
                               const Vector<T, N> & v)
{
  return v.dot(a);
};

template< typename T, unsigned N>  /// Vector / scalar
inline Vector<T,N> operator / (const Vector<T, N> & v,
                               const double a)
{
  return v.dot(1/a);
};

template< typename T, unsigned N>  /// cross product --> Vector % Vector
inline Vector<T,N> operator % (const Vector<T,N> & a,
                               const Vector<T,N> & b)  // cross-product
{
  return a.cross(b);
};


//! \brief output stream used on vector object
/*! @param os output stream
    @param v vector  */
template< typename T, unsigned N>  /// print Vector
inline std::ostream & operator << (std::ostream & os,
                                   const Vector<T,N> & v)
{
  os << "(" ;
  for (unsigned i=0; i<(N-1); ++i)
  {
    os << v[i] << " , ";
  }
  os << v[N-1] << ")";
  return os;
};


/*! \var typedef Vector<double, 3> R_vec
  \brief A type definition for a 3D double vector.

  Because in physics this is the most widely used vector, there is a special
  typedef for it.
*/

//#include "Vector.cpp"

typedef Vector<double, 3> R_vec;

#endif // VECTOR_H