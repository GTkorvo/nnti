#include "cmm.h"
#include "geom.h"

/*

I need to modify this for Doxygen.



The following is an experiment at defining xml documentation for the
code.

@ <component name=cmath/base/geom>
@   <author>Mark P. Sears</author>
@   <description>
@   <p>  This component contains functions for manipulation of
@    ordinary three dimensional vectors. Such a vector is
@    normally declared as an array of doubles of length 3:
@     <code>
@       r8 v[3];
@     </code>
@   </p>
@   </description>
@   <source file="geom.h">The include file.</source>
@   <source file="geom.c">Source code for the functions.</source>
@ </component>

*/

/*
@ <function name=axpy>
@  <declaration><code>void axpy3(r8 a, r8 * x, r8 * y)</code></declaration>
@
@   <argument><code>r8 a</code> Scalar factor</argument>
@   <argument><code>r8 *x</code> Vector</argument>
@   <argument><code>r8 *y</code> Vector</argument>
@
@  <description>
@    This function provides the operation y = y + a*x, where x and y are
@    vectors and a is a scalar.
@  </description>
@
@ </function>
*/

void axpy3(r8 a, r8 * x, r8 * y)
{
  y[0] += a * x[0];
  y[1] += a * x[1];
  y[2] += a * x[2];
}

void drot3(r8 c, r8 * x, r8 s, r8 * y, r8 * z)
{
  z[0] = c * x[0] + s * y[0];
  z[1] = c * x[1] + s * y[1];
  z[2] = c * x[2] + s * y[2];
}

void zero3(r8 *x)
{
  x[0] = 0.;
  x[1] = 0.;
  x[2] = 0.;
}

void set3(r8 tx, r8 ty, r8 tz, r8 *x)
{
  x[0] = tx;
  x[1] = ty;
  x[2] = tz;
}

void copy3(r8 * x, r8 * y)
{
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
}

r8 dist3(r8 * x, r8 * y)
{
  r8 xn;
  r8 t;

  t = x[0] - y[0];
  xn = t * t;
  t = x[1] - y[1];
  xn += t * t;
  t = x[2] - y[2];
  xn += t * t;
  return sqrt(xn);
}

r8 dot3(r8 * x, r8 * y)
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

void cross3(r8 * x, r8 * y, r8 * z)
{
  z[0] = x[1] * y[2] - x[2] * y[1];
  z[1] = x[2] * y[0] - x[0] * y[2];
  z[2] = x[0] * y[1] - x[1] * y[0];
}

void add3(r8 * x, r8 * y, r8 * z)
{
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
}

void diff3(r8 * x, r8 * y, r8 * z)
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
}

r8 norm3(r8 * x)
{
  return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

r8 normalize3(r8 * x)
{
  r8 s = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  r8 s1 = 1./s;
  x[0] *= s1;
  x[1] *= s1;
  x[2] *= s1;
  return s;
}

void scale3(r8 s, r8 * x, r8 * y)
{
  y[0] = s * x[0];
  y[1] = s * x[1];
  y[2] = s * x[2];
}

// compute angle in degrees between two vectors

r8 angle3(r8 * a, r8 * b)
{
  r8 t;

  t = dot3(a, b) / (norm3(a) * norm3(b));
  return acos(t) * (180. / CMM_PI);

}

