#ifndef GEOM_H
#define GEOM_H 1

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

  void zero3(r8 *x);                                // x = (0, 0, 0)
  void set3(r8 tx, r8 ty, r8 tz, r8 *x);            // x = (tx, ty, tz)
  void copy3(r8 * x, r8 * y);                       // y = x
  void add3(r8 * x, r8 * y, r8 * z);                // z = x + y
  void diff3(r8 * x, r8 * y, r8 * z);               // z = x - y
  void axpy3(r8 a, r8 * x, r8 * y);                 // y = y + a*x
  void scale3(r8 s, r8 * x, r8 * y);                // y = s*x
  void drot3(r8 c, r8 * x, r8 s, r8 * y, r8 * z);   // ?
  r8 dist3(r8 * x, r8 * y);                         // return |x-y|
  r8 dot3(r8 * x, r8 * y);                          // return x.y
  void cross3(r8 * x, r8 * y, r8 * z);              // z = x cross y
  r8 norm3(r8 * x);                                 // return |x|
  r8 normalize3(r8 * x);                            // return |x| and normalize x.
  r8 angle3(r8 * a, r8 * b);                        // return angle between two vectors

#ifdef __cplusplus
}
#endif                          /* __cplusplus */

#endif
