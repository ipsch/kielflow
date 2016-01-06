#ifndef O_MATH_HPP_
#define O_MATH_HPP_



template <class T> const T& max (const T& a, const T& b)
{
  return (b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

template <class T> const T& min (const T& a, const T& b)
{
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}





#endif
