#if defined ( LONGREAL )

extern "C" long double sqrtl	  (long double x);
extern "C" long double fabsl	  (long double x);
extern "C" long double logl 	  (long double x);
extern "C" long double expl 	  (long double x);
extern "C" long double acosl    (long double x);
extern "C" long double asinl    (long double x);
extern "C" long double atanl    (long double x);
extern "C" long double atan2l   (long double x, long double y);
extern "C" long double ceill    (long double x);
extern "C" long double cosl     (long double x);
extern "C" long double coshl    (long double x);
extern "C" long double floorl   (long double x);
extern "C" long double fmodl    (long double x, long double y);
extern "C" long double log10l   (long double x);
extern "C" long double powl     (long double x, long double y);
extern "C" long double sinl     (long double x);
extern "C" long double sinhl    (long double x);
extern "C" long double tanl     (long double x);
extern "C" long double tanhl    (long double x);

inline real sqrt( real x )
{
	return sqrtl( x );
}

inline real fabs( real x )
{
	return fabsl( x );
}

inline real log ( real x )
{
	return logl( x );
}

inline real exp ( real x )
{
	return expl( x );
}

inline real acos    (real x)
{
   return acosl( x );
}

inline real asin    (real x)
{
   return asinl( x );
}

inline real atan    (real x)
{
   return atanl( x );
}

inline real atan2   (real x, real y)
{
   return atan2l( x, y );
}

inline real ceil    (real x)
{
   return ceill( x );
}

inline real cos     (real x)
{
   return cosl( x );
}

inline real cosh    (real x)
{
   return coshl( x );
}

inline real floor   (real x)
{
   return floorl( x );
}

inline real fmod    (real x, real y)
{
   return fmodl( x, y );
}

inline real log10   (real x)
{
   return log10l( x );
}

inline real pow     (real x, real y)
{
   return powl( x,y );
}

inline real sin     (real x)
{
   return sinl( x );
}

inline real sinh    (real x)
{
   return sinhl( x );
}

inline real tan     (real x)
{
   return tanl( x );
}

inline real tanh    (real x)
{
   return tanhl( x );
}

#if defined ( CC_SOURCE )
inline ostream &operator<<( ostream &out, const real &x )
	{ out << double( x ); return out; }

inline istream &operator>>( istream &in, real &x )
	{ double d; in >> d; x = real( d ); return in; }
#endif

 
#endif
