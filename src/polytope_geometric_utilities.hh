#ifndef POLYTOPE_GEOMETRIC_UTILITIES_HH
#define POLYTOPE_GEOMETRIC_UTILITIES_HH
//------------------------------------------------------------------------------
// A semi-random collection of stuff related to geometric computations for use
// internally in polytope.
//------------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>

#include "Tessellation.hh"
#include "ReducedPLC.hh"

#ifdef HAVE_MPI
#include <mpi.h>
#include "polytope_parallel_utilities.hh"
#endif

using std::abs;
using std::min;
using std::max;

// --- These live in predicates.cc
// Compute the orientation of point c relative to points a and b
extern double orient2d(double* a, double* b, double* c);
// Sets h = b*e
extern int scale_expansion(int elen, double* e, double b, double* h);
// Sets h = e + b
extern int grow_expansion(int elen, double* e, double b, double* h);

//------------------------------------------------------------------------------
// It seems there is a missing specialization for abs(long unsigned int), so 
// fill it in.
// This is necessary for the collinear method below to compile.  It seems evil
// to insert something into namespace std:: like this, by the way.
//------------------------------------------------------------------------------
namespace std {
  // inline long unsigned int abs(long unsigned int x) { return x; }
  inline uint64_t          abs(uint64_t x)          { return x; }
}

namespace polytope {
namespace geometry {

//------------------------------------------------------------------------------
// The sgn function.
//------------------------------------------------------------------------------
template<typename T> 
int sgn(T val) {
  return (val >= T(0) ? 1 : -1);
}

//------------------------------------------------------------------------------
// A handy method for computing a hash of a position to a 64 bit quantized
// integer value.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType> struct Hasher;

// 2D
template<typename RealType> struct Hasher<2, RealType> {

  // typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef KeyTraits::Key CoordHash;

  static unsigned  num1dbits()                { return 31U; }
  //static unsigned  num1dbits()                { return 30U; }
  static CoordHash coordMax()                 { return (1ULL << num1dbits()) - 1ULL; }
  static CoordHash outerFlag()                { return (1ULL << 63); }
  static CoordHash xmask()                    { return (1ULL << num1dbits()) - 1ULL; }
  static CoordHash ymask()                    { return xmask() << num1dbits(); }
  static CoordHash qxval(const CoordHash val) { return (val & xmask()); }
  static CoordHash qyval(const CoordHash val) { return (val & ymask()) >> num1dbits(); }
  static CoordHash hash(const CoordHash x, const CoordHash y) { return x + (y << num1dbits()); }

  // Hash a 2 position
  static CoordHash hashPosition(const RealType* pos,
				const RealType* xlow_inner,
				const RealType* xhigh_inner,
				const RealType* xlow_outer,
				const RealType* xhigh_outer,
				const RealType minTol = RealType(0)) {
    POLY_ASSERT(xlow_outer[0] <= xlow_inner[0] and
                xlow_outer[1] <= xlow_inner[1]);
    POLY_ASSERT(xhigh_outer[0] >= xhigh_inner[0] and
                xhigh_outer[1] >= xhigh_inner[1]);
    POLY_ASSERT(xlow_inner[0] <= xhigh_inner[0] and
                xlow_inner[1] <= xhigh_inner[1]);
    POLY_ASSERT2(pos[0] >= xlow_outer[0] and pos[0] <= xhigh_outer[0] and
                 pos[1] >= xlow_outer[1] and pos[1] <= xhigh_outer[1],
                 "(" << pos[0] << " " << pos[1] << ") ("
                 << xlow_outer[0] << " " << xlow_outer[1] << ") ("
                 << xhigh_outer[0] << " " << xhigh_outer[1] << ")");

    // Decide the bounding box we're using.
    CoordHash result = 0ULL;
    const RealType *xlow, *xhigh;
    if (pos[0] < xlow_inner[0] or pos[0] > xhigh_inner[0] or
        pos[1] < xlow_inner[1] or pos[1] > xhigh_inner[1]) {
      xlow = xlow_outer;
      xhigh = xhigh_outer;
      result += (1ULL << 63);
    } else {
      xlow = xlow_inner;
      xhigh = xhigh_inner;
    }

    // Quantize away.
    const RealType dx[2] = {std::max(RealType((xhigh[0] - xlow[0])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon())),
                            std::max(RealType((xhigh[1] - xlow[1])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon()))};
    // const RealType delta = std::min(dx[0], dx[1]);
    result += hash(CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[0] - xlow[0])/dx[0]))),
                   CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[1] - xlow[1])/dx[1]))));
    // result += hash(CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[0] - xlow[0])/delta))),
    //                CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[1] - xlow[1])/delta))));
    return result;
  }

  // Unhash a 2 position.
  static void unhashPosition(RealType* pos,
                             const RealType* xlow_inner,
                             const RealType* xhigh_inner,
                             const RealType* xlow_outer,
                             const RealType* xhigh_outer,
                             const CoordHash hashedPosition,
                             const RealType minTol = RealType(0)) {
    POLY_ASSERT(xlow_outer[0] <= xlow_inner[0] and
                xlow_outer[1] <= xlow_inner[1]);
    POLY_ASSERT(xhigh_outer[0] >= xhigh_inner[0] and
                xhigh_outer[1] >= xhigh_inner[1]);
    POLY_ASSERT(xlow_inner[0] <= xhigh_inner[0] and
                xlow_inner[1] <= xhigh_inner[1]);

    // Decide the bounding box we're using.
    const RealType *xlow, *xhigh;
    if (hashedPosition >= (1ULL << 63)) {
      xlow = xlow_outer;
      xhigh = xhigh_outer;
    } else {
      xlow = xlow_inner;
      xhigh = xhigh_inner;
    }

    // Extract the position (for the center of the cell).
    const RealType dx[2] = {std::max(RealType((xhigh[0] - xlow[0])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon())),
                            std::max(RealType((xhigh[1] - xlow[1])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon()))};
    // const RealType delta = std::min(dx[0], dx[1]);
    // pos[0] = std::max(xlow[0], std::min(xhigh[0], RealType(xlow[0] + RealType(( (hashedPosition & xmask())        + 0.5)*dx[0]))));
    // pos[1] = std::max(xlow[1], std::min(xhigh[1], RealType(xlow[1] + RealType((((hashedPosition & ymask()) >> num1dbits()) + 0.5)*dx[1]))));
    pos[0] = std::max(xlow[0], std::min(xhigh[0], RealType(xlow[0] + RealType((qxval(hashedPosition) + 0.5)*dx[0]))));
    pos[1] = std::max(xlow[1], std::min(xhigh[1], RealType(xlow[1] + RealType((qyval(hashedPosition) + 0.5)*dx[1]))));
    // pos[0] = std::max(xlow[0], std::min(xhigh[0], RealType(xlow[0] + RealType((qxval(hashedPosition) + 0.5)*delta))));
    // pos[1] = std::max(xlow[1], std::min(xhigh[1], RealType(xlow[1] + RealType((qyval(hashedPosition) + 0.5)*delta))));

    // Post-conditions.
    POLY_ASSERT2(pos[0] >= xlow[0] and pos[0] <= xhigh[0] and
                 pos[1] >= xlow[1] and pos[1] <= xhigh[1],
                 "(" << pos[0] << " " << pos[1] << ") ("
                 << xlow[0] << " " << xlow[1] << ") ("
                 << xhigh[0] << " " << xhigh[1] << ")");
  }

  // Return hashed integer as a hashed 2 position.
  static void hashedPosition(CoordHash* pos,
			     const CoordHash hashedPosition) {
    pos[0] = qxval(hashedPosition);
    pos[1] = qyval(hashedPosition);
  }
};

// 3D
template<typename RealType> struct Hasher<3, RealType> {

  // typedef DimensionTraits<3, RealType>::CoordHash CoordHash;
  typedef KeyTraits::Key CoordHash;

  static unsigned  num1dbits()                { return 21U; }
  static CoordHash coordMax()                 { return (1ULL << num1dbits()) - 1ULL; }
  static CoordHash outerFlag()                { return (1ULL << 63); }
  static CoordHash xmask()                    { return (1ULL << num1dbits()) - 1ULL; }
  static CoordHash ymask()                    { return xmask() << num1dbits(); }
  static CoordHash zmask()                    { return xmask() << (2*num1dbits()); }
  static CoordHash qxval(const CoordHash val) { return (val & xmask()); }
  static CoordHash qyval(const CoordHash val) { return (val & ymask()) >> num1dbits(); }
  static CoordHash qzval(const CoordHash val) { return (val & zmask()) >> (2*num1dbits()); }
  static CoordHash hash(const CoordHash x, const CoordHash y, const CoordHash z) { return x + (y << num1dbits()) + (z << (2*num1dbits())); }

  // Hash a 3 position
  static CoordHash hashPosition(const RealType* pos,
				const RealType* xlow_inner,
				const RealType* xhigh_inner,
				const RealType* xlow_outer,
				const RealType* xhigh_outer,
				const RealType minTol = RealType(0)) {
    POLY_ASSERT(xlow_outer[0] <= xlow_inner[0] and
                xlow_outer[1] <= xlow_inner[1] and
                xlow_outer[2] <= xlow_inner[2]);
    POLY_ASSERT(xhigh_outer[0] >= xhigh_inner[0] and
                xhigh_outer[1] >= xhigh_inner[1] and
                xhigh_outer[2] >= xhigh_inner[2]);
    POLY_ASSERT(xlow_inner[0] < xhigh_inner[0] and
                xlow_inner[1] < xhigh_inner[1] and
                xlow_inner[2] < xhigh_inner[2]);
    POLY_ASSERT2(pos[0] >= xlow_outer[0] and pos[0] <= xhigh_outer[0] and
                 pos[1] >= xlow_outer[1] and pos[1] <= xhigh_outer[1] and
                 pos[2] >= xlow_outer[2] and pos[2] <= xhigh_outer[2],
                 "(" << pos[0] << " " << pos[1] << " " << pos[2] << ") ("
                 << xlow_outer[0] << " " << xlow_outer[1] << " " << xlow_outer[2] << ") ("
                 << xhigh_outer[0] << " " << xhigh_outer[1] << " " << xhigh_outer[2] << ")");

    // Decide the bounding box we're using.
    CoordHash result = 0ULL;
    const RealType *xlow, *xhigh;
    if (pos[0] < xlow_inner[0] or pos[0] > xhigh_inner[0] or
        pos[1] < xlow_inner[1] or pos[1] > xhigh_inner[1] or
        pos[2] < xlow_inner[2] or pos[2] > xhigh_inner[2]) {
      xlow = xlow_outer;
      xhigh = xhigh_outer;
      result += (1ULL << 63);
    } else {
      xlow = xlow_inner;
      xhigh = xhigh_inner;
    }

    // Quantize away.
    const RealType dx[3] = {std::max(RealType((xhigh[0] - xlow[0])/coordMax()), std::max(std::numeric_limits<RealType>::epsilon(), minTol)),
                            std::max(RealType((xhigh[1] - xlow[1])/coordMax()), std::max(std::numeric_limits<RealType>::epsilon(), minTol)),
                            std::max(RealType((xhigh[2] - xlow[2])/coordMax()), std::max(std::numeric_limits<RealType>::epsilon(), minTol))};
    result += hash(CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[0] - xlow[0])/dx[0]))),
                   CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[1] - xlow[1])/dx[1]))),
                   CoordHash(std::min(coordMax(), CoordHash(std::max(RealType(0), pos[2] - xlow[2])/dx[2]))));
    return result;
  }

  // Unhash a 3 position.
  static void unhashPosition(RealType* pos,
                             const RealType* xlow_inner,
                             const RealType* xhigh_inner,
                             const RealType* xlow_outer,
                             const RealType* xhigh_outer,
                             const CoordHash hashedPosition,
                             const RealType minTol = RealType(0)) {
    POLY_ASSERT(xlow_outer[0] <= xlow_inner[0] and
                xlow_outer[1] <= xlow_inner[1] and
                xlow_outer[2] <= xlow_inner[2]);
    POLY_ASSERT(xhigh_outer[0] >= xhigh_inner[0] and
                xhigh_outer[1] >= xhigh_inner[1] and
                xhigh_outer[2] >= xhigh_inner[2]);
    POLY_ASSERT(xlow_inner[0] < xhigh_inner[0] and
                xlow_inner[1] < xhigh_inner[1] and
                xlow_inner[2] < xhigh_inner[2]);

    // Decide the bounding box we're using.
    const RealType *xlow, *xhigh;
    if (hashedPosition >= (1ULL << 63)) {
      xlow = xlow_outer;
      xhigh = xhigh_outer;
    } else {
      xlow = xlow_inner;
      xhigh = xhigh_inner;
    }

    // Extract the position (for the center of the cell).
    const RealType dx[3] = {std::max(RealType((xhigh[0] - xlow[0])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon())),
                            std::max(RealType((xhigh[1] - xlow[1])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon())),
                            std::max(RealType((xhigh[2] - xlow[2])/coordMax()), std::max(minTol, std::numeric_limits<RealType>::epsilon()))};
    pos[0] = std::max(xlow[0], std::min(xhigh[0], RealType(xlow[0] + RealType(( (hashedPosition & xmask())        + 0.5)*dx[0]))));
    pos[1] = std::max(xlow[1], std::min(xhigh[1], RealType(xlow[1] + RealType((((hashedPosition & ymask()) >> num1dbits()) + 0.5)*dx[1]))));
    pos[2] = std::max(xlow[2], std::min(xhigh[2], RealType(xlow[2] + RealType((((hashedPosition & zmask()) >> (2*num1dbits())) + 0.5)*dx[2]))));

    // Post-conditions.
    POLY_ASSERT2(pos[0] >= xlow[0] and pos[0] <= xhigh[0] and
                 pos[1] >= xlow[1] and pos[1] <= xhigh[1] and
                 pos[2] >= xlow[2] and pos[2] <= xhigh[2],
                 "(" << pos[0] << " " << pos[1] << " " << pos[2] << ") ("
                 << xlow[0] << " " << xlow[1] << " " << xlow[2] << ") ("
                 << xhigh[0] << " " << xhigh[1] << " " << xhigh[2] << ")");
  }
};

//------------------------------------------------------------------------------
// Distance between points.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct DistanceFunctor;

// 2D
template<typename RealType> 
struct DistanceFunctor<2, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    const RealType dx = a[0] - b[0];
    const RealType dy = a[1] - b[1];
    //return std::sqrt(dx*dx + dy*dy);
    return sqrt(dx*dx + dy*dy);
  }
};

// 3D
template<typename RealType> 
struct DistanceFunctor<3, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    const RealType dx = a[0] - b[0];
    const RealType dy = a[1] - b[1];
    const RealType dz = a[2] - b[2];
    //return std::sqrt(dx*dx + dy*dy + dz*dz);
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
RealType
distance(const RealType* a, const RealType* b) {
  return DistanceFunctor<Dimension, RealType>::impl(a, b);
}

//------------------------------------------------------------------------------
// Compute the unit vector.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct UnitVectorFunctor;

// 2D
template<typename RealType> 
struct UnitVectorFunctor<2, RealType> {
  static void impl(RealType* a) {
    const RealType mag = std::max(1.0e-100, sqrt(a[0]*a[0] + a[1]*a[1]));
    a[0] /= mag;
    a[1] /= mag;
  }
};

// 3D
template<typename RealType> 
struct UnitVectorFunctor<3, RealType> {
  static void impl(RealType* a) {
    const RealType mag = std::max(1.0e-100, sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
void
unitVector(RealType* a) {
  UnitVectorFunctor<Dimension, RealType>::impl(a);
}

//------------------------------------------------------------------------------
// Dot product.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct DotFunctor;

// 2D
template<typename RealType> 
struct DotFunctor<2, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    return a[0]*b[0] + a[1]*b[1];
  }
};

// 3D
template<typename RealType> 
struct DotFunctor<3, RealType> {
  static RealType impl(const RealType* a, const RealType* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
RealType
dot(const RealType* a, const RealType* b) {
  return DotFunctor<Dimension, RealType>::impl(a, b);
}

//------------------------------------------------------------------------------
// Cross product.
//------------------------------------------------------------------------------
// Functor's for use in partial specializations.
template<int Dimension, typename RealType> struct CrossFunctor;

// 2D
template<typename RealType> 
struct CrossFunctor<2, RealType> {
  static void impl(const RealType* a, const RealType* b, RealType* c) {
    c[0] = 0.0; c[1] = 0.0; c[2] = a[0]*b[1] - a[1]*b[0];
  }
};

// 3D
template<typename RealType> 
struct CrossFunctor<3, RealType> {
  static void impl(const RealType* a, const RealType* b, RealType* c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
};

//****************************************
// Functional interface.
// This is the one you call!
//****************************************
template<int Dimension, typename RealType>
void
cross(const RealType* a, const RealType* b, RealType* c) {
  CrossFunctor<Dimension, RealType>::impl(a, b, c);
}

//------------------------------------------------------------------------------
// Determine if the given points are collinear to some accuracy.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
collinear(const RealType* a, const RealType* b, const RealType* c, const RealType tol) {
  double ab[Dimension], ac[Dimension], abmag = 0.0, acmag = 0.0;
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] = b[j] - a[j];
    ac[j] = c[j] - a[j];
    abmag += ab[j]*ab[j];
    acmag += ac[j]*ac[j];
  }
  if (abmag < tol or acmag < tol) return true;
  abmag = sqrt(abmag);
  acmag = sqrt(acmag);
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] /= abmag;
    ac[j] /= acmag;
  }
  return std::abs(std::abs(dot<Dimension, double>(ab, ac)) - 1.0) < tol;
}

//------------------------------------------------------------------------------
// Test if point c is between a & b.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
between(const RealType* a, const RealType* b, const RealType* c, const RealType tol) {

  // Compute a bunch of distances.
  RealType ab[Dimension], ac[Dimension], bc[Dimension], ac_mag2 = 0, bc_mag2 = 0, ab_mag2 = 0;
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] = b[j] - a[j];
    ac[j] = c[j] - a[j];
    bc[j] = c[j] - b[j];
    ab_mag2 += ab[j]*ab[j];
    ac_mag2 += ac[j]*ac[j];
    bc_mag2 += bc[j]*bc[j];
  }

  // If c is equal to either endpoint, we count that as between.
  if (std::min(ac_mag2, bc_mag2) <= tol*ab_mag2) return true;

  // If a & b are the same point, but not c it's outside.
  if (ab_mag2 <= tol) return false;

  // The points are distinct.
  const RealType thpt = dot<Dimension, RealType>(ab, ac);
  return ((thpt > 0) and (std::abs(thpt*thpt - ab_mag2*ac_mag2) < tol*ab_mag2) and (ac_mag2 <= ab_mag2));
}


//------------------------------------------------------------------------------
// Determine if a cloud of points is collinear to some accuracy
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
collinear(const std::vector<RealType> points, const RealType tol) {
  POLY_ASSERT(points.size() % Dimension == 0);
  bool isCollinear = true;
  int i;
  if (points.size()/Dimension > 1) {
    i = Dimension;
    while (isCollinear and i != points.size()/Dimension) {
      isCollinear *= collinear<Dimension, RealType>(&points[0],
						    &points[Dimension],
						    &points[Dimension*i],
						    tol);
      ++i;
    }
  }
  return isCollinear;
}


//------------------------------------------------------------------------------
// Find the closest point on a line segment (2D).
//------------------------------------------------------------------------------
template<typename RealType> 
void
closestPointOnSegment2D(const RealType* point, 
                        const RealType* s1,
                        const RealType* s2,
                        RealType* result) {
  RealType shat[2] = {s2[0] - s1[0], s2[1] - s1[1]};
  const RealType seglength = sqrt(shat[0]*shat[0] + shat[1]*shat[1]);
  if (seglength < 1.0e-10) {
    result[0] = 0.5*(s1[0] + s2[0]);
    result[1] = 0.5*(s1[1] + s2[1]);
  } else {
    shat[0] /= seglength;
    shat[1] /= seglength;
    const RealType s1p[2] = {point[0] - s1[0], point[1] - s1[1]};
    const RealType ptest = std::max(0.0, std::min(seglength, geometry::dot<2, RealType>(s1p, shat)));
    result[0] = s1[0] + ptest*shat[0];
    result[1] = s1[1] + ptest*shat[1];
  }
}


//------------------------------------------------------------------------------
// Determine if the point lies inside unclosed polygon determined by vertices
//------------------------------------------------------------------------------
template<typename RealType> 
bool
pointInPolygon(const RealType* point,
	       const unsigned numVertices,
	       const RealType* vertices) {
  unsigned i,j;
  bool result = false;
  for (i = 0, j = numVertices-1; i < numVertices; j = i++) {
    POLY_ASSERT(i < numVertices and j < numVertices);
    if ( ((vertices[2*i+1] > point[1]) != (vertices[2*j+1] > point[1])) and
	 ( point[0] < (vertices[2*i] +
		       (vertices[2*j  ] - vertices[2*i  ]) *
		       (point[1]        - vertices[2*i+1]) /
		       (vertices[2*j+1] - vertices[2*i+1]))) )
      result = not result;
  }
  return result;
}


//------------------------------------------------------------------------------
// Determine if the point lies on a polygon boundary
//------------------------------------------------------------------------------
template<typename RealType> 
bool
pointOnPolygon(const RealType* point,
	       const unsigned numVertices,
	       const RealType* vertices) {
  unsigned i=0, j;
  bool result = false;
  while (i < numVertices and not result) {
    j = (i+1) % numVertices;
    POLY_ASSERT(i < numVertices and j < numVertices);
    result = ( (std::min(vertices[2*i  ],vertices[2*j  ]) <= point[0]) and
	       (std::max(vertices[2*i  ],vertices[2*j  ]) >= point[0]) and
	       (std::min(vertices[2*i+1],vertices[2*j+1]) <= point[1]) and
	       (std::max(vertices[2*i+1],vertices[2*j+1]) >= point[1]) and
	       collinear<2, RealType>(&vertices[2*i], &vertices[2*j], point, 1.0e-10) );
    ++i;
  }
  return result;
}


//------------------------------------------------------------------------------
// // OLD IMPLEMENTATION
// template<typename RealType> 
// bool
// withinPolygon2D(const RealType* point,
//                 const unsigned numVertices,
//                 const RealType* vertices) {
//    unsigned i,j;
//    bool result = false;
//    for (i = 0, j = numVertices; i < numVertices; j = i++) {
//       // if (collinear<2,RealType>(&vertices[2*i], &vertices[2*j], point, 1.0e-10)){
//       //    return true;
//       // }
//       if (((vertices[2*i+1]<= point[1] and vertices[2*j+1] >= point[1])  or
//            (vertices[2*j+1]<= point[1] and vertices[2*i+1] >= point[1])) and
//           (vertices[2*i]   <= point[0] or  vertices[2*j]   <= point[0]) ) {
//          if (collinear<2,RealType>(&vertices[2*i], &vertices[2*j], point, 1.0e-10)){
//             return true;
//          }else{
//             result ^= ( vertices[2*i] + 
//                         ( point[1]        - vertices[2*i+1] ) /
//                         ( vertices[2*j+1] - vertices[2*i+1] ) *
//                         ( vertices[2*j]   - vertices[2*i]   ) < point[0] );
//          }
//       }
//    }
//    return result;
// }

//------------------------------------------------------------------------------
// Find the global bounding box for a set of coordinates.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
computeBoundingBox(const RealType* pos,
                   const unsigned pos_size,
                   const bool globalReduce,
                   RealType xmin[Dimension],
                   RealType xmax[Dimension]) {
   POLY_ASSERT(pos_size % Dimension == 0);
   for (unsigned j = 0; j != Dimension; ++j) {
      xmin[j] = std::numeric_limits<RealType>::max();
      xmax[j] = (std::numeric_limits<RealType>::is_signed ? -xmin[j] : std::numeric_limits<RealType>::min());
   }
   const unsigned n = pos_size/Dimension;
   for (unsigned i = 0; i != n; ++i) {
     for (unsigned j = 0; j != Dimension; ++j) {
       xmin[j] = std::min(xmin[j], pos[Dimension*i + j]);
       xmax[j] = std::max(xmax[j], pos[Dimension*i + j]);
     }
   }
#ifdef HAVE_MPI
   if (globalReduce) {
     for (unsigned j = 0; j != Dimension; ++j) {
       xmin[j] = allReduce(xmin[j], MPI_MIN, MPI_COMM_WORLD);
       xmax[j] = allReduce(xmax[j], MPI_MAX, MPI_COMM_WORLD);
     }
   }
#endif
}




template<int Dimension, typename RealType>
inline
void
computeBoundingBox(const std::vector<RealType>& pos,
                   const bool globalReduce,
                   RealType xmin[Dimension],
                   RealType xmax[Dimension]) {
  computeBoundingBox<Dimension, RealType>(&pos.front(),
                                          pos.size(),
                                          globalReduce,
                                          xmin,
                                          xmax);
}



template<int Dimension, typename RealType>
inline
void
expandBoundingBox(const RealType* pos,
		  const unsigned pos_size,
		  const bool globalReduce,
		  RealType xmin[Dimension],
		  RealType xmax[Dimension]) {
  RealType low[Dimension], high[Dimension];
  computeBoundingBox<Dimension, RealType>(pos,
					  pos_size,
					  globalReduce,
					  low,
					  high);
  const unsigned n = pos_size/Dimension;
  for (unsigned i = 0; i != n; ++i) {
    for (unsigned j = 0; j != Dimension; ++j) {
      xmin[j] = std::min(xmin[j], low [j]);
      xmax[j] = std::max(xmax[j], high[j]);
    }
  }
#ifdef HAVE_MPI
  if (globalReduce) {
    for (unsigned j = 0; j != Dimension; ++j) {
      xmin[j] = allReduce(xmin[j], MPI_MIN, MPI_COMM_WORLD);
      xmax[j] = allReduce(xmax[j], MPI_MAX, MPI_COMM_WORLD);
    }
  }
#endif
}


//------------------------------------------------------------------------------
// Find the ray-plane intersection, return true if there's a valid answer and
// false if not.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// p_plane : point on plane.
// n_plane : unit normal of plane.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point (if one is possible).
//------------------------------------------------------------------------------
template<typename RealType> 
bool
rayPlaneIntersection(const RealType* p_ray,
                     const RealType* n_ray,
                     const RealType* p_plane,
                     const RealType* n_plane,
                     const RealType& tol,
                     RealType* result) {
  const RealType ndots = dot<3, RealType>(n_ray, n_plane);
  const RealType delta[3] = {p_plane[0] - p_ray[0], p_plane[1] - p_ray[1], p_plane[2] - p_ray[2]};
  if (std::abs(ndots) < tol) {
    // Ray is parallel to the plane.
    if (std::abs(dot<3, RealType>(p_ray, delta)) < tol) {
      // Ray is in the plane, just choose the origin point.
      result[0] = p_ray[0]; result[1] = p_ray[1]; result[2] = p_ray[2];
      return true;
    } else {
      // Ray is out of the plane.
      return false;
    }
  }

  // The ray is not parallel to the plane, so check for intersection.
  const RealType d = dot<3, RealType>(delta, n_plane)/ndots;
  if (d >= 0.0) {
    // Yep, they intersect.
    result[0] = p_ray[0] + d*n_ray[0];
    result[1] = p_ray[1] + d*n_ray[1];
    result[2] = p_ray[2] + d*n_ray[2];
    return true;
  } else {
    // The intersection is the wrong direction from the origin of the ray,
    // so they don't intersect.
    return false;
  }
}

//------------------------------------------------------------------------------
// Find the ray-sphere intersection, assuming the ray's origin is inside the 
// sphere.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// p_sphere: center of the sphere.
// r_sphere: radius of the sphere.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point
//------------------------------------------------------------------------------
template<typename RealType> 
bool
raySphereIntersection(const RealType* p_ray,
                      const RealType* n_ray,
                      const RealType* p_sphere,
                      const RealType r_sphere,
                      const RealType& tol,
                      RealType* result) {
  POLY_ASSERT2((distance<3, RealType>(p_ray, p_sphere) <= r_sphere), 
               "(" << p_ray[0] << " " << p_ray[1] << " " << p_ray[2] 
               << ") (" << p_sphere[0] << " " << p_sphere[1] << " " << p_sphere[2] << ")" << " : " 
               << (distance<3, RealType>(p_ray, p_sphere)) << " " << r_sphere);
  POLY_ASSERT(std::abs(n_ray[0]*n_ray[0] + n_ray[1]*n_ray[1] + n_ray[2]*n_ray[2] - 1.0) < 1.0e-10);
  const RealType rs0[3] = {p_ray[0] - p_sphere[0],
                           p_ray[1] - p_sphere[1],
                           p_ray[2] - p_sphere[2]};
  const RealType b = 2.0*dot<3, RealType>(n_ray, rs0);
  const RealType c = dot<3, RealType>(rs0, rs0) - r_sphere*r_sphere;
  RealType d = b*b - 4.0*c;
  POLY_ASSERT(d >= -tol);
  if (d < tol) {
    // Glancing intersection at the edge of the sphere -- it must be the origin of the ray!
    result[0] = p_ray[0]; result[1] = p_ray[1]; result[2] = p_ray[2];
    return true;
  } else {
    //d = std::sqrt(d);
    d = sqrt(d);
    const RealType t = 0.5*std::max(-b - d, -b + d);
    POLY_ASSERT(t >= 0.0);
    result[0] = p_ray[0] + t*n_ray[0];
    result[1] = p_ray[1] + t*n_ray[1];
    result[2] = p_ray[2] + t*n_ray[2];
  }
  POLY_ASSERT(std::abs(distance<3, RealType>(result, p_sphere) - r_sphere) < tol);
  return true;
}

//------------------------------------------------------------------------------
// Find the ray-circle intersection, assuming the ray's origin is inside the 
// circle.
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// p_circle: center of the circle.
// r_circle: radius of the circle.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point
//------------------------------------------------------------------------------
template<typename RealType> 
bool
rayCircleIntersection(const RealType* p_ray,
                      const RealType* n_ray,
                      const RealType* p_circle,
                      const RealType r_circle,
                      const RealType& tol,
                      RealType* result) {
  if (distance<2, RealType>(p_ray, p_circle) > r_circle) {
    result[0] = p_ray[0];  result[1] = p_ray[1];
    return true;
  }
  // POLY_ASSERT2((distance<2, RealType>(p_ray, p_circle) <= r_circle), 
  //              "(" << p_ray[0] << " " << p_ray[1] << ") (" 
  //              << p_circle[0] << " " << p_circle[1] << ")" << " : " 
  //              << (distance<2, RealType>(p_ray, p_circle)) << " " << r_circle);
  POLY_ASSERT(std::abs(n_ray[0]*n_ray[0] + n_ray[1]*n_ray[1] - 1.0) < 1.0e-10);
  const RealType rs0[2] = {p_ray[0] - p_circle[0],
                           p_ray[1] - p_circle[1]};
  const RealType b = dot<2, RealType>(n_ray, rs0);
  RealType d = b*b - dot<2, RealType>(rs0, rs0) + r_circle*r_circle;
  POLY_ASSERT(d >= -tol);
  if (d < tol) {
    // Glancing intersection at the edge of the circle -- it must be the origin of the ray!
    result[0] = p_ray[0]; result[1] = p_ray[1];
    return true;
  } else {
    d = sqrt(d);
    const RealType t = std::max(-b - d, -b + d);
    POLY_ASSERT(t >= 0.0);
    result[0] = p_ray[0] + t*n_ray[0];
    result[1] = p_ray[1] + t*n_ray[1];
  }
  //POLY_ASSERT(std::abs(distance<2, RealType>(result, p_circle) - r_circle) < tol);
  return true;
}

//------------------------------------------------------------------------------
// Find the closest point of intersection for ray and a box.
// Returns an int indicating which plane is intersected:
// -1 : no intersection
//  0 : x = xmin plane
//  1 : x = xmax plane
//  2 : y = ymin plane
//  3 : y = ymax plane
//  4 : z = zmin plane
//  5 : z = zmax plane
//
// Arguments:
//   p_ray : point origin of ray.
//   n_ray : unit normal in direction of ray.
// box_min : minimum coordinate of box.
// box_max : maximum coordinate of box.
//     tol : the tolerance for zero (check for parallel lines and such)
//  result : the intersection point (if one is possible).
//------------------------------------------------------------------------------
template<typename RealType> 
int
rayBoxIntersection(RealType* p_ray,
                   RealType* n_ray,
                   RealType* box_min,
                   RealType* box_max,
                   const RealType& tol,
                   RealType* result) {

  // Prepare.
  int iplane = -1;
  RealType dist, minDist = std::numeric_limits<RealType>::max(), candidate[3], n_plane[3];

  // x=xmin plane.
  n_plane[0] = 1.0; n_plane[1] = 0.0; n_plane[2] = 0.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    iplane = 0;
    dist = distance<3, RealType>(p_ray, candidate);
    minDist = dist;
    result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
  }

  // x=xmax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      iplane = 1;
      minDist = dist;
      result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
    }
  }

  // y=ymin plane.
  n_plane[0] = 0.0; n_plane[1] = 1.0; n_plane[2] = 0.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      iplane = 2;
      minDist = dist;
      result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
    }
  }

  // y=ymax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      iplane = 3;
      minDist = dist;
      result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
    }
  }

  // z=zmin plane.
  n_plane[0] = 0.0; n_plane[1] = 0.0; n_plane[2] = 1.0;
  if (rayPlaneIntersection(p_ray, n_ray, box_min, n_plane, tol, candidate)) {
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      iplane = 4;
      minDist = dist;
      result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
    }
  }

  // z=zmax plane.
  if (rayPlaneIntersection(p_ray, n_ray, box_max, n_plane, tol, candidate)) {
    dist = distance<3, RealType>(p_ray, candidate);
    if (dist < minDist) {
      iplane = 5;
      minDist = dist;
      result[0] = candidate[0]; result[1] = candidate[1]; result[2] = candidate[2];
    }
  }

  return iplane;
}

//------------------------------------------------------------------------------
// Compute the centroid of the mesh cell.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
void
computeCellCentroid(const Tessellation<Dimension, RealType>& mesh,
                    const unsigned ci,
                    RealType* ccent) {
  POLY_ASSERT(ci < mesh.cells.size());
  const unsigned nf = mesh.cells[ci].size();
  unsigned i, j, k;
  for (j = 0; j != Dimension; ++j) ccent[j] = 0.0;
  std::vector<unsigned> uniqueNodes;
  for (k = 0; k != nf; ++k) {
    const std::vector<unsigned>& faceNodes = mesh.faces[internal::positiveID(mesh.cells[ci][k])];
    std::copy(faceNodes.begin(), faceNodes.end(), std::back_inserter(uniqueNodes));
  }

  // Now reduce to the unique vertices.
  std::sort(uniqueNodes.begin(), uniqueNodes.end());
  std::vector<unsigned>::iterator endItr = std::unique(uniqueNodes.begin(), uniqueNodes.end());
  for (std::vector<unsigned>::iterator itr = uniqueNodes.begin(); itr != endItr; ++itr) {
    i = *itr;
    POLY_ASSERT2(i < mesh.nodes.size()/Dimension, i << " " << mesh.nodes.size()/Dimension);
    for (j = 0; j != Dimension; ++j) ccent[j] += mesh.nodes[Dimension*i + j];
  }
  const unsigned n = std::distance(uniqueNodes.begin(), endItr);
  POLY_ASSERT(n > 0);
  for (j = 0; j != Dimension; ++j) ccent[j] /= n;
}

//------------------------------------------------------------------------------
// Compute the cell centroid and signed area.
// Taken from http://www.wikipedia.org/wiki/Centroid
//------------------------------------------------------------------------------
template<typename RealType>
void
computeCellCentroidAndSignedArea(const Tessellation<2, RealType>& mesh,
				 const unsigned ci,
				 const RealType& tol,
				 RealType* ccent,
				 RealType& area) {
  POLY_ASSERT(ci < mesh.cells.size());
  unsigned iface, n0, n1;
  RealType d, x0, x1, y0, y1;
  ccent[0] = 0.0; ccent[1] = 0.0;  area = 0.0;
  for (std::vector<int>::const_iterator itr = mesh.cells[ci].begin();
       itr != mesh.cells[ci].end(); ++itr) {
    iface = (*itr < 0) ? ~(*itr) : *itr;
    POLY_ASSERT(iface < mesh.faces.size());
    POLY_ASSERT(mesh.faces[iface].size() == 2);
    n0 = (*itr < 0) ? mesh.faces[iface][1] : mesh.faces[iface][0];
    n1 = (*itr < 0) ? mesh.faces[iface][0] : mesh.faces[iface][1];
    POLY_ASSERT(n0 < mesh.nodes.size()/2 and n1 < mesh.nodes.size()/2);
    POLY_ASSERT(n0 != n1);
    x0 = mesh.nodes[2*n0];  y0 = mesh.nodes[2*n0+1];
    x1 = mesh.nodes[2*n1];  y1 = mesh.nodes[2*n1+1];
    d = x0*y1 - y0*x1;
    area     += d;
    ccent[0] += d*(x0+x1);
    ccent[1] += d*(y0+y1);
  }
  POLY_ASSERT(std::abs(area) > tol);
  area     /= 2.0;
  ccent[0] /= (6*area);
  ccent[1] /= (6*area);
}

//------------------------------------------------------------------------------
// Compute the centroid and unit normal of a Tessellation face.
//------------------------------------------------------------------------------
template<typename RealType>
void
computeFaceCentroidAndNormal(const Tessellation<3, RealType>& mesh,
                             const unsigned fi,
                             RealType* fcent,
                             RealType* fhat) {
  POLY_ASSERT(fi < mesh.faces.size());
  const unsigned n = mesh.faces[fi].size();
  POLY_ASSERT(n >= 3);
  unsigned i, ni;
  std::vector<unsigned> verts;
  const double degeneracy = 1.0e-10;

  // Compute the centroid, and look for three vertices in the face that are 
  // not collinear.
  fcent[0] = 0.0; fcent[1] = 0.0; fcent[2] = 0.0;
  for (i = 0; i != n; ++i) {
    ni = mesh.faces[fi][i];
    POLY_ASSERT(ni < mesh.nodes.size()/3);
    fcent[0] += mesh.nodes[3*ni];
    fcent[1] += mesh.nodes[3*ni+1];
    fcent[2] += mesh.nodes[3*ni+2];
    if (verts.size() < 2 or
        (verts.size() == 2 and not collinear<3, RealType>(&mesh.nodes[3*verts[0]],
                                                          &mesh.nodes[3*verts[1]],
                                                          &mesh.nodes[3*ni],
                                                          degeneracy))) verts.push_back(ni);
  }
  POLY_ASSERT(n > 0);
  fcent[0] /= n; fcent[1] /= n; fcent[2] /= n;
  
  // Now we can compute the unit normal.
  POLY_ASSERT2(verts.size() == 3, verts.size());
  RealType ab[3], ac[3];
  ab[0] = mesh.nodes[3*verts[1]  ] - mesh.nodes[3*verts[0]  ];
  ab[1] = mesh.nodes[3*verts[1]+1] - mesh.nodes[3*verts[0]+1];
  ab[2] = mesh.nodes[3*verts[1]+2] - mesh.nodes[3*verts[0]+2];
  ac[0] = mesh.nodes[3*verts[2]  ] - mesh.nodes[3*verts[0]  ];
  ac[1] = mesh.nodes[3*verts[2]+1] - mesh.nodes[3*verts[0]+1];
  ac[2] = mesh.nodes[3*verts[2]+2] - mesh.nodes[3*verts[0]+2];
  cross<3, RealType>(ab, ac, fhat);
  unitVector<3, RealType>(fhat);
}

//------------------------------------------------------------------------------
// computeSquaredNorm
// Computed with Shewchuck's fast predicates
//------------------------------------------------------------------------------
inline
double
computeSquaredNorm(double* A) {
  double axax[2], ayay[2], tmp[2];
  scale_expansion(1, &A[0], A[0], axax);
  scale_expansion(1, &A[1], A[1], ayay);
  grow_expansion(1, &axax[1], ayay[1], tmp);
  return tmp[1];
}

//------------------------------------------------------------------------------
// computeCircumcenter
// Computed with Shewchuck's fast predicates
//------------------------------------------------------------------------------
inline
void
computeCircumcenter(double* A, double* B, double* C, double* X) {
  // double Anorm = computeSquaredNorm(A);
  // double Bnorm = computeSquaredNorm(B);
  // double Cnorm = computeSquaredNorm(C);
  double Anorm = A[0]*A[0] + A[1]*A[1];
  double Bnorm = B[0]*B[0] + B[1]*B[1];
  double Cnorm = C[0]*C[0] + C[1]*C[1];  
  double D = 2*orient2d(A,B,C);
  double a0[2] = {Anorm, A[1]},  a1[2] = {A[0], Anorm};
  double b0[2] = {Bnorm, B[1]},  b1[2] = {B[0], Bnorm};
  double c0[2] = {Cnorm, C[1]},  c1[2] = {C[0], Cnorm};
  X[0] = orient2d(a0,b0,c0)/D;
  X[1] = orient2d(a1,b1,c1)/D;
}


//------------------------------------------------------------------------------
// This function computes the circumcenter of a triangle with vertices
// A = (Ax, Ay), B = (Bx, By), and C = (Cx, Cy), and places the result 
// in X.
//------------------------------------------------------------------------------
inline
void 
computeCircumcenter2d(double* A, double* B, double* C, double* X) {
  // This solution was taken from Wikipedia's entry:
  // http://en.wikipedia.org/wiki/Circumscribed_circle
  double D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]));
  X[0] = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
          (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
  X[1] = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
          (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
}

//------------------------------------------------------------------------
// This function computes the circumcenter of a tetrahedron with vertices
// A = (Ax, Ay, Az), B = (Bx, By, Bz), C = (Cx, Cy, Cz), and D = (Dx, Dy, Dz).
// The result is returned in X.
// This solution was taken from http://en.wikipedia.org/wiki/Tetrahedron
//------------------------------------------------------------------------
inline
void
computeCircumcenter3d(double* A, double* B, double* C, double* D, double* X) {
  const double a[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]},
               b[3] = {C[0] - A[0], C[1] - A[1], C[2] - A[2]},
               c[3] = {D[0] - A[0], D[1] - A[1], D[2] - A[2]};
  const double a2 = dot<3, double>(a, a),
               b2 = dot<3, double>(b, b),
               c2 = dot<3, double>(c, c);
  double b_cross_c[3], c_cross_a[3], a_cross_b[3];
  cross<3, double>(b, c, b_cross_c);
  cross<3, double>(c, a, c_cross_a);
  cross<3, double>(a, b, a_cross_b);
  const double denom = 2.0*dot<3, double>(a, b_cross_c);
  POLY_ASSERT2(denom != 0.0,
               "A : (" << A[0] << " " << A[1] << " " << A[2] << ")\n"
               << "B : (" << B[0] << " " << B[1] << " " << B[2] << ")\n"
               << "C : (" << C[0] << " " << C[1] << " " << C[2] << ")\n"
               << "D : (" << D[0] << " " << D[1] << " " << D[2] << ")\n");
  X[0] = A[0] + (a2*b_cross_c[0] + b2*c_cross_a[0] + c2*a_cross_b[0])/denom;
  X[1] = A[1] + (a2*b_cross_c[1] + b2*c_cross_a[1] + c2*a_cross_b[1])/denom;
  X[2] = A[2] + (a2*b_cross_c[2] + b2*c_cross_a[2] + c2*a_cross_b[2])/denom;
}

//------------------------------------------------------------------------------
// Compute the circumcenter of a triangle in 3D.
// Same argument conventions as above.
// Based on the form presented at
// http://en.wikipedia.org/wiki/Circumscribed_circle
//------------------------------------------------------------------------------
inline
bool
computeTriangleCircumcenter3d(double* A, double* B, double* C, double* X) {
  const double degeneracy = 1.0e-10;
  if (collinear<3, double>(A, B, C, degeneracy)) return false;
  
  const double a[3] = {A[0] - C[0], 
                       A[1] - C[1], 
                       A[2] - C[2]};
  const double b[3] = {B[0] - C[0],
                       B[1] - C[1],
                       B[2] - C[2]};
  const double
    a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2],
    b2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  const double c[3] = {a2*b[0] - b2*a[0],
                       a2*b[1] - b2*a[1],
                       a2*b[2] - b2*a[2]};
  double d[3], e[3];
  cross<3, double>(a, b, d);
  cross<3, double>(c, d, e);
  const double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  POLY_ASSERT(d2 > 0.0);
  X[0] = 0.5*e[0]/d2 + C[0];
  X[1] = 0.5*e[1]/d2 + C[1];
  X[2] = 0.5*e[2]/d2 + C[2];
  return true;
}

//------------------------------------------------------------------------------
// Compute the centroid of a triangle in 2D.
//------------------------------------------------------------------------------
inline
void
computeTriangleCentroid2d(double* A, double* B, double* C, double* X) {
  X[0] = (A[0] + B[0] + C[0])/3.0;
  X[1] = (A[1] + B[1] + C[1])/3.0;
}

//------------------------------------------------------------------------------
// Compute the centroid of a triangle in 3D.
//------------------------------------------------------------------------------
inline
void
computeTriangleCentroid3d(double* A, double* B, double* C, double* X) {
  X[0] = (A[0] + B[0] + C[0])/3.0;
  X[1] = (A[1] + B[1] + C[1])/3.0;
  X[2] = (A[2] + B[2] + C[2])/3.0;
}

//------------------------------------------------------------------------------
// Compute the centroid of a tetrahedron in 3D.
//------------------------------------------------------------------------------
inline
void
computeTetCentroid(double* A, double* B, double* C, double* D, double* X) {
  X[0] = 0.25*(A[0] + B[0] + C[0] + D[0]);
  X[1] = 0.25*(A[1] + B[1] + C[1] + D[1]);
  X[2] = 0.25*(A[2] + B[2] + C[2] + D[2]);
}

//------------------------------------------------------------------------------
// Compute intersection of two 2D line segments if it exists
//------------------------------------------------------------------------------
template<typename RealType>
bool
segmentIntersection2D(const RealType* a, 
		      const RealType* b, 
		      const RealType* c, 
		      const RealType* d,
		      RealType* result,
                      const RealType tol = 1.0e-8) {
   RealType r1[2] = {b[0]-a[0] , b[1]-a[1]};  //direction vector of first segment
   RealType r2[2] = {d[0]-c[0] , d[1]-c[1]};  //direction vector of second segment
   RealType ca[2] = {c[0]-a[0] , c[1]-a[1]};
   RealType r1_cross_r2[3];
   cross<2,RealType>(r1, r2, r1_cross_r2);
   // Parallel segments
   if( r1_cross_r2[2] == 0 ) return false;
   RealType t1[3], t2[3];
   cross<2,RealType>(ca, r2, t1);
   cross<2,RealType>(ca, r1, t2);
   RealType p1 = t1[2]/r1_cross_r2[2];
   RealType p2 = t2[2]/r1_cross_r2[2];
   // The finite segments intersect
   if( -tol <= p1 and p1 <= 1+tol and -tol <= p2 and p2 <= 1+tol ) {
     if     (p1 < tol  ) {result[0] = a[0];  result[1] = a[1];}
     else if(p1 > 1-tol) {result[0] = b[0];  result[1] = b[1];}
     else if(p2 < tol  ) {result[0] = c[0];  result[1] = c[1];}
     else if(p2 > 1-tol) {result[0] = d[0];  result[1] = d[1];}
     else {
      result[0] = c[0] + p2*r2[0];
      result[1] = c[1] + p2*r2[1];
     }
     return true;
   }
   return false;
}

//------------------------------------------------------------------------------
// Compute the volume of a triangle given 3 node positions
// Note for efficiency this method does not do the necessary division by 2,
// so for the true volume you need to do this!
//------------------------------------------------------------------------------
template<typename RealType>
inline
double
triangleVolume2(const RealType* n1,
		const RealType* n2,
		const RealType* n3) {
  RealType n12[2] = {n2[0] - n1[0], n2[1] - n1[1]};
  RealType n13[2] = {n3[0] - n1[0], n3[1] - n1[1]};
  RealType A[3];
  cross<2, RealType> (n13, n12, A);
  return A[2];
}

//------------------------------------------------------------------------------
// Compute the volume of a tet given the four node positions.
// Note for efficiency this method does not do the necessary division by 6,
// so for the true tet volume you need to do this!
//------------------------------------------------------------------------------
template<typename RealType>
inline
double
tetrahedralVolume6(const RealType* n1,
                   const RealType* n2,
                   const RealType* n3,
                   const RealType* n4) {
  RealType n12[3] = {n2[0] - n1[0], n2[1] - n1[1], n2[2] - n1[2]};
  RealType n13[3] = {n3[0] - n1[0], n3[1] - n1[1], n3[2] - n1[2]};
  RealType n14[3] = {n4[0] - n1[0], n4[1] - n1[1], n4[2] - n1[2]};
  RealType A[3];
  cross<3, RealType> (n13, n12, A);
  return dot<3>(A, n14);
}

//------------------------------------------------------------------------------
// Same thing to simultaneously return the tet volume and centroid.
//------------------------------------------------------------------------------
template<typename RealType>
inline
void
tetrahedralVolumeAndCentroid6(const RealType* n1,
                              const RealType* n2,
                              const RealType* n3,
                              const RealType* n4,
                              RealType& vol,
                              RealType* centroid) {
  RealType n12[3] = {n2[0] - n1[0], n2[1] - n1[1], n2[2] - n1[2]};
  RealType n13[3] = {n3[0] - n1[0], n3[1] - n1[1], n3[2] - n1[2]};
  RealType n14[3] = {n4[0] - n1[0], n4[1] - n1[1], n4[2] - n1[2]};
  RealType A[3];
  cross<3, RealType> (n13, n12, A);
  vol = dot<3>(A, n14);
  centroid[0] = 0.25*(n1[0] + n2[0] + n3[0] + n4[0]);
  centroid[1] = 0.25*(n1[1] + n2[1] + n3[1] + n4[1]);
  centroid[2] = 0.25*(n1[2] + n2[2] + n3[2] + n4[2]);
}

//------------------------------------------------------------------------------
// Compute the volume of a cell in a tessellation.
//------------------------------------------------------------------------------
template<typename RealType>
void
computeCellCentroidAndSignedVolume(const Tessellation<3, RealType>& mesh,
                                   const unsigned ci,
                                   RealType* ccent,
                                   RealType& cvol) {
  POLY_ASSERT(ci < mesh.cells.size());
  ccent[0] = 0.0; ccent[1] = 0.0; ccent[2] = 0.0; 
  cvol = 0.0;
  const unsigned nf = mesh.cells[ci].size();
  int i, j, k, fi, nn, end, delta;
  RealType cmid[3], fcent[3], fhat[3], tetcent[3], tetvol;

  // Make a first stab at the center of the cell based solely on the vertices.
  computeCellCentroid(mesh, ci, cmid);

  // Walk the faces.
  for (k = 0; k != nf; ++k) {
    fi = mesh.cells[ci][k];
    const std::vector<unsigned>& faceNodes = mesh.faces[internal::positiveID(fi)];
    nn = faceNodes.size();
    computeFaceCentroidAndNormal(mesh, internal::positiveID(fi), fcent, fhat);
    if (fi < 0) {
      i = nn - 1;
      end = -1;
      delta = -1;
    } else {
      i = 0; 
      end = nn;
      delta = 1;
    }
    while (i != end) {
      j = (i + delta + nn) % nn;
      POLY_ASSERT(i >= 0 and i < nn);
      POLY_ASSERT(j >= 0 and j < nn);
      tetrahedralVolumeAndCentroid6(cmid, 
                                    &mesh.nodes[3*faceNodes[j]],
                                    &mesh.nodes[3*faceNodes[i]],
                                    fcent,
                                    tetvol,
                                    tetcent);
      cvol += tetvol;
      ccent[0] += tetvol*tetcent[0];
      ccent[1] += tetvol*tetcent[1];
      ccent[2] += tetvol*tetcent[2];
      i += delta;
    }
  }

  // Normalize and we're done.
  POLY_ASSERT(cvol != 0.0);
  ccent[0] /= cvol; ccent[1] /= cvol; ccent[2] /= cvol;
  cvol /= 6.0;
}

//------------------------------------------------------------------------------
// Return a unique set of points.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
inline
void
uniquePoints(const std::vector<RealType>& points,
             const RealType* xmin,
             const RealType* xmax,
             const RealType tol,
             std::vector<RealType>& uniquePointSet,
             std::vector<unsigned>& indexMap) {
  POLY_ASSERT(points.size() % Dimension == 0);
  typedef geometry::Hasher<Dimension, RealType> HasherType;
  // typedef DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  typedef KeyTraits::Key CoordHash;

  // // Compute the bounding box.
  // RealType xmin[Dimension], xmax[Dimension];
  // computeBoundingBox<Dimension, RealType>(points, true, xmin, xmax);

  const unsigned n = points.size()/Dimension;
  std::map<CoordHash, unsigned> uniqueHashes;
  uniquePointSet = std::vector<RealType>();
  indexMap = std::vector<unsigned>();
  unsigned j = 0;
  RealType pos[Dimension];
  for (unsigned i = 0; i != n; ++i) {
    const CoordHash hashi = HasherType::hashPosition(&points[Dimension*i], xmin, xmax, xmin, xmax, tol);
    std::map<CoordHash, unsigned>::const_iterator itr = uniqueHashes.find(hashi);
    if (itr == uniqueHashes.end()) {
      uniqueHashes[hashi] = j;
      HasherType::unhashPosition(pos, xmin, xmax, xmin, xmax, hashi, tol);
      std::copy(&points[Dimension*i], &points[Dimension*(i+1)], std::back_inserter(uniquePointSet));
      indexMap.push_back(j++);
    } else {
      indexMap.push_back(itr->second);
    }
  }
  POLY_ASSERT(uniquePointSet.size() % Dimension == 0);
  POLY_ASSERT(uniquePointSet.size() <= points.size());
  POLY_ASSERT(indexMap.size() == n);
}


//------------------------------------------------------------------------------
// Express a tessellation cell as a ReducedPLC
//------------------------------------------------------------------------------
// Functor to handle the partial dimension specialization
template<int Dimension, typename RealType> struct CellToReducedPLCFunctor;

// 2D
template<typename RealType>
struct CellToReducedPLCFunctor<2, RealType> {
  static ReducedPLC<2, RealType> impl(const Tessellation<2, RealType>& mesh, 
                                      const unsigned icell) {
    POLY_ASSERT(not mesh.empty());
    POLY_ASSERT(icell < mesh.cells.size());
    ReducedPLC<2, RealType> result;
    const unsigned nFaces = mesh.cells[icell].size();
    POLY_ASSERT(nFaces >= 3);
    result.facets.resize(nFaces, std::vector<int>(2));
    for (unsigned i = 0; i != nFaces; ++i) {
      const bool flip = mesh.cells[icell][i] < 0;
      const unsigned iface = flip ? ~mesh.cells[icell][i] : mesh.cells[icell][i];
      POLY_ASSERT(iface < mesh.faces.size());
      POLY_ASSERT(mesh.faces[iface].size() == 2);
      const int ip = flip ? mesh.faces[iface][0] : mesh.faces[iface][1];
      POLY_ASSERT(ip < mesh.nodes.size()/2);
      result.points.push_back(mesh.nodes[2*ip  ]);
      result.points.push_back(mesh.nodes[2*ip+1]);
      result.facets[i][0] = i;
      result.facets[i][1] = (i+1)%nFaces;
    }
    POLY_ASSERT(result.points.size()/2 == nFaces);
    return result;
  }
};

// 3D
template<typename RealType>
struct CellToReducedPLCFunctor<3, RealType> {
  static ReducedPLC<3, RealType> impl(const Tessellation<3, RealType>& mesh, 
                                      const unsigned icell) {
    POLY_ASSERT(icell < mesh.cells.size());
    ReducedPLC<3, double> result;
    std::map<int, int> old2new;
    for (unsigned i = 0; i != mesh.cells[icell].size(); ++i) {
      result.facets.push_back(std::vector<int>());
      if (mesh.cells[icell][i] < 0) {
        const unsigned iface = ~mesh.cells[icell][i];
        const unsigned nnodes = mesh.faces[iface].size();
        for (int j = nnodes - 1; j != -1; --j) {
          const int ip = mesh.faces[iface][j];
          POLY_ASSERT(ip < mesh.nodes.size()/3);
          if (old2new.find(ip) == old2new.end()) {
            old2new[ip] = result.points.size()/3;
            result.points.push_back(mesh.nodes[3*ip  ]);
            result.points.push_back(mesh.nodes[3*ip+1]);
            result.points.push_back(mesh.nodes[3*ip+2]);
          }
          result.facets.back().push_back(old2new[ip]);
        }
        POLY_ASSERT(result.facets.back().size() == nnodes);
      } else {
        const unsigned iface = mesh.cells[icell][i];
        const unsigned nnodes = mesh.faces[iface].size();
        for (int j = 0; j != nnodes; ++j) {
          const int ip = mesh.faces[iface][j];
          POLY_ASSERT(ip < mesh.nodes.size()/3);
          if (old2new.find(ip) == old2new.end()) {
            old2new[ip] = result.points.size()/3;
            result.points.push_back(mesh.nodes[3*ip  ]);
            result.points.push_back(mesh.nodes[3*ip+1]);
            result.points.push_back(mesh.nodes[3*ip+2]);
          }
          result.facets.back().push_back(old2new[ip]);
        }
        POLY_ASSERT(result.facets.back().size() == nnodes);
      }
    }
    POLY_ASSERT(result.facets.size() == mesh.cells[icell].size());
    return result;
  }
};

// Interface
template<int Dimension, typename RealType>
ReducedPLC<Dimension, RealType>
cellToReducedPLC(const Tessellation<Dimension, RealType>& mesh,
                 const unsigned icell) {
   return CellToReducedPLCFunctor<Dimension, RealType>::impl(mesh, icell);
}



//------------------------------------------------------------------------------
// Compare a point to a line and determine if the point is inside the interior
// half-plane (ret -1), inside the exterior half-plane (ret +1), or colinear
// with the line (ret 0). (Method assumes line points (l1,l2) are specified in 
// CCW manner.)
//------------------------------------------------------------------------------
template<typename RealType>
int 
aboveBelow(const RealType& l1x, const RealType& l1y,
           const RealType& l2x, const RealType& l2y,
           const RealType& px,  const RealType& py) {
  const double ztest = (double(l2x - l1x)*double(py - l1y) -
                        double(l2y - l1y)*double(px - l1x));
  return -(ztest < 0.0 ? -1 :
           ztest > 0.0 ?  1 :
                          0);
}

//------------------------------------------------------------------------------
// Compare a point to a plane and determine if it's inside the plane's interior
// half-space (ret -1), inside its exterior half-space (ret +1), or it's
// coplanar with the plane (ret 0). The plane is specified by the point
// normal [(ox, oy, oz), (nx, ny, nz)], and the normal points in the direction
// of the exterior half-space.
//------------------------------------------------------------------------------
template<typename RealType>
int 
aboveBelow(const RealType& ox, const RealType& oy, const RealType& oz,
           const double& nx,   const double& ny,   const double& nz,
           const RealType& px, const RealType& py, const RealType& pz) {
  const double ztest = nx*double(px - ox) + ny*double(py - oy) + nz*double(pz - oz);
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
          0);
}

//------------------------------------------------------------------------------
// Check if an entire cloud of points is inside the interior half-plane of a
// line (ret -1) or is inside the exterior half-plane (ret +1). If there is a
// mixture of interior/exterior/collinear points, return 0.
//------------------------------------------------------------------------------
template<typename RealType>
int 
aboveBelow(const RealType& l1x, const RealType& l1y, 
           const RealType& l2x, const RealType& l2y, 
           const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 2 == 0);
  POLY_ASSERT(points.size() > 1);
  const unsigned n = points.size() / 2;
  const int result = aboveBelow(l1x, l1y, l2x, l2y, points[0], points[1]);
  unsigned i = 1;
  while (i < n and result == aboveBelow(l1x, l1y, 
                                        l2x, l2y, 
                                        points[2*i], points[2*i + 1])) ++i;
  if (i == n) {
    return result;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Check if an entire cloud of points is inside the interior half-space of a
// plane (ret -1) or is inside the exterior half-space (ret +1). If there is a
// mixture of interior/exterior/coplanar points, return 0.
//------------------------------------------------------------------------------
template<typename RealType>
int 
aboveBelow(const RealType& ox, const RealType& oy, const RealType& oz, 
           const RealType& nx,    const RealType& ny, const RealType& nz, 
           const std::vector<RealType>& points) {
  POLY_ASSERT(points.size() % 3 == 0);
  POLY_ASSERT(points.size() > 1);
  const unsigned n = points.size() / 3;
  const int result = aboveBelow(ox, oy, oz, nx, ny, nz, points[0], points[1], points[2]);
  unsigned i = 1;
  while (i < n and result == aboveBelow(ox, oy, oz, 
                                        nx, ny, nz,
                                        points[3*i], points[3*i + 1], points[3*i + 2])) ++i;
  if (i == n) {
    return result;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Compute the 3D normal given three points (a, b, c).
//------------------------------------------------------------------------------
template<typename RealType>
void 
computeNormal(const RealType& ax, const RealType& ay, const RealType& az,
              const RealType& bx, const RealType& by, const RealType& bz,
              const RealType& cx, const RealType& cy, const RealType& cz,
              double& nx, double& ny, double& nz) {
  const double dx_ab = bx - ax, dy_ab = by - ay, dz_ab = bz - az;
  const double dx_ac = cx - ax, dy_ac = cy - ay, dz_ac = cz - az;
  nx = dy_ab*dz_ac - dz_ab*dz_ac;
  ny = dz_ab*dx_ac - dx_ab*dz_ac;
  nz = dx_ab*dy_ac - dy_ab*dx_ac;
}

}
}

#endif
