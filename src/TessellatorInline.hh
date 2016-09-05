#include "clipQuantizedTessellation.hh"
#include "makeBoxPLC.hh"

namespace polytope {

//----------------------------------------------------------------------------
// Tessellate (unbounded)
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           Tessellation<nDim, RealType>& mesh) const {

  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantizedTessellation quantmesh(points);
  this->tessellateQuantized(quantmesh);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);
}

//----------------------------------------------------------------------------
// Tessellate in a rectangle.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<nDim, RealType>& mesh) const {
  ReducedPLC<nDim, RealType> geometry;
  internal::makeBoxPLC(geometry, low, high);
  this->tessellate(points, geometry, mesh);
}

//----------------------------------------------------------------------------
// Tessellate in a PLC.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           const std::vector<RealType>& PLCpoints,
           const PLC<nDim, RealType>& geometry,
           Tessellation<nDim, RealType>& mesh) const {
  
  // Pre-conditions
  POLY_ASSERT(mesh.empty());
  POLY_ASSERT(points.size() > 0);
  POLY_ASSERT(points.size() % 2 == 0);

  // Invoke the descendant method to fill the quant mesh.
  QuantizedTessellation quantmesh(points);
  this->tessellateQuantized(quantmesh);

  // Clip against the boundary.
  clipQuantizedTessellation(quantmesh, points, geometry);

  // Copy the QuantTessellation to the output.
  quantmesh.fillTessellation(mesh);
}

//----------------------------------------------------------------------------
// Tessellate in a ReducedPLC.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
void
Tessellator<nDim, RealType>::
tessellate(const std::vector<RealType>& points,
           const ReducedPLC<nDim, RealType>& geometry,
           Tessellation<nDim, RealType>& mesh) const {
  this->tessellate(points, geometry.points, geometry, mesh);
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Unbounded case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[nDim], xmax[nDim];
  geometry::computeBoundingBox<nDim, RealType>(points, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// Bounded by a box.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     RealType* low,
                     RealType* high,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, low, high, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, low, high, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// PLC case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const std::vector<RealType>& PLCpoints,
                     const PLC<nDim, RealType>& geometry,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {

  // Compute the bounding box.
  RealType xmin[nDim], xmax[nDim];
  geometry::computeBoundingBox<nDim, RealType>(PLCpoints, true, xmin, xmax);

  // Hash the input generators to a unique set.
  std::vector<RealType> uniquePoints;
  std::vector<unsigned> result;
  geometry::uniquePoints<nDim, RealType>(points, xmin, xmax, tol, uniquePoints, result);

  // Tessellate the unique set, and return the mapping.
  this->tessellate(uniquePoints, PLCpoints, geometry, mesh);
  return result;
}

//------------------------------------------------------------------------------
// Do a tessellation with potential degenerate generators.
// ReducedPLC case.
//------------------------------------------------------------------------------
template<int nDim, typename RealType>
inline
std::vector<unsigned>
Tessellator<nDim, RealType>::
tessellateDegenerate(const std::vector<RealType>& points,
                     const ReducedPLC<nDim, RealType>& geometry,
                     const RealType tol,
                     Tessellation<nDim, RealType>& mesh) const {
  return this->tessellateDegenerate(points, geometry.points, geometry, tol, mesh);
}

}
