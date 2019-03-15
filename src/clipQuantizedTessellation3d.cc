//------------------------------------------------------------------------------
// 3D implementation of clipQuantizedTessellation.
//
// This method relies upon the Boost.Polygon library for geometric operations.
//------------------------------------------------------------------------------

#include "clipQuantizedTessellation.hh"
#include "polytope_internal.hh"

namespace polytope {

template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation3d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<3, RealType>& geometry,
                               const Tessellator<3, RealType>& tessellator) {
  POLY_ASSERT2(false, "Implement me!");
}

//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
template
void
clipQuantizedTessellation<int, double>(QuantizedTessellation3d<int, double>& qmesh,
                                       const std::vector<double>& PLCpoints,
                                       const PLC<3, double>& geometry,
                                       const Tessellator<3, double>& tessellator);
}
