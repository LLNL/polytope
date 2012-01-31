//---------------------------------Spheral++----------------------------------//
// VoroPP_3d
// 
// An implemenation of the Tessellator interface that uses the 3D Voro++
// library.
//----------------------------------------------------------------------------//
#ifndef __Polytope_VoroPP_3d__
#define __Polytope_VoroPP_3d__

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope {

template<typename Real>
class VoroPP_3d: public Tessellator<3, Real> {

  //--------------------------- Public Interface ---------------------------//
public:

  //! Constructor.
  //! The parameters (nx, ny) are used internally to Voro++ in order to make
  //! the selection of generators that can influence any particular generator
  //! more efficient.  The results of the tessellation should be independent of
  //! of these choices -- they only affect computational expense.
  //! \param nx The number of boxes to carve the volume into in the x direction.
  //! \param ny The number of boxes to carve the volume into in the y direction.
  //! \param degeneracy The tolerance for merging nodes in a cell.
  VoroPP_3d(const unsigned nx = 20,
            const unsigned ny = 20,
            const unsigned nz = 20,
            const Real degeneracy = 1.0e-14);
  ~VoroPP_3d();

  //! Generate a Voronoi tessellation for the given set of generator points
  //! with a bounding box specified by \a low and \a high. Here, low[i]
  //! contains the ith coordinate for the "lower-left-near" corner of the 
  //! bounding box in 3D or 3D, and high[i] contains the corresponding 
  //! opposite corner. The coordinates of these points are stored in 
  //! point-major order and the 0th component of the ith point appears in 
  //! points[Dimension*i].
  //! \param points A (Dimension*numPoints) array containing point coordinates.
  //! \param low The coordinates of the "lower-left-near" bounding box corner.
  //! \param high The coordinates of the "upper-right-far" bounding box corner.
  //! \param mesh This will store the resulting tessellation.
  virtual void tessellate(std::vector<Real>& points,
                          Real* low,
                          Real* high,
                          Tessellation<3, Real>& mesh) const;


  // This Tessellator does not handle PLCs... yet.
  bool handlesPLCs() const { return false; }

  // Access our attributes.
  unsigned nx() const { return mNx; }
  unsigned ny() const { return mNy; }
  unsigned nz() const { return mNz; }
  Real degeneracy() const { return std::sqrt(mDegeneracy2); }

private:
  unsigned mNx, mNy, mNz;
  Real mDegeneracy2;

  // Forbidden methods.
  VoroPP_3d(const VoroPP_3d&);
  VoroPP_3d& operator=(const VoroPP_3d&);
};

}

#endif
