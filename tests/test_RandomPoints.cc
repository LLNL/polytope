// test_RandomPoints
//
// Stress test for meshing complicated PLC boundaries with/without holes.
// Iterate over each of the default boundaries defined in Boundary2D.hh
// and tessellate using N randomly-distributed generators for N=10,100,1000.
// Can test both Triangle and Voro++ 2D tessellators. Voro++ has been
// commented out since it currently lacks PLC capabilities.

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D<double>& boundary,
                  Tessellator<2,double>& tessellator,
                  int boundaryID) {
  // output name
  ostringstream os;
  os << "RandomPoints_" << tessellator.name();
  string testName = os.str();

  Generators<2,double> generators(boundary);
  unsigned nPoints = 10;
  Tessellation<2,double> mesh;
  for( unsigned n = 0; n < 3; ++n ){
    POLY_ASSERT(mesh.empty());
    nPoints = nPoints * 10;
    int plotIndex = 3*boundaryID + n;

    cout << nPoints << " points..." << endl;
    generators.randomPoints( nPoints );      
    tessellate2D(generators.mPoints,boundary,tessellator,mesh);    
    outputMesh(mesh, testName, generators.mPoints, plotIndex);
    mesh.clear();   
    plotIndex++;
  }
}

// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
void testAllBoundaries(Tessellator<2,double>& tessellator) {
  for (int bid = 0; bid < 10; ++bid){
    cout << "Testing boundary type " << bid << endl;
    Boundary2D<double> boundary;
    boundary.setDefaultBoundary(bid);
    testBoundary(boundary, tessellator, bid);
  }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef HAVE_MPI
   MPI_Init(&argc, &argv);
#endif


#ifdef HAVE_TRIANGLE
   {
     cout << "\nTriangle Tessellator:\n" << endl;
     TriangleTessellator<double> tessellator;
     testAllBoundaries(tessellator);
   }
#endif
   
#ifdef HAVE_BOOST_VORONOI
   {
     cout << "\nBoost Tessellator:\n" << endl;
     BoostTessellator<double> tessellator;
     testAllBoundaries(tessellator);
   }
#endif      

   cout << "PASS" << endl;

#ifdef HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
