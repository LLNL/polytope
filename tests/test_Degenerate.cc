// test_Degenerate
//
// Tessellate a unit square with 50x50 equally-spaced Cartesian generators,
// then randomly perturb their positions by +/- epsilon for epsilon increasing
// from 1e-16 to 1e-7 by factors of 10. Check to see if the resulting
// tessellation is Cartesian. If not, compute the minimum face length.
// Test applied to both Triangle and Voro++ 2D tessellators.

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

#if HAVE_MPI
#include "mpi.h"
#endif

#define POLY_CHECK_BOOL(x) if (!(x)) { return false; }

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// minLength
// -----------------------------------------------------------------------
double minLength(Tessellation<2,double>& mesh)
{
   double faceLength = FLT_MAX;
   for (unsigned iface = 0; iface != mesh.faces.size(); ++iface)
   {
      POLY_ASSERT( mesh.faces[iface].size() == 2 );
      const unsigned inode0 = mesh.faces[iface][0];
      const unsigned inode1 = mesh.faces[iface][1];
      double x0 = mesh.nodes[2*inode0], y0 = mesh.nodes[2*inode0+1];
      double x1 = mesh.nodes[2*inode1], y1 = mesh.nodes[2*inode1+1];
      double len = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
      faceLength = min( faceLength, sqrt(len) );
   }
   return faceLength;
}

// -----------------------------------------------------------------------
// checkIfCartesian
// -----------------------------------------------------------------------
bool checkIfCartesian(Tessellation<2,double>& mesh, unsigned nx, unsigned ny)
{
   POLY_CHECK_BOOL(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
   POLY_CHECK_BOOL(mesh.cells.size()   == nx*ny );
   POLY_CHECK_BOOL(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
   for (unsigned i = 0; i != nx*ny; ++i) POLY_CHECK_BOOL(mesh.cells[i].size() == 4);
   
   std::vector<std::set<unsigned> > nodeCells = mesh.computeNodeCells();
   for (unsigned i = 0; i != (nx+1)*(ny+1); ++i)
   {
      POLY_CHECK_BOOL( (nodeCells[i].size() == 4) ||
                       (nodeCells[i].size() == 2) ||
                       (nodeCells[i].size() == 1) );
   }
   return true;
}

// -----------------------------------------------------------------------
// generateMesh
// -----------------------------------------------------------------------
void generateMesh(Tessellator<2,double>& tessellator)
{
   // Set the boundary
   Boundary2D<double> boundary;
   boundary.setUnitSquare();
   Generators<2,double> generators(boundary);

   unsigned nx = 10;
   std::vector<unsigned> nxny(2,nx);
   
   // Create generators
   cout << "Generator locations randomly perturbed by" << endl;
   double epsilon = 2.0e-12;
   for ( int i = 0; i != 8; ++i, epsilon *= 10)
   {
      cout << "+/- " << epsilon/2 << "...";
      generators.cartesianPoints( nxny );         // reset locations
      generators.perturb( epsilon );              // perturb
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints,boundary,tessellator,mesh);
      bool isCartesian = checkIfCartesian(mesh,nx,nx);
      if( isCartesian ){ 
         cout << "PASS" << endl; 
      }else{
         cout << "Degeneracy reached! Minimum face length = " << minLength(mesh) << endl;
      }

#if HAVE_SILO
      vector<double> index( mesh.cells.size());
      for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
      map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
      cellFields["cell_index"] = &index[0];
      ostringstream os;
      os << "test_Degenerate_" << epsilon;
      polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                             faceFields, cellFields, os.str());
#endif
   }
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

#if HAVE_TRIANGLE
   cout << "\nTriangle Tessellator:\n" << endl;
   TriangleTessellator<double> triangle;
   generateMesh(triangle);
#endif   

   // NOTE: Voro does not give sensible results at this time
   //
   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // generateMesh(voro);

   cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
   return 0;
}
