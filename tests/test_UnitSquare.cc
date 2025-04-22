// test_UnitSquare
//
// Meshes a unit square with (nx-by-nx) Cartesian generators with nx increasing
// by a given factor between iterations. Performs checks on the resulting
// tessellation to see if it is indeed Cartesian.
// Triangle and Boost tessellators are tested here.
// -----------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"

#include "timingUtilities.hh"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// checkCartesianMesh
// -----------------------------------------------------------------------
void checkCartesianMesh(Tessellation<2,double>& mesh, unsigned nx, unsigned ny)
{
   POLY_CHECK(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
   POLY_CHECK(mesh.cells.size()   == nx*ny );
   POLY_CHECK(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
   for (unsigned i = 0; i != nx*ny; ++i) POLY_CHECK(mesh.cells[i].size() == 4);

   std::vector<std::set<unsigned> > nodeCells = mesh.computeNodeCells();
   for (unsigned i = 0; i != (nx+1)*(ny+1); ++i)
   {
      POLY_CHECK( (nodeCells[i].size() == 4) ||
                  (nodeCells[i].size() == 2) ||
                  (nodeCells[i].size() == 1) );
   }
}

// -----------------------------------------------------------------------
// generateMesh
// -----------------------------------------------------------------------
void generateMeshes(Tessellator<2,double>& tessellator,
                    int Nmin, int Nmax, int increment)
{
  // Set the boundary
  Boundary2D<double> boundary;
  boundary.setUnitSquare();
  Generators<2,double> generators( boundary );

  for (unsigned nx = Nmin; nx <= Nmax; nx += increment) {
    cout << "Testing nx=" << nx << endl;

    // Create generators
    std::vector<unsigned> nxny(2,nx);
    generators.cartesianPoints(nxny);
    Tessellation<2,double> mesh;

    Timing::Time start = Timing::currentTime();
    tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
    double time = Timing::difference(start, Timing::currentTime());

    cout << "Mesh construction took " << time << " seconds." << endl;

    // CHECKS:
    cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
    cout << "   num mesh cells : " << mesh.cells.size()   << endl;
    cout << "   num mesh faces : " << mesh.faces.size()   << endl;
    checkCartesianMesh(mesh,nx,nx);

    // ofstream os;
    // os.open("oneMillionRandomGenerators.txt");
    // os << "# x \t\t y \t\t z" << endl;
    // for (unsigned j = 0; j != generators.mPoints.size()/2; ++j) {
    //   os << generators.mPoints[2*j] << " \t "
    //      << generators.mPoints[2*j+1] << " \t " << 0.0 << endl;
    // }
    // os.close();
  }
}

int positiveIntFromString(const char* s)
{
  int i = atoi(s);
  if (i < 1)
  {
    cout << "Invalid argument: " << i << " (must be positive)." << endl;
    cout << "FAIL" << endl;
    exit(-1);
  }
  return i;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Accepts 3 input parameters:
  // 1. The minimum value of nx (defaults to 2).
  // 2. The maximum value of nx (defaults to 100).
  // 3. The positive (integer) increment by which nx is increased after each
  //    mesh generation (defaults to 1).
  int increment = 1;
  int Nmin = 2;
  int Nmax = 100;

  if (argc >= 2)
    Nmin = positiveIntFromString(argv[1]);
  if (argc >= 3)
  {
    Nmax = positiveIntFromString(argv[2]);
    if (Nmax < Nmin)
    {
      cout << "Invalid Nmax: " << Nmax << " (must be >= Nmin)." << endl;
      cout << "FAIL" << endl;
      exit(-1);
    }
  }
  if (argc >= 4)
  {
    increment = positiveIntFromString(argv[3]);
    if (increment > (Nmax - Nmin))
    {
      cout << "Invalid increment: " << increment << " (must be less than"
           << Nmax - Nmin << ")." << endl;
      cout << "FAIL" << endl;
      exit(-1);
    }
  }

#ifdef HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    generateMeshes(tessellator, Nmin, Nmax, increment);
  }
#endif

#ifdef HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    generateMeshes(tessellator, Nmin, Nmax, increment);
  }
#endif

   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // generateMesh(voro);
   // cout << "Voro 2D: PASS" << endl;

  cout << "PASS" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
