
#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"
#include "polytope_geometric_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------------
// Compute the square of the distance.
//------------------------------------------------------------------------------
double distance2(const double x1, const double y1,
                 const double x2, const double y2) {
  return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
}

// -----------------------------------------------------------------------
// lloyd
// -----------------------------------------------------------------------
void lloyd(Tessellation<2,double>& mesh,
           vector<double>& points) {
   for (unsigned i = 0; i < mesh.cells.size(); ++i){
      double cent[2], area;
      geometry::computeCellCentroidAndSignedArea(mesh, i, 1.0e-12, cent, area);
      points[2*i  ] = 0.5*(points[2*i  ] + cent[0]);
      points[2*i+1] = 0.5*(points[2*i+1] + cent[1]);
   }
}

// -----------------------------------------------------------------------
// lloydTestDistributed
// -----------------------------------------------------------------------
void lloydTestDistributed(Tessellator<2,double>& tessellator) {
  const unsigned nPoints     = 2000;     // Number of generators
  const unsigned nIter       = 100;     // Number of iterations
  const unsigned outputEvery = 5;        // Output frequency
  const int btype = 9;

  // Seed the random number generator the same on all processes.
  srand(10489592);
  
  // Test name
  string testName = "Distributed_LloydTest_" + tessellator.name();

  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2,double> generators(boundary);
  generators.randomPoints(nPoints);

  // Assign points to processors in quasi-Voronoi fashion
  vector<double> points, xproc, yproc;
  double p[2];
  xproc.reserve(numProcs);
  yproc.reserve(numProcs);
  for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
    boundary.getPointInside(p);
    xproc.push_back(p[0]);
    yproc.push_back(p[1]);
  }
  for (unsigned i = 0; i < nPoints; ++i){
    unsigned owner = 0;
    double minDist2 = distance2(generators.mPoints[2*i], 
                                generators.mPoints[2*i+1], 
                                xproc[0], yproc[0]);
    for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
      const double d2 = distance2(generators.mPoints[2*i], 
                                  generators.mPoints[2*i+1], 
                                  xproc[iproc], yproc[iproc]);
      if( d2 < minDist2 ){
        owner = iproc;
        minDist2 = d2;
      }
    }
    if (rank == owner) {
      points.push_back(generators.mPoints[2*i  ]);
      points.push_back(generators.mPoints[2*i+1]);
    }
  }
  
  // Initialize mesh and tessellator
  Tessellation<2,double> mesh;
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);

  // Do the Lloyd iteration thang
  unsigned iter = 0;
  outputMesh(mesh, testName, points, iter);
  while (iter != nIter) {
    lloyd(mesh,points);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    if (iter % outputEvery == 0) 
       outputMesh(mesh, testName, points, iter);
  }
}


// -----------------------------------------------------------------------
// lloydTest
// -----------------------------------------------------------------------
void lloydTest(Tessellator<2,double>& tessellator) {
  const unsigned nPoints = 1000;     // Number of generators
  const unsigned nIter   = 100;      // Number of iterations
  const int btype = 2;
  
  string testName = "Centroidal_LloydTest_" + tessellator.name();

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);
  Generators<2,double> generators(boundary);
  generators.randomPoints(nPoints);
  std::vector<double> points;
  for (unsigned i = 0; i != nPoints; ++i) {
    if (boundary.testInside(&generators.mPoints[2*i])) {
      std::copy(&generators.mPoints[2*i], &generators.mPoints[2*i+2], std::back_inserter(points));
    }
  }
   
  // Initialize mesh and tessellator
  Tessellation<2,double> mesh;
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);

  unsigned iter = 0;
  outputMesh(mesh, testName, points, iter);
  while (iter != nIter) {
    lloyd(mesh,points);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, points, iter);
  }
}

// -----------------------------------------------------------------------
// cleaningTest
// -----------------------------------------------------------------------
void cleaningTest(Tessellator<2,double>& tessellator) {
  const unsigned nPoints = 100;     // Number of generators
  const unsigned nIter   = 100;     // Number of iterations
  const double edgeTol = 0.001;     // Relative small-edge tolerance
  const int btype = 3;

  string testName = "Centroidal_CleaningTest_" + tessellator.name();

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);  
  Generators<2,double> generators(boundary);
  generators.randomPoints(nPoints);
  std::vector<double> points;
  for (unsigned i = 0; i != nPoints; ++i) {
    if (boundary.testInside(&generators.mPoints[2*i])) {
      std::copy(&generators.mPoints[2*i], &generators.mPoints[2*i+2], std::back_inserter(points));
    }
  }
   
  // Initialize mesh and tessellator
  Tessellation<2,double> mesh;
  MeshEditor<2, double> meshEditor(mesh);
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
  
  unsigned iter = 0;
  outputMesh(mesh, testName, points, iter);
  while (iter != nIter) {
    meshEditor.cleanEdges(edgeTol);
    lloyd(mesh,points);
    // meshEditor.cleanEdges(edgeTol);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, points, iter);
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif
   
#if HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    lloydTest(tessellator);
    cleaningTest(tessellator);
  }
#endif   

#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    lloydTest(tessellator);
    cleaningTest(tessellator);
  }
#endif

#if HAVE_MPI
  {
    cout << "\nDistributed Triangle:\n" << endl;
    DistributedTessellator<2,double> tessellator(new TriangleTessellator<double>(), true, true);
    lloydTestDistributed(tessellator);
  }
#endif

  cout << "PASS" << endl;
   
#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
