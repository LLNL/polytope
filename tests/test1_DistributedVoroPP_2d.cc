// Try tessellating a simple lattice of generators in a box in parallel.
// We use randomly chosen seed locations to divide up the generators
// between processors.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include "mpi.h"

#include "polytope.hh"

#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

//------------------------------------------------------------------------------
// Compute the square of the distance.
//------------------------------------------------------------------------------
double distance2(const double x1, const double y1,
                 const double x2, const double y2) {
  return (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Initialize MPI.
  MPI_Init(&argc, &argv);

  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;

  // Figure out our parallel configuration.
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Seed the random number generator the same on all processes.
  srand(10489592);

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = 10; nx != 50; ++nx) {
    if (rank == 0) cout << "Testing nx=" << nx << endl;

    // Create the seed positions for each domain.  Note we rely on this sequence
    // being the same for all processors and therefore don't need to communicate
    // this information.
    vector<double> xproc, yproc;
    xproc.reserve(numProcs);
    yproc.reserve(numProcs);
    for (unsigned iproc = 0; iproc != numProcs; ++iproc) {
      xproc.push_back(x1 + random01()*(x2 - x1));
      yproc.push_back(y1 + random01()*(y2 - y1));
    }

    // Create the local generators.  Note this is not efficient in a couple of ways!  
    // All processes are walking all generators and checking which ones belong to them,
    // and the processor search process is N^2 in the number of processors.  But crimine,
    // this is just supposed to be a little unit test!
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
    unsigned ix, iy;
    double xi, yi;
    for (iy = 0; iy != nx; ++iy) {
      yi = std::max(y1, std::min(y2, y1 + (iy + 0.5)*dy));
      for (ix = 0; ix != nx; ++ix) {
        xi = std::max(x1, std::min(x2, x1 + (ix + 0.5)*dx));
        unsigned owner = 0;
        double minDist2 = distance2(xi, yi, xproc[0], yproc[0]);
        for (unsigned iproc = 1; iproc < numProcs; ++iproc) {
          const double d2 = distance2(xi, yi, xproc[iproc], yproc[iproc]);
          if (d2 < minDist2) {
            owner = iproc;
            minDist2 = d2;
          }
        }
        if (rank == owner) {
          generators.push_back(xi);
          generators.push_back(yi);
        }
      }
    }

    // Create the tessellation.
    double xmin[2] = { x1, y1 };
    double xmax[2] = { x2, y2 };
    polytope::Tessellation<2, double> mesh;
    polytope::VoroPP_2d<double> voro;
    polytope::DistributedTessellator<2, double> distVoro(voro);
    distVoro.tessellate(generators, xmin, xmax, mesh);

    // Gather some global statistics.
    unsigned ncellsGlobal, nnodesGlobal, nfacesGlobal;
    unsigned ncells = mesh.cells.size();
    unsigned nnodes = mesh.nodes.size()/2;
    unsigned nfaces = mesh.faces.size();
    MPI_Allreduce(&ncells, &ncellsGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nnodes, &nnodesGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nfaces, &nfacesGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    // Spew the mesh statistics.
    if (rank == 0) {
      cout << "   num mesh cells : " << ncells << " " << ncellsGlobal << endl;
      cout << "   num mesh nodes : " << nnodes << " " << nnodesGlobal << endl;
      cout << "   num mesh faces : " << nfaces << " " << nfacesGlobal << endl;
    }
//     cout << "Node positions: " << endl;
//     for (unsigned i = 0; i != mesh.nodes.size()/2; ++i) {
//       cout << "   Node " << i << " @ (" << mesh.nodes[2*i] << " " << mesh.nodes[2*i + 1] << ")" << endl;
//     }
//     cout << "Face node sets: " << endl;
//     for (unsigned i = 0; i != nx*nx; ++i) {
//       cout << "   FACES for mesh cell " << i << " :";
//       for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
//       cout << endl;
//     }
//     for (unsigned i = 0; i != mesh.faces.size(); ++i) {
//       double xf = 0.0, yf = 0.0;
//       cout << "   NODES for mesh face " << i << " :";
//       for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
//         unsigned k = mesh.faces[i][j];
//         cout << " " << k;
//         xf += mesh.nodes[2*k];
//         yf += mesh.nodes[2*k + 1];
//       }
//       xf /= mesh.faces[i].size();
//       yf /= mesh.faces[i].size();
//       cout << " @ (" << xf << " " << yf << ")"  << endl;
//     }

//     // Now do the checks.
//     CHECK(mesh.nodes.size()/2 == (nx + 1)*(nx + 1));
//     CHECK(mesh.cells.size() == nx*nx);
//     for (unsigned i = 0; i != nx*nx; ++i) CHECK(mesh.cells[i].size() == 4);
//     CHECK(mesh.faces.size() == 2*nx*(nx + 1));
  }

  cout << "PASS" << endl;
  MPI_Finalize();
  return 0;
}