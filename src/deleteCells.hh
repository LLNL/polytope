//----------------------------------------------------------------------------//
// deleteCells
// A method which deletes a set of cells from a Tessellation, along with
// any nodes and faces which no longer are part of a cell.  The remaining
// elements are renumbered accordingly.
// The inputs are:
// 1.  mesh : the mesh to be edited.
// 2.  cellMask : an array of length mesh.cells.size of either 0 or 1:
//     0 => delete cell
//     1 => keep cell
//----------------------------------------------------------------------------//
#ifndef __polytope_deleteCells__
#define __polytope_deleteCells__
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>

#include "polytope.hh"

namespace polytope {

template<int Dimension, typename RealType>
void
deleteCells(Tessellation<Dimension, RealType>& mesh,
            const std::vector<unsigned>& cellMask) {

  // Pre-conditions.
  const unsigned ncells0 = mesh.cells.size();
  const unsigned nfaces0 = mesh.faces.size();
  const unsigned nnodes0 = mesh.nodes.size()/Dimension;
  ASSERT(cellMask.size() == ncells0);
  ASSERT(ncells0 == 0 or *max_element(cellMask.begin(), cellMask.end()) == 1);

  // Create masks for the nodes and faces.
  std::vector<int> nodeMask(nnodes0, 0), faceMask(nfaces0, 0);
  for (unsigned icell = 0; icell != ncells0; ++icell) {
    if (cellMask[icell] == 1) {
      for (std::vector<int>::const_iterator faceItr = mesh.cells[icell].begin();
           faceItr != mesh.cells[icell].end();
           ++faceItr) {
        const unsigned iface = (*faceItr >= 0 ? *faceItr : ~(*faceItr));
        ASSERT(iface < mesh.faces.size());
        faceMask[iface] = 1;
        for (std::vector<unsigned>::const_iterator nodeItr = mesh.faces[iface].begin();
             nodeItr != mesh.faces[iface].end();
             ++nodeItr) {
          const unsigned inode = *nodeItr;
          ASSERT(inode < mesh.nodes.size()/Dimension);
          nodeMask[inode] = 1;
        }
      }
    }
  }

  // Determine the new cell, face, and node numberings.
  unsigned nnodes1 = 0, nfaces1 = 0, ncells1 = 0;
  std::map<unsigned, unsigned> old2new_nodes, old2new_faces, old2new_cells;
  {
    for (unsigned i = 0; i != nnodes0; ++i) {
      if (nodeMask[i] == 1) old2new_nodes[i] = nnodes1++;
    }
    for (unsigned i = 0; i != nfaces0; ++i) {
      if (faceMask[i] == 1) old2new_faces[i] = nfaces1++;
    }
    for (unsigned i = 0; i != ncells0; ++i) {
      if (cellMask[i] == 1) old2new_cells[i] = ncells1++;
    }
  }

  // Reconstruct the nodes.
  {
    std::vector<RealType> newNodes;
    newNodes.reserve(nnodes1);
    for (unsigned i = 0; i != nnodes0; ++i) {
      if (nodeMask[i] == 1) {
        std::copy(&mesh.nodes[Dimension*i], &mesh.nodes[Dimension*(i + 1)], back_inserter(newNodes));
      }
    }
    mesh.nodes = newNodes;
  }

  // Reconstruct the faces.
  {
    std::vector<std::vector<unsigned> > newFaces, newFaceCells;
    newFaces.reserve(nfaces1);
    newFaceCells.reserve(nfaces1);
    for (unsigned i = 0; i != nfaces0; ++i) {
      if (faceMask[i] == 1) {
        newFaces.push_back(mesh.faces[i]);
        for (std::vector<unsigned>::iterator itr = newFaces.back().begin();
             itr != newFaces.back().end();
             ++itr) {
          ASSERT(old2new_nodes.find(*itr) != old2new_nodes.end());
          *itr = old2new_nodes[*itr];
        }
        ASSERT(mesh.faceCells[i].size() == 1 or
               mesh.faceCells[i].size() == 2);
        newFaceCells.push_back(std::vector<unsigned>());
        if (cellMask[mesh.faceCells[i][0]] == 1) newFaceCells.back().push_back(old2new_cells[mesh.faceCells[i][0]]);
        if (mesh.faceCells[i].size() == 2 and
            cellMask[mesh.faceCells[i][1]] == 1) newFaceCells.back().push_back(old2new_cells[mesh.faceCells[i][1]]);
        ASSERT(newFaceCells.back().size() == 1 or
               newFaceCells.back().size() == 2);
      }
    }
    mesh.faces = newFaces;
    mesh.faceCells = newFaceCells;
  }

  // Reconstruct the cells.
  {
    std::vector<std::vector<int> > newCells;
    newCells.reserve(ncells1);
    for (unsigned i = 0; i != ncells0; ++i) {
      if (cellMask[i] == 1) {
        newCells.push_back(mesh.cells[i]);
        for (std::vector<int>::iterator itr = newCells.back().begin();
             itr != newCells.back().end();
             ++itr) {
          const int iface = (*itr >= 0 ? *itr : ~(*itr));
          ASSERT(old2new_faces.find(iface) != old2new_faces.end());
          *itr = (*itr >= 0 ? old2new_faces[iface] : ~old2new_faces[iface]);
        }
      }
    }
    mesh.cells = newCells;
  }

  // If there was a convex hull in the mesh, it's probably no longer valid.
  mesh.convexHull = PLC<Dimension, RealType>();

  // Post-conditions.
  ASSERT(mesh.nodes.size() == Dimension*nnodes1);
  ASSERT(mesh.faces.size() == nfaces1);
  ASSERT(mesh.faceCells.size() == nfaces1);
  ASSERT(mesh.cells.size() == ncells1);
}

}

#endif
