#include "polytope_c.h"
#include "polytope.hh"
#include "BoostTessellator.hh"
#include "TriangleTessellator.hh"
#include "TetgenTessellator.hh"

using namespace std;
using namespace polytope;

// Helper function for constructing C PLCs from C++ ones.
template <int Dimension>
polytope_plc_t* polytope_plc_from_PLC(const polytope::PLC<Dimension, polytope_real_t>& plc);

namespace 
{

// This helper crafts a C polytope_tessellation object from a C++ Tessellation.
template <int Dimension>
void fill_tessellation(const Tessellation<Dimension, polytope_real_t>& t, polytope_tessellation_t* tess)
{
  // Copy node coordinates.
  tess->num_nodes = t.nodes.size()/Dimension;
  tess->nodes = (polytope_real_t*)malloc(sizeof(polytope_real_t) * Dimension * tess->num_nodes);
  copy(t.nodes.begin(), t.nodes.end(), tess->nodes);

  // Copy cell-face data.
  {
    tess->num_cells = t.cells.size();
    tess->cell_offsets = (int*)malloc(sizeof(int) * tess->num_cells);
    int size = 0;
    for (size_t i = 0; i < t.cells.size(); ++i)
      size += t.cells[i].size();
    tess->cell_faces = (int*)malloc(sizeof(int) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.cells.size(); ++i)
    {
      tess->cell_offsets[i] = offset;
      for (size_t j = 0; j < t.cells[i].size(); ++j, ++offset)
        tess->cell_faces[offset] = t.cells[i][j];
    }
    tess->cell_offsets[t.cells.size()] = offset;
  }

  // Copy face-node data.
  {
    tess->num_faces = t.faces.size();
    tess->face_offsets = (int*)malloc(sizeof(int) * tess->num_faces);
    int size = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
      size += t.faces[i].size();
    tess->face_nodes = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
    {
      tess->face_offsets[i] = offset;
      for (size_t j = 0; j < t.faces[i].size(); ++j, ++offset)
        tess->face_nodes[offset] = t.faces[i][j];
    }
    tess->face_offsets[t.faces.size()] = offset;
  }

  // Inf nodes and faces.
  tess->num_inf_nodes = (int)t.infNodes.size();
  tess->inf_nodes = (unsigned*)malloc(sizeof(unsigned) * tess->num_inf_nodes);
  copy(t.infNodes.begin(), t.infNodes.end(), tess->inf_nodes);

  tess->num_inf_faces = (int)t.infFaces.size();
  tess->inf_faces = (unsigned*)malloc(sizeof(unsigned) * tess->num_inf_faces);
  copy(t.infFaces.begin(), t.infFaces.end(), tess->inf_faces);

  // Cells attached to faces.
  POLY_ASSERT(tess->num_faces == (int)t.faceCells.size());
  tess->face_cells = (int*)malloc(2 * tess->num_faces * sizeof(int));
  for (int f = 0; f < t.faceCells.size(); ++f)
  {
    tess->face_cells[2*f] = t.faceCells[f][0];
    tess->face_cells[2*f+1] = (t.faceCells[f].size() == 2) ? t.faceCells[f][1] : -1;
  }

  // Convex hull.
  if (!t.convexHull.empty())
    tess->convex_hull = polytope_plc_from_PLC(t.convexHull);

  // Neighbor domain information.
  tess->num_neighbor_domains = (int)t.neighborDomains.size();
  tess->neighbor_domains = (unsigned*)malloc(sizeof(unsigned) * tess->num_neighbor_domains);
  copy(t.neighborDomains.begin(), t.neighborDomains.end(), tess->neighbor_domains);

  {
    tess->shared_node_domain_offsets = (int*)malloc(sizeof(int) * (tess->num_neighbor_domains + 1));
    int size = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
      size += t.sharedNodes[i].size();
    tess->shared_nodes = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.sharedNodes.size(); ++i)
    {
      tess->shared_node_domain_offsets[i] = offset;
      for (size_t j = 0; j < t.sharedNodes[i].size(); ++j, ++offset)
        tess->shared_nodes[offset] = t.sharedNodes[i][j];
    }
    tess->shared_node_domain_offsets[t.sharedNodes.size()] = offset;
  }

  {
    tess->shared_face_domain_offsets = (int*)malloc(sizeof(int) * (tess->num_neighbor_domains + 1));
    int size = 0;
    for (size_t i = 0; i < t.faces.size(); ++i)
      size += t.sharedFaces[i].size();
    tess->shared_faces = (unsigned*)malloc(sizeof(unsigned) * (size+1));
    int offset = 0;
    for (size_t i = 0; i < t.sharedFaces.size(); ++i)
    {
      tess->shared_face_domain_offsets[i] = offset;
      for (size_t j = 0; j < t.sharedFaces[i].size(); ++j, ++offset)
        tess->shared_faces[offset] = t.sharedFaces[i][j];
    }
    tess->shared_face_domain_offsets[t.sharedFaces.size()] = offset;
  }


  // Node->cell connectivity.
  // FIXME
}

}

extern "C"
{

//------------------------------------------------------------------------
struct polytope_tessellator_t
{
  Tessellator<2, polytope_real_t>* tess2;
  Tessellator<3, polytope_real_t>* tess3;
};
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_free(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    delete tessellator->tess2;
  else
    delete tessellator->tess3;
  free(tessellator);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_unbounded(polytope_tessellator_t* tessellator,
                                               polytope_real_t* points, int num_points,
                                               polytope_tessellation_t* mesh)
{
  if (tessellator->tess2 != NULL)
  {
    vector<polytope_real_t> pts(points, points + 2*num_points);
    Tessellation<2, polytope_real_t> t;
    tessellator->tess2->tessellate(pts, t);
    fill_tessellation(t, mesh);
  }
  else
  {
    POLY_ASSERT(tessellator->tess3 != NULL);
    vector<polytope_real_t> pts(points, points + 3*num_points);
    Tessellation<3, polytope_real_t> t;
    tessellator->tess3->tessellate(pts, t);
    fill_tessellation(t, mesh);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_box(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* low,
                                            polytope_real_t* high,
                                            polytope_tessellation_t* mesh)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellator_tessellate_in_plc(polytope_tessellator_t* tessellator,
                                            polytope_real_t* points, int num_points,
                                            polytope_real_t* plc_points, int num_plc_points,
                                            polytope_plc_t* piecewise_linear_complex,
                                            polytope_tessellation_t* mesh)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
bool polytope_tessellator_handles_plcs(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->handlesPLCs();
  else
    return tessellator->tess3->handlesPLCs();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
const char* polytope_tessellator_name(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->name().c_str();
  else
    return tessellator->tess3->name().c_str();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_real_t polytope_tessellator_degeneracy(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return tessellator->tess2->degeneracy();
  else
    return tessellator->tess3->degeneracy();
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int polytope_tessellator_dimension(polytope_tessellator_t* tessellator)
{
  if (tessellator->tess2 != NULL)
    return 2;
  else
    return 3;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#if HAVE_BOOST
polytope_tessellator_t* boost_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = new BoostTessellator<polytope_real_t>();
  t->tess3 = NULL;
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#if HAVE_TRIANGLE
polytope_tessellator_t* triangle_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = new TriangleTessellator<polytope_real_t>();
  t->tess3 = NULL;
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
#if HAVE_TETGEN
polytope_tessellator_t* tetgen_tessellator_new()
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  t->tess2 = NULL;
  t->tess3 = new TetgenTessellator();
  return t;
}
#endif
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_tessellator_t* voroplusplus_tessellator_new(int dimension)
{
  polytope_tessellator_t* t = (polytope_tessellator_t*)malloc(sizeof(polytope_tessellator_t));
  if (dimension == 2)
  {
    t->tess2 = new VoroPP_2d<polytope_real_t>();
    t->tess3 = NULL;
  }
  else
  {
    t->tess2 = NULL;
    t->tess3 = new VoroPP_3d<polytope_real_t>();
  }
  return t;
}
//------------------------------------------------------------------------

}

