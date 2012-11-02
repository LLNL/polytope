//------------------------------------------------------------------------
// TriangleTessellator
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include <map>
#include <limits>
#include "float.h"

#include "polytope.hh" // Pulls in ASSERT and TriangleTessellator.hh.
#include "convexHull_2d.hh"

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>

// Since triangle isn't built to work out-of-the-box with C++, we 
// slurp in its source here, bracketing it with the necessary dressing.
#define TRILIBRARY
#define REAL double
#define ANSI_DECLARATORS 
#define CDT_ONLY // Conforming Delaunay triangulations only! 

// Because Hang Si has messed with some of Shewchuk's predicates and
// included them with his own Tetgen library, we need to rename some of 
// the symbols therein to prevent duplicate symbols from confusing the 
// linker. Gross.
#define exactinit triangle_exactinit
#define fast_expansion_sum_zeroelim triangle_fast_expansion_sum_zeroelim
#define scale_expansion_zeroelim triangle_scale_expansion_zeroelim
#define estimate triangle_estimate
#define orient3dadapt triangle_orient3dadapt
#define orient3d triangle_orient3d
#define incircleadapt triangle_incircleadapt
#define incircle triangle_incircle
#include "triangle.c"

//------------------------------------------------------------------------------
// Teach Boost.Geometry how to handle our Point2 class with appropriate traits.
//------------------------------------------------------------------------------
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y)

// Fast predicate for determining colinearity of points.
extern double orient2d(double* pa, double* pb, double* pc);

namespace polytope {

using namespace std;
using std::min;
using std::max;
using std::abs;

namespace {

//------------------------------------------------------------------------
// This function computes the circumcenter of a triangle with vertices
// A = (Ax, Ay), B = (Bx, By), and C = (Cx, Cy), and places the result 
// in X.
//------------------------------------------------------------------------
void 
computeCircumcenter(double* A, double* B, double* C, double* X)
{
  // This solution was taken from Wikipedia's entry:
  // http://en.wikipedia.org/wiki/Circumscribed_circle
  double D = 2.0*(A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]));
  X[0] = ((A[0]*A[0] + A[1]*A[1])*(B[1]-C[1]) + (B[0]*B[0] + B[1]*B[1])*(C[1]-A[1]) + 
          (C[0]*C[0] + C[1]*C[1])*(A[1]-B[1]))/D;
  X[1] = ((A[0]*A[0] + A[1]*A[1])*(C[0]-B[0]) + (B[0]*B[0] + B[1]*B[1])*(A[0]-C[0]) + 
          (C[0]*C[0] + C[1]*C[1])*(B[0]-A[0]))/D;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned> {
public:
  CounterMap(): std::map<Key, unsigned>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned>::operator[](key) = 0U;
      itr = this->find(key);
    }
    ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// An implementation of the map specialized for testing true/false.
// If tested with a key that does not exist, it is initialized as false.
//------------------------------------------------------------------------------
template<typename Key, 
         typename Comparator = std::less<Key> >
class BoolMap: public std::map<Key, bool> {
public:
  BoolMap(): std::map<Key, bool>() {}
  virtual ~BoolMap() {}
  bool operator[](const Key& key) const {
    typename std::map<Key, bool>::const_iterator itr = this->find(key);
    if (itr == this->end()) return false;
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// Hash two node indices uniquely to represent an edge.
//------------------------------------------------------------------------------
std::pair<int, int>
hashEdge(const int i, const int j) {
  ASSERT(i != j);
  return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

//------------------------------------------------------------------------------
// Predicate to check if either element of a std::pair matches the given value.
//------------------------------------------------------------------------------
template<typename T>
struct MatchEitherPairValue {
  T mval;
  MatchEitherPairValue(const T& x): mval(x) {}
  bool operator()(const std::pair<T, T>& x) const { return (x.first == mval or x.second == mval); }
};

//------------------------------------------------------------------------------
// Update the map of thingies to unique indices.
//------------------------------------------------------------------------------
template<typename Key>
int
addKeyToMap(const Key& key, std::map<Key, int>& key2id) {
  const typename map<Key, int>::const_iterator itr = key2id.find(key);
  int result;
  if (itr == key2id.end()) {
    result = key2id.size();
    key2id[key] = result;
  } else {
    result = itr->second;
  }
  return result;
}

} // end anonymous namespace

//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
TriangleTessellator():
  Tessellator<2, RealType>() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
TriangleTessellator<RealType>::
~TriangleTessellator() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           Tessellation<2, RealType>& mesh) const {
  // Use the PLC method with an empty geometry.
  PLC<2, RealType> geometry;
  tessellate(points, points, geometry, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           RealType* low,
           RealType* high,
           Tessellation<2, RealType>& mesh) const 
{
  // Build a PLC with the bounding box, and then use the PLC method.
  ReducedPLC<2, RealType> box = this->boundingBox(low, high);
  this->tessellate(points, box.points, box, mesh);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RealType>
void
TriangleTessellator<RealType>::
tessellate(const vector<RealType>& points,
           const vector<RealType>& PLCpoints,
           const PLC<2, RealType>& geometry,
           Tessellation<2, RealType>& mesh) const 
{
  ASSERT(!points.empty());

  // Make sure we're not modifying an existing tessellation.
  ASSERT(mesh.empty());

  typedef std::pair<int, int> EdgeHash;
  typedef uint64_t CoordHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<double> RealPoint;
  typedef boost::geometry::model::polygon<RealPoint, // point type
                                          false>     // clockwise
    BGpolygon;
  typedef boost::geometry::model::ring<RealPoint,    // point type
                                       false>        // clockwise
    BGring;
  typedef boost::geometry::model::multi_point<RealPoint> BGmulti_point;
  const CoordHash coordMax = numeric_limits<CoordHash>::max() / 2U;

  triangulateio in, delaunay;

  // Find the range of the generator points.
  const unsigned numGenerators = points.size()/2;
  const unsigned numPLCpoints = PLCpoints.size()/2;
  RealType  low[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType high[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  int i, j;
  for (i = 0; i != numGenerators; ++i) {
    low[0] = min(low[0], points[2*i]);
    low[1] = min(low[1], points[2*i+1]);
    high[0] = max(high[0], points[2*i]);
    high[1] = max(high[1], points[2*i+1]);
  }
  for (i = 0; i != numPLCpoints; ++i) {
    low[0] = min(low[0], PLCpoints[2*i]);
    low[1] = min(low[1], PLCpoints[2*i+1]);
    high[0] = max(high[0], PLCpoints[2*i]);
    high[1] = max(high[1], PLCpoints[2*i+1]);
  }
  ASSERT(low[0] < high[0] and low[1] < high[1]);
  RealType box[2] = {high[0] - low[0], 
                     high[1] - low[1]};
  const double boxsize = 2.0*max(box[0], box[1]);
  const double dx = max(1e-12, max(box[0], box[1])/coordMax);

  // Define input points, including our false external generators.
  in.numberofpoints = numGenerators + 4;
  in.pointlist = new RealType[2*in.numberofpoints];
  copy(points.begin(), points.end(), in.pointlist);
  const double xmin = 0.5*(low[0] + high[0]) - boxsize;
  const double ymin = 0.5*(low[1] + high[1]) - boxsize;
  const double xmax = 0.5*(low[0] + high[0]) + boxsize;
  const double ymax = 0.5*(low[1] + high[1]) + boxsize;
  in.pointlist[2*numGenerators  ] = xmin; in.pointlist[2*numGenerators+1] = ymin;
  in.pointlist[2*numGenerators+2] = xmax; in.pointlist[2*numGenerators+3] = ymin;
  in.pointlist[2*numGenerators+4] = xmax; in.pointlist[2*numGenerators+5] = ymax;
  in.pointlist[2*numGenerators+6] = xmin; in.pointlist[2*numGenerators+7] = ymax;
  // cerr << "Chose bounding range : (" << xmin << " " << ymin << ") (" << xmax << " " << ymax << ")" << endl;

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // Fill in Triangle's boundary info.  We use our imposed outer box of fake
  // generators.
  in.numberofsegments = 4;
  in.segmentlist = new int[2*in.numberofsegments];
  j = 0;
  for (i = 0; i != 4; ++i) {
    in.segmentlist[j++] = numGenerators + i;
    in.segmentlist[j++] = numGenerators + ((i + 1) % 4);
  }
  in.segmentmarkerlist = 0;
  in.numberofholes = 0;
  in.holelist = 0;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;

  // Set up the structure for the triangulation.
  delaunay.pointlist = 0;
  delaunay.pointattributelist = 0;
  delaunay.pointmarkerlist = 0;
  delaunay.trianglelist = 0;
  delaunay.triangleattributelist = 0;
  delaunay.neighborlist = 0;
  delaunay.segmentlist = 0;
  delaunay.segmentmarkerlist = 0;
  delaunay.edgelist = 0;
  delaunay.edgemarkerlist = 0;
  delaunay.holelist = 0;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -c : Generates convex hull and places it in out.segmentlist.
  // -p : Uses the given PLC information.
  // if (geometry.empty())
  //   triangulate((char*)"Qzec", &in, &delaunay, 0);
  // else
  //   triangulate((char*)"Qzep", &in, &delaunay, 0);
  triangulate((char*)"Qzep", &in, &delaunay, 0);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    error("TriangleTessellator: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != numGenerators + 4) {
    char err[1024];
    snprintf(err, 1024, "TriangleTessellator: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, (int)numGenerators);
    error(err);
  }

  //--------------------------------------------------------
  // Create the Voronoi tessellation from the triangulation.
  //--------------------------------------------------------

  // Create the Voronoi nodes from the list of triangles. Each triangle 
  // has 3 nodes p, q, r, and corresponds to a Voronoi node at (X,Y), say.
  // The Voronoi node is located at the center of the triangle, though things
  // get a little squirrely at boundaries.  On boundary edges we create
  // a vertex at the edge center.

  // Find the circumcenters of each triangle, and build the set of triangles
  // associated with each generator.
  RealType  clow[2] = { numeric_limits<RealType>::max(),  numeric_limits<RealType>::max()};
  RealType chigh[2] = {-numeric_limits<RealType>::max(), -numeric_limits<RealType>::max()};
  CounterMap<EdgeHash> edgeCounter;
  vector<RealPoint> circumcenters(delaunay.numberoftriangles);
  map<int, set<int> > gen2tri;
  map<int, set<int> > neighbors;
  int k, pindex, qindex, rindex, iedge;
  for (i = 0; i != delaunay.numberoftriangles; ++i) {
    pindex = delaunay.trianglelist[3*i];
    qindex = delaunay.trianglelist[3*i + 1];
    rindex = delaunay.trianglelist[3*i + 2];
    if (pindex < numGenerators and qindex < numGenerators and rindex < numGenerators) {
      ++edgeCounter[hashEdge(pindex, qindex)];
      ++edgeCounter[hashEdge(qindex, rindex)];
      ++edgeCounter[hashEdge(rindex, pindex)];
    }
    computeCircumcenter(&delaunay.pointlist[2*pindex],
                        &delaunay.pointlist[2*qindex],
                        &delaunay.pointlist[2*rindex],
                        &circumcenters[i].x);
    gen2tri[pindex].insert(i);
    gen2tri[qindex].insert(i);
    gen2tri[rindex].insert(i);
    neighbors[pindex].insert(qindex);
    neighbors[pindex].insert(rindex);
    neighbors[qindex].insert(pindex);
    neighbors[qindex].insert(rindex);
    neighbors[rindex].insert(pindex);
    neighbors[rindex].insert(qindex);
    // cerr << "circumcenter : " << circumcenters[i] << endl;
    clow[0] = min(clow[0], circumcenters[i].x);
    clow[1] = min(clow[1], circumcenters[i].y);
    chigh[0] = max(chigh[0], circumcenters[i].x);
    chigh[1] = max(chigh[1], circumcenters[i].y);
  }
  ASSERT(circumcenters.size() == delaunay.numberoftriangles);
  ASSERT(clow[0] < chigh[0] and clow[1] < chigh[1]);
  RealType cbox[2] = {chigh[0] - clow[0], 
                      chigh[1] - clow[1]};
  const double cboxsize = 2.0*max(cbox[0], cbox[1]);
  const double cdx = max(1e-10, max(cbox[0], cbox[1])/coordMax);

  // Flag any generators on the edge of the tessellation.  Here we mean the actual
  // generators, not our added boundary ones.
  vector<bool> exteriorGenerators(numGenerators, false);
  BoolMap<EdgeHash> exteriorEdgeTest;
  list<EdgeHash> exteriorEdges;
  vector<vector<EdgeHash> > exteriorEdgesOfGen(numGenerators);
  for (typename CounterMap<EdgeHash>::const_iterator itr = edgeCounter.begin();
       itr != edgeCounter.end();
       ++itr) {
    ASSERT(itr->second == 1 or itr->second == 2);
    if (itr->second == 1) {
      i = itr->first.first;
      j = itr->first.second;
      exteriorGenerators[i] = true;
      exteriorGenerators[j] = true;
      exteriorEdges.push_back(itr->first);
      exteriorEdgeTest.insert(make_pair(itr->first, true));
      exteriorEdgesOfGen[i].push_back(itr->first);
      exteriorEdgesOfGen[j].push_back(itr->first);
    }
  }

  // Build the polygon representing our boundaries.
  BGpolygon boundary;
  if (geometry.empty()) {
    // The user did not provide a boundary, so use our local edge use count to 
    // find the bounding edges.  Note in this case there will be no holes.
    const unsigned numBoundaryPoints = exteriorEdges.size();
    vector<RealPoint> boundaryPoints;
    vector<EdgeHash> boundaryEdgeOrder;
    boundaryEdgeOrder.push_back(exteriorEdges.front());
    exteriorEdges.pop_front();
    i = boundaryEdgeOrder.front().first;
    boundaryPoints.reserve(numBoundaryPoints + 1);
    boundaryPoints.push_back(RealPoint(points[2*i], points[2*i+1]));
    i = boundaryEdgeOrder.front().second;
    boundaryPoints.push_back(RealPoint(points[2*i], points[2*i+1]));
    while (exteriorEdges.size() > 0) {
      list<EdgeHash>::iterator itr = find_if(exteriorEdges.begin(), exteriorEdges.end(),
                                             MatchEitherPairValue<int>(i));
      ASSERT(itr != exteriorEdges.end());
      i = (itr->first == i ? itr->second : itr->first);
      boundaryEdgeOrder.push_back(*itr);
      exteriorEdges.erase(itr);
      boundaryPoints.push_back(RealPoint(points[2*i], points[2*i+1]));
    }
    ASSERT(boundaryEdgeOrder.size() == numBoundaryPoints);
    ASSERT(MatchEitherPairValue<int>(boundaryEdgeOrder.front().first)(boundaryEdgeOrder.back()));
    ASSERT(boundaryPoints.size() == numBoundaryPoints + 1);
    ASSERT(boundaryPoints.front() == boundaryPoints.back());
    boost::geometry::assign(boundary, BGring(boundaryPoints.begin(), boundaryPoints.end()));

  } else {
    // Copy the PLC provided boundary information into a Boost.Geometry polygon.
    vector<RealPoint> boundaryPoints;
    boundaryPoints.reserve(geometry.facets.size() + 1);
    i = geometry.facets[0][0];
    boundaryPoints.push_back(RealPoint(PLCpoints[2*i], PLCpoints[2*i+1]));
    for (j = 0; j != geometry.facets.size(); ++j) {
      ASSERT(geometry.facets[j].size() == 2);
      i =  geometry.facets[j][1];
      boundaryPoints.push_back(RealPoint(PLCpoints[2*i], PLCpoints[2*i+1]));
    }
    ASSERT(boundaryPoints.size() == geometry.facets.size() + 1);
    ASSERT(boundaryPoints.front() == boundaryPoints.back());
    boost::geometry::assign(boundary, BGring(boundaryPoints.begin(), boundaryPoints.end()));

    // Add any interior holes.
    const unsigned numHoles = geometry.holes.size();
    if (numHoles > 0) {
      typename BGpolygon::inner_container_type& holes = boundary.inners();
      holes.resize(numHoles);
      for (k = 0; k != numHoles; ++k) {
        boundaryPoints = vector<RealPoint>();
        boundaryPoints.reserve(geometry.holes[k].size() + 1);
        i = geometry.holes[k][0][0];
        boundaryPoints.push_back(RealPoint(PLCpoints[2*i], PLCpoints[2*i+1]));
        for (j = 0; j != geometry.holes[k].size(); ++j) {
          ASSERT(geometry.holes[k][j].size() == 2);
          i =  geometry.holes[k][j][1];
          boundaryPoints.push_back(RealPoint(PLCpoints[2*i], PLCpoints[2*i+1]));
        }
        ASSERT(boundaryPoints.size() == geometry.holes[k].size() + 1);
        ASSERT(boundaryPoints.front() == boundaryPoints.back());
        boost::geometry::assign(holes[k], BGring(boundaryPoints.begin(), boundaryPoints.end()));
      }
    }
  }

  // Walk each generator and build up it's unique nodes and faces.
  mesh.cells.resize(numGenerators);
  bool inside;
  double minR, thpt;
  RealPoint X;
  IntPoint pX1, pX2;
  map<IntPoint, int> point2node;
  map<EdgeHash, int> edgeHash2id;
  map<int, vector<int> > edgeCells;
  map<int, BGring> cellRings;
  for (i = 0; i != numGenerators; ++i) {

    // Add the circumcenters as points for the cell.
    set<IntPoint> cellPointSet;
    for (set<int>::const_iterator triItr = gen2tri[i].begin();
         triItr != gen2tri[i].end();
         ++triItr) {
      cellPointSet.insert(IntPoint(circumcenters[*triItr].x - clow[0],
                                   circumcenters[*triItr].y - clow[1],
                                   cdx));
      // cellPoints.push_back(circumcenters[*triItr]);
      // cerr << "Cell " << i << " adding circumcenter " << cellPoints.back() << endl;
    }
    ASSERT(cellPointSet.size() >= 3);

    // // Build the convex hull of the cell points.
    // vector<double> cellPointCoords;
    // for (j = 0; j != cellPoints.size(); ++j) {
    //   cellPointCoords.push_back(cellPoints[j].x);
    //   cellPointCoords.push_back(cellPoints[j].y);
    // }
    // ASSERT(cellPointCoords.size() == 2*cellPoints.size());
    // PLC<2, double> hull = convexHull_2d<double>(cellPointCoords, low, dx);
    // ASSERT(hull.facets.size() >= 3);
    // ASSERT(hull.facets[0][0] < cellPoints.size());
    // vector<RealPoint> ringPoints;
    // ringPoints.push_back(cellPoints[hull.facets[0][0]]);
    // for (j = 0; j != hull.facets.size(); ++j) {
    //   ASSERT(hull.facets[j].size() == 2);
    //   ASSERT(hull.facets[j][1] < cellPoints.size());
    //   ringPoints.push_back(cellPoints[hull.facets[j][1]]);
    // }
    // cellRings[i] = BGring(ringPoints.begin(), ringPoints.end());

    // Build the convex hull of the cell points.
    BGmulti_point mpoints;
    for (typename set<IntPoint>::const_iterator itr = cellPointSet.begin();
         itr != cellPointSet.end();
         ++itr) mpoints.push_back(RealPoint(itr->realx(clow[0], cdx),
                                            itr->realy(clow[1], cdx)));
    boost::geometry::convex_hull(mpoints, cellRings[i]);

    // Intersect with the boundary to get the bounded cell.
    // Since for complex boundaries this may return more than one polygon, we find
    // the one that contains the generator.
    vector<BGring> cellIntersections;
    boost::geometry::intersection(boundary, cellRings[i], cellIntersections);
    if (cellIntersections.size() == 0) {
      cerr << points[2*i] << " " << points[2*i+1] << endl 
           << boost::geometry::dsv(cellRings[i]) << endl
           << boost::geometry::dsv(mpoints) << endl
           << boost::geometry::dsv(boundary) << endl;
    }
    ASSERT(cellIntersections.size() > 0);    
    if (cellIntersections.size() == 1) {
      cellRings[i] = cellIntersections[0];
    } else {
      X = RealPoint(points[2*i], points[2*i+1]);
      inside = boost::geometry::within(X, cellIntersections[j]);
      minR = numeric_limits<double>::max();
      j = 0;
      for (j = 0; j != cellIntersections.size(); ++j) {
        thpt = boost::geometry::distance(X, cellIntersections[j]);
        if (thpt < minR) {
          k = j;
          minR = thpt;
        }
      }
      ASSERT(k < cellIntersections.size());
      cellRings[i] = cellIntersections[k];
    }

    // Now build the unique mesh nodes and cell info.
    for (typename BGring::const_iterator itr = cellRings[i].begin();
         itr != cellRings[i].end() - 1;
         ++itr) {
      pX1 = IntPoint(itr->x - clow[0],
                     itr->y - clow[1], cdx);
      pX2 = IntPoint((itr+1)->x - clow[0],
                     (itr+1)->y - clow[1], cdx);
      if (pX1 != pX2) {
        j = addKeyToMap(pX1, point2node);
        k = addKeyToMap(pX2, point2node);
        ASSERT(j != k);
        iedge = addKeyToMap(hashEdge(j, k), edgeHash2id);
        edgeCells[iedge].push_back(j < k ? i : ~i);
        mesh.cells[i].push_back(j < k ? iedge : ~iedge);
        // cerr << "Cell " << i << " adding edge " << iedge << " : " << *itr << " " << *(itr + 1) << endl;
      }
    }
    ASSERT(mesh.cells[i].size() >= 3);
  }
  ASSERT(edgeCells.size() == edgeHash2id.size());

  // Fill in the mesh nodes.
  mesh.nodes = vector<RealType>(2*point2node.size());
  for (typename map<IntPoint, int>::const_iterator itr = point2node.begin();
       itr != point2node.end();
       ++itr) {
    const IntPoint& p = itr->first;
    i = itr->second;
    ASSERT(i < mesh.nodes.size()/2);
    mesh.nodes[2*i]   = p.realx(clow[0], cdx);
    mesh.nodes[2*i+1] = p.realy(clow[1], cdx);
  }

  // Fill in the mesh faces.
  mesh.faces = vector<vector<unsigned> >(edgeHash2id.size());
  for (typename map<EdgeHash, int>::const_iterator itr = edgeHash2id.begin();
       itr != edgeHash2id.end();
       ++itr) {
    const EdgeHash& ehash = itr->first;
    i = itr->second;
    ASSERT(i < mesh.faces.size());
    ASSERT(mesh.faces[i].size() == 0);
    mesh.faces[i].push_back(ehash.first);
    mesh.faces[i].push_back(ehash.second);
  }

  // Fill in the mesh faceCells.
  mesh.faceCells = vector<vector<int> >(mesh.faces.size());
  for (i = 0; i != mesh.faces.size(); ++i) {
    // if (not(edgeCells[i].size() == 1 or edgeCells[i].size() == 2)) {
    //   cerr << "Blago! " << i << " " << edgeCells[i].size() << endl;
    //   for (j = 0; j != edgeCells[i].size(); ++j) cerr << " --> " << edgeCells[i][j] << " " << points[2*edgeCells[i][j]] << " " << points[2*edgeCells[i][j]+1] << endl;
    // }
    ASSERT(edgeCells[i].size() == 1 or edgeCells[i].size() == 2);
    mesh.faceCells[i] = edgeCells[i];
  }

  // Clean up.
  delete [] in.pointlist;
  trifree((VOID*)delaunay.pointlist);
  trifree((VOID*)delaunay.pointmarkerlist);
  trifree((VOID*)delaunay.trianglelist);
  trifree((VOID*)delaunay.edgelist);
  trifree((VOID*)delaunay.edgemarkerlist);
  if (geometry.empty())
  {
    trifree((VOID*)delaunay.segmentlist);
    trifree((VOID*)delaunay.segmentmarkerlist);
  }
  else
  {
    delete [] in.segmentlist;
    delete [] in.holelist;
    trifree((VOID*) delaunay.segmentlist);
    trifree((VOID*) delaunay.segmentmarkerlist);
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class TriangleTessellator<double>;

}
