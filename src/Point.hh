//----------------------------------------------------------------------------//
// 2D and 3D integral Point types used internally in polytope.  Not really
// for external consumption!.
//----------------------------------------------------------------------------//
#ifndef __polytope_Point__
#define __polytope_Point__

#include <iostream>
#include <iterator>

#include "polytope_serialize.hh"
#include "polytope.hh"

namespace polytope {

//------------------------------------------------------------------------------
// A integer version of the simple 2D point.
//------------------------------------------------------------------------------
template<typename UintType>
struct Point2 {
  UintType x, y;
  unsigned index;
  Point2(): x(0), y(0), index(0) {}
  Point2(const UintType& xi, const UintType& yi, const unsigned i = 0): x(xi), y(yi), index(i) {}
  Point2& operator=(const Point2& rhs) { x = rhs.x; y = rhs.y; index = rhs.index; return *this; }
  bool operator==(const Point2& rhs) const { return (x == rhs.x and y == rhs.y); }
  bool operator<(const Point2& rhs) const {
    return (x < rhs.x                ? true :
            x == rhs.x and y < rhs.y ? true :
            false);
  }
  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, const RealType& dx, const unsigned i = 0): 
    x(static_cast<UintType>(xi/dx + 0.5)),
    y(static_cast<UintType>(yi/dx + 0.5)),
    index(i) {}
  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) { return static_cast<RealType>(y*dy) + ymin; }
  Point2& operator+=(const Point2& rhs) { x += rhs.x; y += rhs.y; return *this; }
  Point2& operator/=(const UintType& rhs) { x /= rhs; y /= rhs; return *this; }
};

// It's nice being able to print these things.
template<typename UintType>
std::ostream&
operator<<(std::ostream& os, const Point2<UintType>& p) {
  os << "(" << p.x << " " << p.y << ")(" << p.index << ")";
  return os;
}

// Serialization.
template<typename UintType>
struct Serializer<Point2<UintType> > {

  static void serializeImpl(const Point2<UintType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point2<UintType>& value,
                              std::vector<char>::const_iterator bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.x, bufItr, endItr);
    deserialize(value.y, bufItr, endItr);
    deserialize(value.index, bufItr, endItr);
  }
};

//------------------------------------------------------------------------------
// A integer version of the simple 3D point.
//------------------------------------------------------------------------------
template<typename UintType>
struct Point3 {
  UintType x, y, z;
  unsigned index;
  Point3(): x(0), y(0), z(0), index(0) {}
  Point3(const UintType& xi, const UintType& yi, const UintType& zi, const unsigned i = 0): x(xi), y(yi), z(zi), index(i) {}
  Point3& operator=(const Point3& rhs) { x = rhs.x; y = rhs.y; z = rhs.z; index = rhs.index; return *this; }
  bool operator==(const Point3& rhs) const { return (x == rhs.x and y == rhs.y and z == rhs.z); }
  bool operator<(const Point3& rhs) const {
    return (x < rhs.x                               ? true :
            x == rhs.x and y < rhs.y                ? true :
            x == rhs.x and y == rhs.y and z < rhs.z ? true :
            false);
  }
  template<typename RealType>
  Point3(const RealType& xi, const RealType& yi, const RealType& zi, const RealType& dx, const unsigned i = 0): 
    x(static_cast<UintType>(xi/dx + 0.5)),
    y(static_cast<UintType>(yi/dx + 0.5)),
    z(static_cast<UintType>(zi/dx + 0.5)),
    index(i) {}
  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) { return static_cast<RealType>(y*dy) + ymin; }
  template<typename RealType> RealType realz(const RealType& zmin, const RealType& dz) { return static_cast<RealType>(z*dz) + zmin; }
  Point3& operator+=(const Point3& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  Point3& operator/=(const UintType& rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }
};

// It's nice being able to print these things.
template<typename UintType>
std::ostream&
operator<<(std::ostream& os, const Point3<UintType>& p) {
  os << "(" << p.x << " " << p.y << " " << p.z <<  ")(" << p.index << ")";
  return os;
}

// Serialization.
template<typename UintType>
struct Serializer<Point3<UintType> > {

  static void serializeImpl(const Point3<UintType>& value,
                            std::vector<char>& buffer) {
    serialize(value.x, buffer);
    serialize(value.y, buffer);
    serialize(value.z, buffer);
    serialize(value.index, buffer);
  }

  static void deserializeImpl(Point3<UintType>& value,
                              std::vector<char>::const_iterator bufItr,
                              const std::vector<char>::const_iterator endItr) {
    deserialize(value.x, bufItr, endItr);
    deserialize(value.y, bufItr, endItr);
    deserialize(value.z, bufItr, endItr);
    deserialize(value.index, bufItr, endItr);
  }
};

}

#endif