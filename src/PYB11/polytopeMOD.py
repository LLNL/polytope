"Master module for polytope"

from PYB11Generator import *

PYB11includes = ['"polytope.hh"',
                 '"polytope_write_OOGL.hh"',
                 '"polytope_pybind11_helpers.hh"']

PYB11namespaces = ["polytope"]

#-------------------------------------------------------------------------------
# Add the polytope classes
#-------------------------------------------------------------------------------
from PLC import *
from ReducedPLC import *
from Tessellation import *
from Tessellator import *
from Tessellators import *
from MeshEditor import *

#-------------------------------------------------------------------------------
# Polytope functions
#-------------------------------------------------------------------------------
@PYB11template("RealType")
def writePLCtoOFF(plc = "const PLC<%(Dimension)s, %(RealType)s>&",
                  coords = "const std::vector<%(RealType)s>&",
                  filename = "const std::string"):
    "Write a %(Dimension)sD OFF polylist file."
    return "void"

writePLCtoOFF2d = PYB11TemplateFunction(writePLCtoOFF, pyname="writePLCtoOFF", template_parameters={"RealType": "double",
                                                                                                    "Dimension": "2"})
writePLCtoOFF3d = PYB11TemplateFunction(writePLCtoOFF, pyname="writePLCtoOFF", template_parameters={"RealType": "double",
                                                                                                    "Dimension": "3"})

#...............................................................................
@PYB11template("int Dimension", "RealType")
@PYB11implementation("""[](const Tessellation<%(Dimension)s, %(RealType)s>& mesh,
                           std::string filePrefix,
                           py::dict nodeFieldsDict,
                           py::dict edgeFieldsDict,
                           py::dict faceFieldsDict,
                           py::dict cellFieldsDict,
                           int cycle,
                           %(RealType)s time) {
#ifdef HAVE_SILO
                               std::map<std::string, %(RealType)s*> nodeFields, edgeFields, faceFields, cellFields;
                               auto nodeVals = pybind11_helpers::copyDictToMap<std::string, double>(nodeFieldsDict);
                               auto edgeVals = pybind11_helpers::copyDictToMap<std::string, double>(edgeFieldsDict);
                               auto faceVals = pybind11_helpers::copyDictToMap<std::string, double>(faceFieldsDict);
                               auto cellVals = pybind11_helpers::copyDictToMap<std::string, double>(cellFieldsDict);
                               for (auto& kv: nodeVals) nodeFields[kv.first] = &(kv.second.front());
                               for (auto& kv: edgeVals) edgeFields[kv.first] = &(kv.second.front());
                               for (auto& kv: faceVals) faceFields[kv.first] = &(kv.second.front());
                               for (auto& kv: cellVals) cellFields[kv.first] = &(kv.second.front());
                               SiloWriter<%(Dimension)s, %(RealType)s>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, cycle, time);
#else
                               throw std::runtime_error("Polytope built without SILO support");
#endif
                           }""")
def writeTessellation(mesh = "const Tessellation<%(Dimension)s, %(RealType)s>&",
                      filePrefix = "std::string",
                      nodeFieldsDict = "py::dict",
                      edgeFieldsDict = "py::dict",
                      faceFieldsDict = "py::dict",
                      cellFieldsDict = "py::dict",
                      cycle = "int",
                      time = "%(RealType)s"):
    "Write a tessellation (and arbitrary fields) to a silo file."
    return "void"

writeTessellation2d = PYB11TemplateFunction(writeTessellation, template_parameters=("2", "double"))
writeTessellation3d = PYB11TemplateFunction(writeTessellation, template_parameters=("3", "double"))

#...............................................................................
@PYB11template("int Dimension", "RealType")
@PYB11implementation("""[](const Tessellation<%(Dimension)s, %(RealType)s>& mesh,
                           std::string filePrefix,
                           std::string directory,
                           py::dict nodeFieldsDict,
                           py::dict edgeFieldsDict,
                           py::dict faceFieldsDict,
                           py::dict cellFieldsDict,
                           int cycle,
                           %(RealType)s time) {
#ifdef HAVE_SILO
                               std::map<std::string, %(RealType)s*> nodeFields, edgeFields, faceFields, cellFields;
                               auto nodeVals = pybind11_helpers::copyDictToMap<std::string, double>(nodeFieldsDict);
                               auto edgeVals = pybind11_helpers::copyDictToMap<std::string, double>(edgeFieldsDict);
                               auto faceVals = pybind11_helpers::copyDictToMap<std::string, double>(faceFieldsDict);
                               auto cellVals = pybind11_helpers::copyDictToMap<std::string, double>(cellFieldsDict);
                               for (auto& kv: nodeVals) nodeFields[kv.first] = &(kv.second.front());
                               for (auto& kv: edgeVals) edgeFields[kv.first] = &(kv.second.front());
                               for (auto& kv: faceVals) faceFields[kv.first] = &(kv.second.front());
                               for (auto& kv: cellVals) cellFields[kv.first] = &(kv.second.front());
                               SiloWriter<%(Dimension)s, %(RealType)s>::write(mesh, nodeFields, edgeFields, faceFields, cellFields, filePrefix, directory, cycle, time);
#else
                               throw std::runtime_error("Polytope built without SILO support");
#endif
                           }""")
def writeTessellationDir(mesh = "const Tessellation<%(Dimension)s, %(RealType)s>&",
                         filePrefix = "std::string",
                         directory = "std::string",
                         nodeFieldsDict = "py::dict",
                         edgeFieldsDict = "py::dict",
                         faceFieldsDict = "py::dict",
                         cellFieldsDict = "py::dict",
                         cycle = "int",
                         time = "%(RealType)s"):
    "Write a tessellation (and arbitrary fields) to a silo file.  This version supports specifying the output directory."
    return "void"

writeTessellationDir2d = PYB11TemplateFunction(writeTessellationDir, template_parameters=("2", "double"), pyname="writeTessellation2d")
writeTessellationDir3d = PYB11TemplateFunction(writeTessellationDir, template_parameters=("3", "double"), pyname="writeTessellation3d")

#...............................................................................
@PYB11template("RealType")
@PYB11implementation("""[](py::list& points,
                           py::tuple& low,
                           %(RealType)s& dx) {
                               const auto coords = pybind11_helpers::copyCoords<2, %(RealType)s>(points);
                               std::vector<%(RealType)s> vlow = {low[0].cast<%(RealType)s>(), low[1].cast<%(RealType)s>()};
                               auto result = convexHull_2d(coords, &vlow.front(), dx);
                               return result;
                           }""")
def convexHull2d(points = "py::list",
                 low = "py::tuple",
                 dx = "%(RealType)s&"):
    """Return the convex hull of the points in a list of tuple positions of the form
  points = [(x,y), (x,y), ...]
  low: a tuple of coordinates (x,y) which should be less than any of the coordinates in "points"
  dx:  the precision for distinguishing between coordinate values
"""
    return "PLC<2, %(RealType)s>"

constructConvexHull2d = PYB11TemplateFunction(convexHull2d, template_parameters="double")

#...............................................................................
@PYB11template("RealType")
@PYB11implementation("""[](py::list& points,
                           py::tuple& low,
                           %(RealType)s& dx) {
                               const auto coords = pybind11_helpers::copyCoords<3, %(RealType)s>(points);
                               std::vector<%(RealType)s> vlow = {low[0].cast<%(RealType)s>(), 
                                                                 low[1].cast<%(RealType)s>(),
                                                                 low[2].cast<%(RealType)s>()};
                               auto result = convexHull_3d(coords, &vlow.front(), dx);
                               return result;
                           }""")
def convexHull3d(points = "py::list",
                 low = "py::tuple",
                 dx = "%(RealType)s&"):
    """Return the convex hull of the points in a list of tuple positions of the form
  points = [(x,y,z), (x,y,z), ...]
  low: a tuple of coordinates (x,y,z) which should be less than any of the coordinates in "points"
  dx:  the precision for distinguishing between coordinate values
"""
    return "PLC<3, %(RealType)s>"

constructConvexHull3d = PYB11TemplateFunction(convexHull3d, template_parameters="double")
