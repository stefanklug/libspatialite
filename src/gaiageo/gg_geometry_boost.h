#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/recursive_wrapper.hpp>
#include <boost/variant/variant.hpp>
#include <vector>



typedef ::boost::geometry::model::point<double, 3, ::boost::geometry::cs::cartesian> point;
typedef ::boost::geometry::model::box<point> box;
typedef ::boost::geometry::model::ring<point> linearring;
typedef ::boost::geometry::model::linestring<point> linestring;
typedef ::boost::geometry::model::polygon<point> polygon;
typedef ::boost::geometry::model::multi_point<point> multi_point;
typedef ::boost::geometry::model::multi_linestring<linestring> multi_linestring;
typedef ::boost::geometry::model::multi_polygon<polygon> multi_polygon;

struct geometry_collection {
	multi_point points;
	multi_linestring linestrings;
	multi_polygon polygons;
};

typedef ::boost::variant<
  point,
  linestring,
  polygon,
  multi_point,
  multi_linestring,
  multi_polygon,
  geometry_collection
> geometry;


