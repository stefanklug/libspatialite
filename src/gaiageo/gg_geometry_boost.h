#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/recursive_wrapper.hpp>
#include <boost/variant/variant.hpp>
#include <vector>


//::boost::geometry::model::d2::point_xy

class point {
public:
	point(double _x=0.0, double _y=0.0, double _z=0.0, double _m=0.0):x(_x), y(_y), z(_z) {}

	double x,y,z,m;
	int DimensionModel;
};

BOOST_GEOMETRY_REGISTER_POINT_2D(point, double, ::boost::geometry::cs::cartesian, x, y)

typedef ::boost::geometry::model::box<point> box;

class linearring: public ::boost::geometry::model::ring<point> {
public:
	int DimensionModel;
};

BOOST_GEOMETRY_REGISTER_RING(::linearring)

class linestring: public ::boost::geometry::model::linestring<point> {
public:
	int DimensionModel;
};

BOOST_GEOMETRY_REGISTER_LINESTRING(linestring)


class polygon {
public:
	// Member types
	typedef ::linearring ring_type;
	typedef std::vector< ::linearring> inner_container_type;

	inline ring_type const& outer() const { return m_outer; }
	inline inner_container_type const& inners() const { return m_inners; }

	inline ring_type& outer() { return m_outer; }
	inline inner_container_type & inners() { return m_inners; }

	/// Utility method, clears outer and inner rings
	inline void clear()
	{
		m_outer.clear();
		m_inners.clear();
	}

	int DimensionModel;

private:
	ring_type m_outer;
	inner_container_type m_inners;

};

namespace boost { namespace geometry { namespace traits
{
	template<> struct tag< ::polygon > { typedef polygon_tag type; };
	template<> struct ring_const_type< ::polygon > { typedef typename ::polygon::ring_type const& type; };
	template<> struct ring_mutable_type< ::polygon > { typedef typename ::polygon::ring_type & type; };
	template<> struct interior_const_type< ::polygon > { typedef typename ::polygon::inner_container_type const& type; };
	template<> struct interior_mutable_type< ::polygon > { typedef typename ::polygon::inner_container_type & type; };

	template<> struct exterior_ring< ::polygon> {
		static inline typename ::polygon::ring_type& get(::polygon& p) { return p.outer(); }
		static inline typename ::polygon::ring_type const& get(::polygon const& p) { return p.outer(); }
	};

	template<> struct interior_rings< ::polygon> {
		static inline typename ::polygon::inner_container_type& get(::polygon& p) { return p.inners(); }
		static inline typename ::polygon::inner_container_type const& get(::polygon const& p) { return p.inners(); }
	};

}}} // namespace traits

typedef ::boost::geometry::model::multi_point< ::point> multi_point;
typedef ::boost::geometry::model::multi_linestring< ::linestring> multi_linestring;
typedef ::boost::geometry::model::multi_polygon< ::polygon> multi_polygon;

struct geometry_collection {

	geometry_collection() {
		Srid = 0;
		DimensionModel = 0;
		DeclaredType = 0;
	}

	multi_point points;
	multi_linestring linestrings;
	multi_polygon polygons;

	int Srid;		/* the SRID value for this GEOMETRY */
	int DimensionModel;	/* (x,y), (x,y,z), (x,y,m) or (x,y,z,m) */
	int DeclaredType; /** any valid Geometry Class type */
};

/*typedef ::boost::variant<
  point,
  linestring,
  polygon,
  multi_point,
  multi_linestring,
  multi_polygon,
  geometry_collection
> geometry;
*/

