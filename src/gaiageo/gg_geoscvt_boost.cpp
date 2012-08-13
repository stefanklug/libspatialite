/*

 gg_geoscvt.c -- Gaia / GEOS conversion [Geometry]
    
 version 3.0, 2011 July 20

 Author: Sandro Furieri a.furieri@lqt.it

 ------------------------------------------------------------------------------
 
 Version: MPL 1.1/GPL 2.0/LGPL 2.1
 
 The contents of this file are subject to the Mozilla Public License Version
 1.1 (the "License"); you may not use this file except in compliance with
 the License. You may obtain a copy of the License at
 http://www.mozilla.org/MPL/
 
Software distributed under the License is distributed on an "AS IS" basis,
WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
for the specific language governing rights and limitations under the
License.

The Original Code is the SpatiaLite library

The Initial Developer of the Original Code is Alessandro Furieri
 
Portions created by the Initial Developer are Copyright (C) 2008
the Initial Developer. All Rights Reserved.

Contributor(s):

Alternatively, the contents of this file may be used under the terms of
either the GNU General Public License Version 2 or later (the "GPL"), or
the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
in which case the provisions of the GPL or the LGPL are applicable instead
of those above. If you wish to allow use of your version of this file only
under the terms of either the GPL or the LGPL, and not to allow others to
use your version of this file under the terms of the MPL, indicate your
decision by deleting the provisions above and replace them with the notice
and other provisions required by the GPL or the LGPL. If you do not delete
the provisions above, a recipient may use your version of this file under
the terms of any one of the MPL, the GPL or the LGPL.
 
*/

#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef OMIT_BOOSTGEOMETRY		/* including Boost */
#include "gg_geometry_boost.h"
#endif

#ifdef SPL_AMALGAMATION		/* spatialite-amalgamation */
#include <spatialite/sqlite3ext.h>
#else
#include <sqlite3ext.h>
#endif

#include <spatialite/gaiageo.h>

#ifndef OMIT_BOOSTGEOMETRY		/* including Boost */

namespace bg = boost::geometry;

static inline point toBoostPoint(const gaiaPointPtr pt) {
	switch(pt->DimensionModel) {
	case GAIA_XY_Z:
    case GAIA_XY_Z_M:
    	return point(pt->X, pt->Y, pt->Z);
    default:
    	return point(pt->X, pt->Y, 0.0);
	}
}

static linestring toBoostLinestring(const gaiaLinestringPtr ln) {
	double x;
	double y;
	double z;
	double m;
	double* v;
	int iv;

	linestring res;

	v = ln->Coords;
	for (iv = 0; iv < ln->Points; iv++)
	{
		switch(ln->DimensionModel) {
		{
		  case GAIA_XY_Z:
			  res.push_back(point(v[0], v[1], v[2]));
			  v += 3;
			  break;
		  case GAIA_XY_M:
			  res.push_back(point(v[0], v[1], 0.0));
			  v += 3;
			  break;
		  case GAIA_XY_Z_M:
			  res.push_back(point(v[0], v[1], v[2]));
			  v += 4;
			  break;
		  default:
			  res.push_back(point(v[0], v[1], 0.0));
			  v += 2;
			  break;
		  };
		}
	}

	return res;
}

static linearring toBoostRing(const gaiaRingPtr ln) {
	double x;
	double y;
	double z;
	double m;
	double* v;
	int iv;

	linearring res;

	v = ln->Coords;
	for (iv = 0; iv < ln->Points; iv++)
	{
		switch(ln->DimensionModel) {
		{
		  case GAIA_XY_Z:
			  res.push_back(point(v[0], v[1], v[2]));
			  v += 3;
			  break;
		  case GAIA_XY_M:
			  res.push_back(point(v[0], v[1], 0.0));
			  v += 3;
			  break;
		  case GAIA_XY_Z_M:
			  res.push_back(point(v[0], v[1], v[2]));
			  v += 4;
			  break;
		  default:
			  res.push_back(point(v[0], v[1], 0.0));
			  v += 2;
			  break;
		  };
		}
	}

	return res;
}

static polygon toBoostPolygon(const gaiaPolygonPtr pg) {
	polygon res;
	linearring ring;

	ring = toBoostRing(pg->Exterior);
	res.outer().swap(ring);

	for(i=0; i<pg->NumInteriors; i++) {
		res.inners().push_back(toBoostRing(pg->Interiors[i]));
	}

	return res;
}


static multi_point toBoostMultiPoint(const gaiaPointPtr pt) {
	multi_point res;
	while(pt) {
		res.push_back(toBoostPoint(pt));
		pt = pt->Next;
	}
	return res;
}

static multi_linestring toBoostMultiLinestring(const gaiaLinestringPtr ln) {
	multi_linestring res;
	while(ln) {
		res.push_back(toBoostLinestring(ln));
		ln = ln->Next;
	}
	return res;
}

static multi_polygon toBoostMultiPolygon(const gaiaPolygonPtr pg) {
	multi_polygon res;
	while(pg) {
		res.push_back(toBoostPolygon(pg));
		pg = pg->Next;
	}
	return res;
}

static geometry
toBoostGeometry (const gaiaGeomCollPtr gaia)
{
/* converting a GAIA Geometry into a GEOS Geometry */
    int pts = 0;
    int lns = 0;
    int pgs = 0;
    int type;
    int geos_type;
    unsigned int dims;
    int iv;
    int ib;
    int nItem;
    double x;
    double y;
    double z;
    double m;
    gaiaPointPtr pt;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;
    gaiaRingPtr rng;

    geometry res;

    if (!gaia)
	return res;
    type = gaiaGeometryAliasType(gaia);

    switch (type) {
      case GAIA_POINT:
    	  res = toBoostPoint(gaia->FirstPoint);
	  break;


      case GAIA_LINESTRING:
    	  res = toBoostLinestring(gaia->FirstLinestring);
	  break;

      case GAIA_POLYGON:
    	  res = toBoostPolygon(gaia->FirstPolygon);
	  break;

      case GAIA_MULTIPOINT:
    	  res = toBoostMultiPoint(gaia->FirstPoint);
    	  break;
      case GAIA_MULTILINESTRING:
    	  res = toBoostMultiLinestring(gaia->FirstLinestring);
    	  break;
      case GAIA_MULTIPOLYGON:
    	  res = toBoostMultiPolygon(gaia->FirstPolygon);
    	  break;
      case GAIA_GEOMETRYCOLLECTION:
    	  geometry_collection c;
		  pt = gaia->FirstPoint;
		  while(pt) {
			  c.points.push_back(toBoostPoint(pt));
			  pt = pt->Next;
		  }

		  ln = gaia->FirstLinestring;
		  while(ln) {
			  c.linestrings.push_back(toBoostLinestring(ln));
			  ln = ln->Next;
		  }

		  ln = gaia->FirstPolygon;
		  while(pg) {
			  c.polygons.push_back(toBoostPolygon(pg));
			  pg = pg->Next;
		  }

		  res = c;
		  break;
      };
    //if (geos)
	//GEOSSetSRID (geos, gaia->Srid);
    return res;
}

struct BoostToGaiaCollectionVisitor: public boost::static_visitor<gaiaGeomCollPtr>
{

	BoostToGaiaCollectionVisitor()

    void operator()(geometry_collection& v) const
    {
        i *= 2;
    }

    void operator()(multi_polygon& v) const
    {
        str += str;
    }

    void operator()(polygon& v) const
	{
		str += str;
	}

    void operator()(polygon& v) const
	{
		str += str;
	}

private:


};

static gaiaGeomCollPtr
fromBoostGeometry (const geometry& geom, const int dimension_model)
{
/* converting a GEOS Geometry into a GAIA Geometry */
    int type;
    int itemType;
    unsigned int dims;
    int iv;
    int ib;
    int it;
    int sub_it;
    int nItems;
    int nSubItems;
    int holes;
    unsigned int points;
    double x;
    double y;
    double z;
    const GEOSCoordSequence *cs;
    const GEOSGeometry *geos_ring;
    const GEOSGeometry *geos_item;
    const GEOSGeometry *geos_sub_item;
    gaiaGeomCollPtr gaia = NULL;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;
    gaiaRingPtr rng;
    if (!geos)
	return NULL;
    type = GEOSGeomTypeId (geos);
    switch (type)
      {
      case GEOS_POINT:
	  if (dimension_model == GAIA_XY_Z)
	      gaia = gaiaAllocGeomCollXYZ ();
	  else if (dimension_model == GAIA_XY_M)
	      gaia = gaiaAllocGeomCollXYM ();
	  else if (dimension_model == GAIA_XY_Z_M)
	      gaia = gaiaAllocGeomCollXYZM ();
	  else
	      gaia = gaiaAllocGeomColl ();
	  gaia->DeclaredType = GAIA_POINT;
	  gaia->Srid = GEOSGetSRID (geos);
	  cs = GEOSGeom_getCoordSeq (geos);
	  GEOSCoordSeq_getDimensions (cs, &dims);
	  if (dims == 3)
	    {
		GEOSCoordSeq_getX (cs, 0, &x);
		GEOSCoordSeq_getY (cs, 0, &y);
		GEOSCoordSeq_getZ (cs, 0, &z);
	    }
	  else
	    {
		GEOSCoordSeq_getX (cs, 0, &x);
		GEOSCoordSeq_getY (cs, 0, &y);
		z = 0.0;
	    }
	  if (dimension_model == GAIA_XY_Z)
	      gaiaAddPointToGeomCollXYZ (gaia, x, y, z);
	  else if (dimension_model == GAIA_XY_M)
	      gaiaAddPointToGeomCollXYM (gaia, x, y, 0.0);
	  else if (dimension_model == GAIA_XY_Z_M)
	      gaiaAddPointToGeomCollXYZM (gaia, x, y, z, 0.0);
	  else
	      gaiaAddPointToGeomColl (gaia, x, y);
	  break;
      case GEOS_LINESTRING:
	  if (dimension_model == GAIA_XY_Z)
	      gaia = gaiaAllocGeomCollXYZ ();
	  else if (dimension_model == GAIA_XY_M)
	      gaia = gaiaAllocGeomCollXYM ();
	  else if (dimension_model == GAIA_XY_Z_M)
	      gaia = gaiaAllocGeomCollXYZM ();
	  else
	      gaia = gaiaAllocGeomColl ();
	  gaia->DeclaredType = GAIA_LINESTRING;
	  gaia->Srid = GEOSGetSRID (geos);
	  cs = GEOSGeom_getCoordSeq (geos);
	  GEOSCoordSeq_getDimensions (cs, &dims);
	  GEOSCoordSeq_getSize (cs, &points);
	  ln = gaiaAddLinestringToGeomColl (gaia, points);
	  for (iv = 0; iv < (int) points; iv++)
	    {
		if (dims == 3)
		  {
		      GEOSCoordSeq_getX (cs, iv, &x);
		      GEOSCoordSeq_getY (cs, iv, &y);
		      GEOSCoordSeq_getZ (cs, iv, &z);
		  }
		else
		  {
		      GEOSCoordSeq_getX (cs, iv, &x);
		      GEOSCoordSeq_getY (cs, iv, &y);
		      z = 0.0;
		  }
		if (dimension_model == GAIA_XY_Z)
		  {
		      gaiaSetPointXYZ (ln->Coords, iv, x, y, z);
		  }
		else if (dimension_model == GAIA_XY_M)
		  {
		      gaiaSetPointXYM (ln->Coords, iv, x, y, 0.0);
		  }
		else if (dimension_model == GAIA_XY_Z_M)
		  {
		      gaiaSetPointXYZM (ln->Coords, iv, x, y, z, 0.0);
		  }
		else
		  {
		      gaiaSetPoint (ln->Coords, iv, x, y);
		  }
	    }
	  break;
      case GEOS_POLYGON:
	  if (dimension_model == GAIA_XY_Z)
	      gaia = gaiaAllocGeomCollXYZ ();
	  else if (dimension_model == GAIA_XY_M)
	      gaia = gaiaAllocGeomCollXYM ();
	  else if (dimension_model == GAIA_XY_Z_M)
	      gaia = gaiaAllocGeomCollXYZM ();
	  else
	      gaia = gaiaAllocGeomColl ();
	  gaia->DeclaredType = GAIA_POLYGON;
	  gaia->Srid = GEOSGetSRID (geos);
	  /* exterior ring */
	  holes = GEOSGetNumInteriorRings (geos);
	  geos_ring = GEOSGetExteriorRing (geos);
	  cs = GEOSGeom_getCoordSeq (geos_ring);
	  GEOSCoordSeq_getDimensions (cs, &dims);
	  GEOSCoordSeq_getSize (cs, &points);
	  pg = gaiaAddPolygonToGeomColl (gaia, points, holes);
	  rng = pg->Exterior;
	  for (iv = 0; iv < (int) points; iv++)
	    {
		if (dims == 3)
		  {
		      GEOSCoordSeq_getX (cs, iv, &x);
		      GEOSCoordSeq_getY (cs, iv, &y);
		      GEOSCoordSeq_getZ (cs, iv, &z);
		  }
		else
		  {
		      GEOSCoordSeq_getX (cs, iv, &x);
		      GEOSCoordSeq_getY (cs, iv, &y);
		      z = 0.0;
		  }
		if (dimension_model == GAIA_XY_Z)
		  {
		      gaiaSetPointXYZ (rng->Coords, iv, x, y, z);
		  }
		else if (dimension_model == GAIA_XY_M)
		  {
		      gaiaSetPointXYM (rng->Coords, iv, x, y, 0.0);
		  }
		else if (dimension_model == GAIA_XY_Z_M)
		  {
		      gaiaSetPointXYZM (rng->Coords, iv, x, y, z, 0.0);
		  }
		else
		  {
		      gaiaSetPoint (rng->Coords, iv, x, y);
		  }
	    }
	  for (ib = 0; ib < holes; ib++)
	    {
		/* interior rings */
		geos_ring = GEOSGetInteriorRingN (geos, ib);
		cs = GEOSGeom_getCoordSeq (geos_ring);
		GEOSCoordSeq_getDimensions (cs, &dims);
		GEOSCoordSeq_getSize (cs, &points);
		rng = gaiaAddInteriorRing (pg, ib, points);
		for (iv = 0; iv < (int) points; iv++)
		  {
		      if (dims == 3)
			{
			    GEOSCoordSeq_getX (cs, iv, &x);
			    GEOSCoordSeq_getY (cs, iv, &y);
			    GEOSCoordSeq_getZ (cs, iv, &z);
			}
		      else
			{
			    GEOSCoordSeq_getX (cs, iv, &x);
			    GEOSCoordSeq_getY (cs, iv, &y);
			    z = 0.0;
			}
		      if (dimension_model == GAIA_XY_Z)
			{
			    gaiaSetPointXYZ (rng->Coords, iv, x, y, z);
			}
		      else if (dimension_model == GAIA_XY_M)
			{
			    gaiaSetPointXYM (rng->Coords, iv, x, y, 0.0);
			}
		      else if (dimension_model == GAIA_XY_Z_M)
			{
			    gaiaSetPointXYZM (rng->Coords, iv, x, y, z, 0.0);
			}
		      else
			{
			    gaiaSetPoint (rng->Coords, iv, x, y);
			}
		  }
	    }
	  break;
      case GEOS_MULTIPOINT:
      case GEOS_MULTILINESTRING:
      case GEOS_MULTIPOLYGON:
      case GEOS_GEOMETRYCOLLECTION:
	  if (dimension_model == GAIA_XY_Z)
	      gaia = gaiaAllocGeomCollXYZ ();
	  else if (dimension_model == GAIA_XY_M)
	      gaia = gaiaAllocGeomCollXYM ();
	  else if (dimension_model == GAIA_XY_Z_M)
	      gaia = gaiaAllocGeomCollXYZM ();
	  else
	      gaia = gaiaAllocGeomColl ();
	  if (type == GEOS_MULTIPOINT)
	      gaia->DeclaredType = GAIA_MULTIPOINT;
	  else if (type == GEOS_MULTILINESTRING)
	      gaia->DeclaredType = GAIA_MULTILINESTRING;
	  else if (type == GEOS_MULTIPOLYGON)
	      gaia->DeclaredType = GAIA_MULTIPOLYGON;
	  else
	      gaia->DeclaredType = GAIA_GEOMETRYCOLLECTION;
	  gaia->Srid = GEOSGetSRID (geos);
	  nItems = GEOSGetNumGeometries (geos);
	  for (it = 0; it < nItems; it++)
	    {
		/* looping on elementaty geometries */
		geos_item = GEOSGetGeometryN (geos, it);
		itemType = GEOSGeomTypeId (geos_item);
		switch (itemType)
		  {
		  case GEOS_POINT:
		      cs = GEOSGeom_getCoordSeq (geos_item);
		      GEOSCoordSeq_getDimensions (cs, &dims);
		      if (dims == 3)
			{
			    GEOSCoordSeq_getX (cs, 0, &x);
			    GEOSCoordSeq_getY (cs, 0, &y);
			    GEOSCoordSeq_getZ (cs, 0, &z);
			}
		      else
			{
			    GEOSCoordSeq_getX (cs, 0, &x);
			    GEOSCoordSeq_getY (cs, 0, &y);
			    z = 0.0;
			}
		      if (dimension_model == GAIA_XY_Z)
			  gaiaAddPointToGeomCollXYZ (gaia, x, y, z);
		      else if (dimension_model == GAIA_XY_M)
			  gaiaAddPointToGeomCollXYM (gaia, x, y, 0.0);
		      else if (dimension_model == GAIA_XY_Z_M)
			  gaiaAddPointToGeomCollXYZM (gaia, x, y, z, 0.0);
		      else
			  gaiaAddPointToGeomColl (gaia, x, y);
		      break;
		  case GEOS_LINESTRING:
		      cs = GEOSGeom_getCoordSeq (geos_item);
		      GEOSCoordSeq_getDimensions (cs, &dims);
		      GEOSCoordSeq_getSize (cs, &points);
		      ln = gaiaAddLinestringToGeomColl (gaia, points);
		      for (iv = 0; iv < (int) points; iv++)
			{
			    if (dims == 3)
			      {
				  GEOSCoordSeq_getX (cs, iv, &x);
				  GEOSCoordSeq_getY (cs, iv, &y);
				  GEOSCoordSeq_getZ (cs, iv, &z);
			      }
			    else
			      {
				  GEOSCoordSeq_getX (cs, iv, &x);
				  GEOSCoordSeq_getY (cs, iv, &y);
				  z = 0.0;
			      }
			    if (dimension_model == GAIA_XY_Z)
			      {
				  gaiaSetPointXYZ (ln->Coords, iv, x, y, z);
			      }
			    else if (dimension_model == GAIA_XY_M)
			      {
				  gaiaSetPointXYM (ln->Coords, iv, x, y, 0.0);
			      }
			    else if (dimension_model == GAIA_XY_Z_M)
			      {
				  gaiaSetPointXYZM (ln->Coords, iv, x, y, z,
						    0.0);
			      }
			    else
			      {
				  gaiaSetPoint (ln->Coords, iv, x, y);
			      }
			}
		      break;
		  case GEOS_MULTILINESTRING:
		      nSubItems = GEOSGetNumGeometries (geos_item);
		      for (sub_it = 0; sub_it < nSubItems; sub_it++)
			{
			    /* looping on elementaty geometries */
			    geos_sub_item =
				GEOSGetGeometryN (geos_item, sub_it);
			    cs = GEOSGeom_getCoordSeq (geos_sub_item);
			    GEOSCoordSeq_getDimensions (cs, &dims);
			    GEOSCoordSeq_getSize (cs, &points);
			    ln = gaiaAddLinestringToGeomColl (gaia, points);
			    for (iv = 0; iv < (int) points; iv++)
			      {
				  if (dims == 3)
				    {
					GEOSCoordSeq_getX (cs, iv, &x);
					GEOSCoordSeq_getY (cs, iv, &y);
					GEOSCoordSeq_getZ (cs, iv, &z);
				    }
				  else
				    {
					GEOSCoordSeq_getX (cs, iv, &x);
					GEOSCoordSeq_getY (cs, iv, &y);
					z = 0.0;
				    }
				  if (dimension_model == GAIA_XY_Z)
				    {
					gaiaSetPointXYZ (ln->Coords, iv, x, y,
							 z);
				    }
				  else if (dimension_model == GAIA_XY_M)
				    {
					gaiaSetPointXYM (ln->Coords, iv, x, y,
							 0.0);
				    }
				  else if (dimension_model == GAIA_XY_Z_M)
				    {
					gaiaSetPointXYZM (ln->Coords, iv, x, y,
							  z, 0.0);
				    }
				  else
				    {
					gaiaSetPoint (ln->Coords, iv, x, y);
				    }
			      }
			}
		      break;
		  case GEOS_POLYGON:
		      /* exterior ring */
		      holes = GEOSGetNumInteriorRings (geos_item);
		      geos_ring = GEOSGetExteriorRing (geos_item);
		      cs = GEOSGeom_getCoordSeq (geos_ring);
		      GEOSCoordSeq_getDimensions (cs, &dims);
		      GEOSCoordSeq_getSize (cs, &points);
		      pg = gaiaAddPolygonToGeomColl (gaia, points, holes);
		      rng = pg->Exterior;
		      for (iv = 0; iv < (int) points; iv++)
			{
			    if (dims == 3)
			      {
				  GEOSCoordSeq_getX (cs, iv, &x);
				  GEOSCoordSeq_getY (cs, iv, &y);
				  GEOSCoordSeq_getZ (cs, iv, &z);
			      }
			    else
			      {
				  GEOSCoordSeq_getX (cs, iv, &x);
				  GEOSCoordSeq_getY (cs, iv, &y);
				  z = 0.0;
			      }
			    if (dimension_model == GAIA_XY_Z)
			      {
				  gaiaSetPointXYZ (rng->Coords, iv, x, y, z);
			      }
			    else if (dimension_model == GAIA_XY_M)
			      {
				  gaiaSetPointXYM (rng->Coords, iv, x, y, 0.0);
			      }
			    else if (dimension_model == GAIA_XY_Z_M)
			      {
				  gaiaSetPointXYZM (rng->Coords, iv, x, y, z,
						    0.0);
			      }
			    else
			      {
				  gaiaSetPoint (rng->Coords, iv, x, y);
			      }
			}
		      for (ib = 0; ib < holes; ib++)
			{
			    /* interior rings */
			    geos_ring = GEOSGetInteriorRingN (geos_item, ib);
			    cs = GEOSGeom_getCoordSeq (geos_ring);
			    GEOSCoordSeq_getDimensions (cs, &dims);
			    GEOSCoordSeq_getSize (cs, &points);
			    rng = gaiaAddInteriorRing (pg, ib, points);
			    for (iv = 0; iv < (int) points; iv++)
			      {
				  if (dims == 3)
				    {
					GEOSCoordSeq_getX (cs, iv, &x);
					GEOSCoordSeq_getY (cs, iv, &y);
					GEOSCoordSeq_getZ (cs, iv, &z);
				    }
				  else
				    {
					GEOSCoordSeq_getX (cs, iv, &x);
					GEOSCoordSeq_getY (cs, iv, &y);
					z = 0.0;
				    }
				  if (dimension_model == GAIA_XY_Z)
				    {
					gaiaSetPointXYZ (rng->Coords, iv, x, y,
							 z);
				    }
				  else if (dimension_model == GAIA_XY_M)
				    {
					gaiaSetPointXYM (rng->Coords, iv, x, y,
							 0.0);
				    }
				  else if (dimension_model == GAIA_XY_Z_M)
				    {
					gaiaSetPointXYZM (rng->Coords, iv, x, y,
							  z, 0.0);
				    }
				  else
				    {
					gaiaSetPoint (rng->Coords, iv, x, y);
				    }
			      }
			}
		      break;
		  };
	    }
	  break;
      };
    return gaia;
}

GAIAGEO_DECLARE void *
gaiaToGeos (const gaiaGeomCollPtr gaia)
{
/* converting a GAIA Geometry into a GEOS Geometry */
    return toGeosGeometry (gaia);
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaFromGeos_XY (const void *xgeos)
{
/* converting a GEOS Geometry into a GAIA Geometry [XY] */
    const GEOSGeometry *geos = xgeos;
    return fromGeosGeometry (geos, GAIA_XY);
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaFromGeos_XYZ (const void *xgeos)
{
/* converting a GEOS Geometry into a GAIA Geometry [XYZ] */
    const GEOSGeometry *geos = xgeos;
    return fromGeosGeometry (geos, GAIA_XY_Z);
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaFromGeos_XYM (const void *xgeos)
{
/* converting a GEOS Geometry into a GAIA Geometry [XYM] */
    const GEOSGeometry *geos = xgeos;
    return fromGeosGeometry (geos, GAIA_XY_M);
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaFromGeos_XYZM (const void *xgeos)
{
/* converting a GEOS Geometry into a GAIA Geometry [XYZM] */
    const GEOSGeometry *geos = xgeos;
    return fromGeosGeometry (geos, GAIA_XY_Z_M);
}

#endif /* end including GEOS */
