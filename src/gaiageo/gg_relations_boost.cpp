/*

 gg_relations.c -- Gaia spatial relations
    
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
#include <string.h>
#include <float.h>

#ifndef OMIT_BOOSTGEOMETRY		/* including GEOS */
#include "gg_geoscvt_boost.h"
#include <boost/geometry/geometry.hpp>
#endif

#ifdef SPL_AMALGAMATION		/* spatialite-amalgamation */
#include <spatialite/sqlite3ext.h>
#else
#include <sqlite3ext.h>
#endif



#include <spatialite/gaiageo.h>


#ifndef OMIT_BOOSTGEOMETRY		/* including GEOS */

namespace bg=boost::geometry;

GAIAGEO_DECLARE int
gaiaGeomCollEquals (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
	if (!geom1 || !geom2)
	return -1;

	int i;
	geometry_collection g1 = toBoostGeometry(geom1);
	geometry_collection g2 = toBoostGeometry(geom2);

	//multilinestrings and multipoints are not supported by boost right now
	if(g1.points.size() != g2.points.size()) return 0;
	for(i=0; i< g1.points.size(); i++) {
		if(!bg::equals(g1.points[i], g2.points[i])) return 0;
	}

	if(g1.linestrings.size() != g2.linestrings.size()) return 0;
	for(i=0; i< g1.linestrings.size(); i++) {
		if(!bg::equals(g1.linestrings[i], g2.linestrings[i])) return 0;
	}

	if(!bg::equals(g1.polygons, g2.polygons)) return 0;

	return 1;
}



GAIAGEO_DECLARE int
gaiaGeomCollIntersects (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if two Geometries do "spatially intersects" */
	if (!geom1 || !geom2)
	return -1;

	int i,j;
	geometry_collection g1 = toBoostGeometry(geom1);
	geometry_collection g2 = toBoostGeometry(geom2);

	for(i=0; i< g1.linestrings.size(); i++) {
		for(j=0; j< g2.linestrings.size(); j++) {
			if(!bg::intersects(g1.linestrings[i], g2.linestrings[j])) return 0;
		}
	}

	for(i=0; i< g1.polygons.size(); i++) {
		for(j=0; j< g2.polygons.size(); j++) {
			if(!bg::intersects(g1.polygons[i], g2.polygons[j])) return 0;
		}
	}

	for(i=0; i< g1.linestrings.size(); i++) {
		for(j=0; j< g2.polygons.size(); j++) {
			if(!bg::intersects(g1.linestrings[i], g2.polygons[j])) return 0;
		}
	}

	for(i=0; i< g1.polygons.size(); i++) {
		for(j=0; j< g2.linestrings.size(); j++) {
			if(!bg::intersects(g1.polygons[i], g2.linestrings[j])) return 0;
		}
	}
}

#if 0

GAIAGEO_DECLARE int
gaiaGeomCollDisjoint (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if two Geometries are "spatially disjoint" */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSDisjoint (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollOverlaps (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if two Geometries do "spatially overlaps" */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSOverlaps (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollCrosses (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if two Geometries do "spatially crosses" */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSCrosses (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollTouches (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if two Geometries do "spatially touches" */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSTouches (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollWithin (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if GEOM-1 is completely contained within GEOM-2 */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSWithin (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollContains (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* checks if GEOM-1 completely contains GEOM-2 */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSContains (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollRelate (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2,
		    const char *pattern)
{
/* checks if if GEOM-1 and GEOM-2 have a spatial relationship as specified by the pattern Matrix */
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return -1;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSRelatePattern (g1, g2, pattern);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollLength (gaiaGeomCollPtr geom, double *xlength)
{
/* computes the total length for this Geometry */
    double length;
    int ret;
    GEOSGeometry *g;
    if (!geom)
	return 0;
    g = gaiaToGeos (geom);
    ret = GEOSLength (g, &length);
    GEOSGeom_destroy (g);
    if (ret)
	*xlength = length;
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollArea (gaiaGeomCollPtr geom, double *xarea)
{
/* computes the total area for this Geometry */
    double area;
    int ret;
    GEOSGeometry *g;
    if (!geom)
	return 0;
    g = gaiaToGeos (geom);
    ret = GEOSArea (g, &area);
    GEOSGeom_destroy (g);
    if (ret)
	*xarea = area;
    return ret;
}

GAIAGEO_DECLARE int
gaiaGeomCollDistance (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2,
		      double *xdist)
{
/* computes the minimum distance intercurring between GEOM-1 and GEOM-2 */
    double dist;
    int ret;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom1 || !geom2)
	return 0;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    ret = GEOSDistance (g1, g2, &dist);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    if (ret)
	*xdist = dist;
    return ret;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeometryIntersection (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* builds a new geometry representing the "spatial intersection" of GEOM-1 and GEOM-2 */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    GEOSGeometry *g3;
    if (!geom1 || !geom2)
	return NULL;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    g3 = GEOSIntersection (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    if (!g3)
	return NULL;
    if (geom1->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g3);
    else if (geom1->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g3);
    else if (geom1->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g3);
    else
	geo = gaiaFromGeos_XY (g3);
    GEOSGeom_destroy (g3);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom1->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeometryUnion (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* builds a new geometry representing the "spatial union" of GEOM-1 and GEOM-2 */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    GEOSGeometry *g3;
    if (!geom1 || !geom2)
	return NULL;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    g3 = GEOSUnion (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    if (geom1->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g3);
    else if (geom1->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g3);
    else if (geom1->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g3);
    else
	geo = gaiaFromGeos_XY (g3);
    GEOSGeom_destroy (g3);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom1->Srid;
    if (geo->DeclaredType == GAIA_POINT &&
	geom1->DeclaredType == GAIA_MULTIPOINT)
	geo->DeclaredType = GAIA_MULTIPOINT;
    if (geo->DeclaredType == GAIA_LINESTRING &&
	geom1->DeclaredType == GAIA_MULTILINESTRING)
	geo->DeclaredType = GAIA_MULTILINESTRING;
    if (geo->DeclaredType == GAIA_POLYGON &&
	geom1->DeclaredType == GAIA_MULTIPOLYGON)
	geo->DeclaredType = GAIA_MULTIPOLYGON;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeometryDifference (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* builds a new geometry representing the "spatial difference" of GEOM-1 and GEOM-2 */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    GEOSGeometry *g3;
    if (!geom1 || !geom2)
	return NULL;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    g3 = GEOSDifference (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    if (!g3)
	return NULL;
    if (geom1->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g3);
    else if (geom1->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g3);
    else if (geom1->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g3);
    else
	geo = gaiaFromGeos_XY (g3);
    GEOSGeom_destroy (g3);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom1->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeometrySymDifference (gaiaGeomCollPtr geom1, gaiaGeomCollPtr geom2)
{
/* builds a new geometry representing the "spatial symmetric difference" of GEOM-1 and GEOM-2 */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    GEOSGeometry *g3;
    if (!geom1 || !geom2)
	return NULL;
    g1 = gaiaToGeos (geom1);
    g2 = gaiaToGeos (geom2);
    g3 = GEOSSymDifference (g1, g2);
    GEOSGeom_destroy (g1);
    GEOSGeom_destroy (g2);
    if (!g3)
	return NULL;
    if (geom1->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g3);
    else if (geom1->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g3);
    else if (geom1->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g3);
    else
	geo = gaiaFromGeos_XY (g3);
    GEOSGeom_destroy (g3);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom1->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaBoundary (gaiaGeomCollPtr geom)
{
/* builds a new geometry representing the conbinatorial boundary of GEOM */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return NULL;
    g1 = gaiaToGeos (geom);
    g2 = GEOSBoundary (g1);
    GEOSGeom_destroy (g1);
    if (!g2)
	return NULL;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom->Srid;
    return geo;
}

GAIAGEO_DECLARE int
gaiaGeomCollCentroid (gaiaGeomCollPtr geom, double *x, double *y)
{
/* returns a Point representing the centroid for this Geometry */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return 0;
    g1 = gaiaToGeos (geom);
    g2 = GEOSGetCentroid (g1);
    GEOSGeom_destroy (g1);
    if (!g2)
	return 0;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return 0;
    if (geo->FirstPoint)
      {
	  *x = geo->FirstPoint->X;
	  *y = geo->FirstPoint->Y;
	  gaiaFreeGeomColl (geo);
	  return 1;
      }
    gaiaFreeGeomColl (geo);
    return 0;
}

GAIAGEO_DECLARE int
gaiaGetPointOnSurface (gaiaGeomCollPtr geom, double *x, double *y)
{
/* returns a Point guaranteed to lie on the Surface */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return 0;
    g1 = gaiaToGeos (geom);
    g2 = GEOSPointOnSurface (g1);
    GEOSGeom_destroy (g1);
    if (!g2)
	return 0;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return 0;
    if (geo->FirstPoint)
      {
	  *x = geo->FirstPoint->X;
	  *y = geo->FirstPoint->Y;
	  gaiaFreeGeomColl (geo);
	  return 1;
      }
    gaiaFreeGeomColl (geo);
    return 0;
}

GAIAGEO_DECLARE int
gaiaIsSimple (gaiaGeomCollPtr geom)
{
/* checks if this GEOMETRYCOLLECTION is a simple one */
    int ret;
    GEOSGeometry *g;
    if (!geom)
	return -1;
    if (gaiaIsToxic (geom))
	return 0;
    g = gaiaToGeos (geom);
    ret = GEOSisSimple (g);
    GEOSGeom_destroy (g);
    if (ret == 2)
	return -1;
    return ret;
}

GAIAGEO_DECLARE int
gaiaIsRing (gaiaLinestringPtr line)
{
/* checks if this LINESTRING can be a valid RING */
    gaiaGeomCollPtr geo;
    gaiaLinestringPtr line2;
    int ret;
    int iv;
    double x;
    double y;
    double z;
    double m;
    GEOSGeometry *g;
    if (!line)
	return -1;
    if (line->DimensionModel == GAIA_XY_Z)
	geo = gaiaAllocGeomCollXYZ ();
    else if (line->DimensionModel == GAIA_XY_M)
	geo = gaiaAllocGeomCollXYM ();
    else if (line->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaAllocGeomCollXYZM ();
    else
	geo = gaiaAllocGeomColl ();
    line2 = gaiaAddLinestringToGeomColl (geo, line->Points);
    for (iv = 0; iv < line2->Points; iv++)
      {
	  z = 0.0;
	  m = 0.0;
	  if (line->DimensionModel == GAIA_XY_Z)
	    {
		gaiaGetPointXYZ (line->Coords, iv, &x, &y, &z);
	    }
	  else if (line->DimensionModel == GAIA_XY_M)
	    {
		gaiaGetPointXYM (line->Coords, iv, &x, &y, &m);
	    }
	  else if (line->DimensionModel == GAIA_XY_Z_M)
	    {
		gaiaGetPointXYZM (line->Coords, iv, &x, &y, &z, &m);
	    }
	  else
	    {
		gaiaGetPoint (line->Coords, iv, &x, &y);
	    }
	  if (line2->DimensionModel == GAIA_XY_Z)
	    {
		gaiaSetPointXYZ (line2->Coords, iv, x, y, z);
	    }
	  else if (line2->DimensionModel == GAIA_XY_M)
	    {
		gaiaSetPointXYM (line2->Coords, iv, x, y, m);
	    }
	  else if (line2->DimensionModel == GAIA_XY_Z_M)
	    {
		gaiaSetPointXYZM (line2->Coords, iv, x, y, z, m);
	    }
	  else
	    {
		gaiaSetPoint (line2->Coords, iv, x, y);
	    }
      }
    g = gaiaToGeos (geo);
    gaiaFreeGeomColl (geo);
    ret = GEOSisRing (g);
    GEOSGeom_destroy (g);
    if (ret == 2)
	return -1;
    return ret;
}

GAIAGEO_DECLARE int
gaiaIsValid (gaiaGeomCollPtr geom)
{
/* checks if this GEOMETRYCOLLECTION is a valid one */
    int ret;
    GEOSGeometry *g;
    gaiaResetGeosMsg ();
    if (!geom)
	return -1;
    if (gaiaIsToxic (geom))
	return 0;
    g = gaiaToGeos (geom);
    ret = GEOSisValid (g);
    GEOSGeom_destroy (g);
    if (ret == 2)
	return -1;
    return ret;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeomCollSimplify (gaiaGeomCollPtr geom, double tolerance)
{
/* builds a simplified geometry using the Douglas-Peuker algorihtm */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return NULL;
    g1 = gaiaToGeos (geom);
    g2 = GEOSSimplify (g1, tolerance);
    GEOSGeom_destroy (g1);
    if (!g2)
	return NULL;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeomCollSimplifyPreserveTopology (gaiaGeomCollPtr geom, double tolerance)
{
/* builds a simplified geometry using the Douglas-Peuker algorihtm [preserving topology] */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return NULL;
    g1 = gaiaToGeos (geom);
    g2 = GEOSTopologyPreserveSimplify (g1, tolerance);
    GEOSGeom_destroy (g1);
    if (!g2)
	return NULL;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaConvexHull (gaiaGeomCollPtr geom)
{
/* builds a geometry that is the convex hull of GEOM */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return NULL;
    g1 = gaiaToGeos (geom);
    g2 = GEOSConvexHull (g1);
    GEOSGeom_destroy (g1);
    if (!g2)
	return NULL;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom->Srid;
    return geo;
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaGeomCollBuffer (gaiaGeomCollPtr geom, double radius, int points)
{
/* builds a geometry that is the GIS buffer of GEOM */
    gaiaGeomCollPtr geo;
    GEOSGeometry *g1;
    GEOSGeometry *g2;
    if (!geom)
	return NULL;
    if (gaiaIsToxic (geom))
	return NULL;
    g1 = gaiaToGeos (geom);
    g2 = GEOSBuffer (g1, radius, points);
    GEOSGeom_destroy (g1);
    if (!g2)
	return NULL;
    if (geom->DimensionModel == GAIA_XY_Z)
	geo = gaiaFromGeos_XYZ (g2);
    else if (geom->DimensionModel == GAIA_XY_M)
	geo = gaiaFromGeos_XYM (g2);
    else if (geom->DimensionModel == GAIA_XY_Z_M)
	geo = gaiaFromGeos_XYZM (g2);
    else
	geo = gaiaFromGeos_XY (g2);
    GEOSGeom_destroy (g2);
    if (geo == NULL)
	return NULL;
    geo->Srid = geom->Srid;
    return geo;
}

static void
test_interior_ring (gaiaDynamicLinePtr dyn1, gaiaDynamicLinePtr dyn2,
		    int *contains, int *within, int *crosses)
{
/* testing if Ring-1 contains Ring-2 */
    gaiaGeomCollPtr geom1;
    gaiaGeomCollPtr geom2;
    gaiaPolygonPtr pg;
    gaiaRingPtr rng;
    int iv;
    int pts;
    gaiaPointPtr pt;

/* creating the Polygon-1 geometry */
    pts = 0;
    pt = dyn1->First;
    while (pt)
      {
	  pts++;
	  pt = pt->Next;
      }
    geom1 = gaiaAllocGeomColl ();
    pg = gaiaAddPolygonToGeomColl (geom1, pts, 0);
    rng = pg->Exterior;
    iv = 0;
    pt = dyn1->First;
    while (pt)
      {
	  /* EXTERIOR RING */
	  gaiaSetPoint (rng->Coords, iv, pt->X, pt->Y);
	  iv++;
	  pt = pt->Next;
      }

/* creating the Polygon-2 geometry */
    pts = 0;
    pt = dyn2->First;
    while (pt)
      {
	  pts++;
	  pt = pt->Next;
      }
    geom2 = gaiaAllocGeomColl ();
    pg = gaiaAddPolygonToGeomColl (geom2, pts, 0);
    rng = pg->Exterior;
    iv = 0;
    pt = dyn2->First;
    while (pt)
      {
	  /* EXTERIOR RING */
	  gaiaSetPoint (rng->Coords, iv, pt->X, pt->Y);
	  iv++;
	  pt = pt->Next;
      }
    *contains = gaiaGeomCollContains (geom1, geom2);
    *within = gaiaGeomCollWithin (geom1, geom2);
    *crosses = gaiaGeomCollCrosses (geom1, geom2);
    gaiaFreeGeomColl (geom1);
    gaiaFreeGeomColl (geom2);
}

static gaiaDynamicLinePtr
build_dyn_ring (gaiaLinestringPtr ln)
{
/* creating a DynamicLine from a Linestring */
    int iv;
    double x;
    double y;
    double m;
    double z;
    gaiaDynamicLinePtr dyn = gaiaAllocDynamicLine ();
    for (iv = 0; iv < ln->Points; iv++)
      {
	  if (ln->DimensionModel == GAIA_XY_Z_M)
	    {
		gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		gaiaAppendPointZMToDynamicLine (dyn, x, y, z, m);
	    }
	  else if (ln->DimensionModel == GAIA_XY_Z)
	    {
		gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		gaiaAppendPointZToDynamicLine (dyn, x, y, z);
	    }
	  else if (ln->DimensionModel == GAIA_XY_M)
	    {
		gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		gaiaAppendPointMToDynamicLine (dyn, x, y, m);
	    }
	  else
	    {
		gaiaGetPoint (ln->Coords, iv, &x, &y);
		gaiaAppendPointToDynamicLine (dyn, x, y);
	    }
      }
    return dyn;
}

static int
is_closed_dyn_ring (gaiaDynamicLinePtr dyn)
{
/* checking if a candidate Ring is already closed */
    gaiaPointPtr pt1;
    gaiaPointPtr pt2;
    if (!dyn)
	return 0;
    pt1 = dyn->First;
    pt2 = dyn->Last;
    if (pt1 == NULL || pt2 == NULL)
	return 0;
    if (pt1 == pt2)
	return 0;
    if (pt1->X == pt2->X && pt1->Y == pt2->Y && pt1->Z == pt2->Z)
	return 1;
    return 0;
}

static int
to_be_appended (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln)
{
/* checks is the Linestring has to be appended to the DynamicLine */
    gaiaPointPtr pt = dyn->Last;
    int iv = 0;
    double x;
    double y;
    double z;
    double m;
    if (ln->DimensionModel == GAIA_XY_Z_M)
      {
	  gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
      }
    else if (ln->DimensionModel == GAIA_XY_Z)
      {
	  gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
      }
    else if (ln->DimensionModel == GAIA_XY_M)
      {
	  gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
      }
    else
      {
	  gaiaGetPoint (ln->Coords, iv, &x, &y);
      }
    if (ln->DimensionModel == GAIA_XY_Z_M || ln->DimensionModel == GAIA_XY_Z)
      {
	  if (pt->X == x && pt->Y == y && pt->Z == z)
	      return 1;
      }
    else
      {
	  if (pt->X == x && pt->Y == y)
	      return 1;
      }
    return 0;
}

static int
to_be_prepended (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln)
{
/* checks is the Linestring has to be prepended to the DynamicLine */
    gaiaPointPtr pt = dyn->First;
    int iv = ln->Points - 1;
    double x;
    double y;
    double z;
    double m;
    if (ln->DimensionModel == GAIA_XY_Z_M)
      {
	  gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
      }
    else if (ln->DimensionModel == GAIA_XY_Z)
      {
	  gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
      }
    else if (ln->DimensionModel == GAIA_XY_M)
      {
	  gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
      }
    else
      {
	  gaiaGetPoint (ln->Coords, iv, &x, &y);
      }
    if (ln->DimensionModel == GAIA_XY_Z_M || ln->DimensionModel == GAIA_XY_Z)
      {
	  if (pt->X == x && pt->Y == y && pt->Z == z)
	      return 1;
      }
    else
      {
	  if (pt->X == x && pt->Y == y)
	      return 1;
      }
    return 0;
}

static int
to_be_appended_reverse (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln)
{
/* checks is the Linestring (reversed) has to be appended to the DynamicLine */
    gaiaPointPtr pt = dyn->Last;
    int iv = ln->Points - 1;
    double x;
    double y;
    double z;
    double m;
    if (ln->DimensionModel == GAIA_XY_Z_M)
      {
	  gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
      }
    else if (ln->DimensionModel == GAIA_XY_Z)
      {
	  gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
      }
    else if (ln->DimensionModel == GAIA_XY_M)
      {
	  gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
      }
    else
      {
	  gaiaGetPoint (ln->Coords, iv, &x, &y);
      }
    if (ln->DimensionModel == GAIA_XY_Z_M || ln->DimensionModel == GAIA_XY_Z)
      {
	  if (pt->X == x && pt->Y == y && pt->Z == z)
	      return 1;
      }
    else
      {
	  if (pt->X == x && pt->Y == y)
	      return 1;
      }
    return 0;
}

static int
to_be_prepended_reverse (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln)
{
/* checks is the Linestring (reversed) has to be prepended to the DynamicLine */
    gaiaPointPtr pt = dyn->First;
    int iv = 0;
    double x;
    double y;
    double z;
    double m;
    if (ln->DimensionModel == GAIA_XY_Z_M)
      {
	  gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
      }
    else if (ln->DimensionModel == GAIA_XY_Z)
      {
	  gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
      }
    else if (ln->DimensionModel == GAIA_XY_M)
      {
	  gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
      }
    else
      {
	  gaiaGetPoint (ln->Coords, iv, &x, &y);
      }
    if (ln->DimensionModel == GAIA_XY_Z_M || ln->DimensionModel == GAIA_XY_Z)
      {
	  if (pt->X == x && pt->Y == y && pt->Z == z)
	      return 1;
      }
    else
      {
	  if (pt->X == x && pt->Y == y)
	      return 1;
      }
    return 0;
}

static void
append_to_ring (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln, int reversed)
{
/* appending a Linestring to a DynamicRing */
    int iv;
    double x;
    double y;
    double z;
    double m;
    if (reversed)
      {
	  for (iv = ln->Points - 2; iv >= 0; iv--)
	    {
		if (ln->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		      gaiaAppendPointZMToDynamicLine (dyn, x, y, z, m);
		  }
		else if (ln->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		      gaiaAppendPointZToDynamicLine (dyn, x, y, z);
		  }
		else if (ln->DimensionModel == GAIA_XY_M)
		  {
		      gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		      gaiaAppendPointMToDynamicLine (dyn, x, y, m);
		  }
		else
		  {
		      gaiaGetPoint (ln->Coords, iv, &x, &y);
		      gaiaAppendPointToDynamicLine (dyn, x, y);
		  }
	    }
      }
    else
      {
	  for (iv = 1; iv < ln->Points; iv++)
	    {
		if (ln->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		      gaiaAppendPointZMToDynamicLine (dyn, x, y, z, m);
		  }
		else if (ln->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		      gaiaAppendPointZToDynamicLine (dyn, x, y, z);
		  }
		else if (ln->DimensionModel == GAIA_XY_M)
		  {
		      gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		      gaiaAppendPointMToDynamicLine (dyn, x, y, m);
		  }
		else
		  {
		      gaiaGetPoint (ln->Coords, iv, &x, &y);
		      gaiaAppendPointToDynamicLine (dyn, x, y);
		  }
	    }
      }
}

static void
prepend_to_ring (gaiaDynamicLinePtr dyn, gaiaLinestringPtr ln, int reversed)
{
/* appending a Linestring to a DynamicRing */
    int iv;
    double x;
    double y;
    double z;
    double m;
    if (reversed)
      {
	  for (iv = 1; iv < ln->Points; iv++)
	    {
		if (ln->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		      gaiaPrependPointZMToDynamicLine (dyn, x, y, z, m);
		  }
		else if (ln->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		      gaiaPrependPointZToDynamicLine (dyn, x, y, z);
		  }
		else if (ln->DimensionModel == GAIA_XY_M)
		  {
		      gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		      gaiaPrependPointMToDynamicLine (dyn, x, y, m);
		  }
		else
		  {
		      gaiaGetPoint (ln->Coords, iv, &x, &y);
		      gaiaPrependPointToDynamicLine (dyn, x, y);
		  }
	    }
      }
    else
      {
	  for (iv = ln->Points - 2; iv >= 0; iv--)
	    {
		if (ln->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		      gaiaPrependPointZMToDynamicLine (dyn, x, y, z, m);
		  }
		else if (ln->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		      gaiaPrependPointZToDynamicLine (dyn, x, y, z);
		  }
		else if (ln->DimensionModel == GAIA_XY_M)
		  {
		      gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		      gaiaPrependPointMToDynamicLine (dyn, x, y, m);
		  }
		else
		  {
		      gaiaGetPoint (ln->Coords, iv, &x, &y);
		      gaiaPrependPointToDynamicLine (dyn, x, y);
		  }
	    }
      }
}

GAIAGEO_DECLARE gaiaGeomCollPtr
gaiaPolygonize (gaiaGeomCollPtr geom, int force_multi)
{
/* attempts to rearrange a generic Geometry into a (multi)polygon */
    int pts = 0;
    int lns = 0;
    int pgs = 0;
    gaiaGeomCollPtr result;
    gaiaPointPtr pt;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;
    gaiaRingPtr rng;
    int dummy;
    int ok;
    int ok2;
    int i;
    int i2;
    int iv;
    int ib;
    double x;
    double y;
    double m;
    double z;
    double x0;
    double y0;
    double z0;
    int contains;
    int within;
    int crosses;
    int num_interiors;
    gaiaLinestringPtr *ln_array = NULL;
    gaiaDynamicLinePtr *dyn_array = NULL;
    gaiaDynamicLinePtr *ext_array = NULL;
    gaiaDynamicLinePtr dyn;
    gaiaDynamicLinePtr dyn2;

    if (!geom)
	return NULL;
    pt = geom->FirstPoint;
    while (pt)
      {
	  pts++;
	  pt = pt->Next;
      }
    pg = geom->FirstPolygon;
    while (pg)
      {
	  pgs++;
	  pg = pg->Next;
      }
    if (pts || pgs)
	return NULL;
    ln = geom->FirstLinestring;
    while (ln)
      {
	  lns++;
	  ln = ln->Next;
      }
    if (!lns)
	return NULL;
/* allocating and initializing aux-arrays */
    ln_array = malloc (sizeof (gaiaLinestringPtr) * lns);
    dyn_array = malloc (sizeof (gaiaDynamicLinePtr) * lns);
    ext_array = malloc (sizeof (gaiaDynamicLinePtr) * lns);
    i = 0;
    ln = geom->FirstLinestring;
    while (ln)
      {
	  ln_array[i] = ln;
	  dyn_array[i] = NULL;
	  ext_array[i] = NULL;
	  i++;
	  ln = ln->Next;
      }

    for (i = 0; i < lns; i++)
      {
	  /* processing closed rings */
	  ln = ln_array[i];
	  iv = ln->Points - 1;
	  if (ln->DimensionModel == GAIA_XY_Z_M)
	    {
		gaiaGetPointXYZM (ln->Coords, 0, &x0, &y0, &z0, &m);
		gaiaGetPointXYZM (ln->Coords, iv, &x, &y, &z, &m);
		if (x0 == x && y0 == y && z0 == z)
		  {
		      dyn_array[i] = build_dyn_ring (ln);
		      ln_array[i] = NULL;
		  }
	    }
	  else if (ln->DimensionModel == GAIA_XY_Z)
	    {
		gaiaGetPointXYZ (ln->Coords, 0, &x0, &y0, &z0);
		gaiaGetPointXYZ (ln->Coords, iv, &x, &y, &z);
		if (x0 == x && y0 == y && z0 == z)
		  {
		      dyn_array[i] = build_dyn_ring (ln);
		      ln_array[i] = NULL;
		  }
	    }
	  else if (ln->DimensionModel == GAIA_XY_M)
	    {
		gaiaGetPointXYM (ln->Coords, 0, &x0, &y0, &m);
		gaiaGetPointXYM (ln->Coords, iv, &x, &y, &m);
		if (x0 == x && y0 == y)
		  {
		      dyn_array[i] = build_dyn_ring (ln);
		      ln_array[i] = NULL;
		  }
	    }
	  else
	    {
		gaiaGetPoint (ln->Coords, 0, &x0, &y0);
		gaiaGetPoint (ln->Coords, iv, &x, &y);
		if (x0 == x && y0 == y)
		  {
		      dyn_array[i] = build_dyn_ring (ln);
		      ln_array[i] = NULL;
		  }
	    }
      }

    ok = 1;
    while (ok)
      {
	  if (dummy == 0)
	      ok = dummy;	/* simply suppressing stupid compiler warnings */
	  ok = 0;
	  for (i = 0; i < lns; i++)
	    {
		/* attempting to create rings */
		ln = ln_array[i];
		if (ln == NULL)
		    continue;
		ok = 1;
		dyn_array[i] = build_dyn_ring (ln);
		ln_array[i] = NULL;
		dyn = dyn_array[i];
		ok2 = 1;
		while (ok2)
		  {
		      ok2 = 0;
		      for (i2 = 0; i2 < lns; i2++)
			{
			    if (is_closed_dyn_ring (dyn) == 1)
				goto ring_done;
			    ln = ln_array[i2];
			    if (ln == NULL)
				continue;
			    if (to_be_appended (dyn, ln) == 1)
			      {
				  append_to_ring (dyn, ln, 0);
				  ln_array[i2] = NULL;
				  ok2 = 1;
				  break;
			      }
			    if (to_be_prepended (dyn, ln) == 1)
			      {
				  prepend_to_ring (dyn, ln, 0);
				  ln_array[i2] = NULL;
				  ok2 = 1;
				  break;
			      }
			    if (to_be_appended_reverse (dyn, ln) == 1)
			      {
				  append_to_ring (dyn, ln, 1);
				  ln_array[i2] = NULL;
				  ok2 = 1;
				  break;
			      }
			    if (to_be_prepended_reverse (dyn, ln) == 1)
			      {
				  prepend_to_ring (dyn, ln, 1);
				  ln_array[i2] = NULL;
				  ok2 = 1;
				  break;
			      }
			}
		  }
	    }
	ring_done:
	  dummy = 0;
      }

    ok = 1;
    for (i = 0; i < lns; i++)
      {
	  /* checking if any ring is closed */
	  dyn = dyn_array[i];
	  if (dyn == NULL)
	      continue;
	  if (is_closed_dyn_ring (dyn) == 0)
	      ok = 0;
      }
    if (ok == 0)
      {
	  /* invalid: quitting */
	  for (i = 0; i < lns; i++)
	    {
		dyn = dyn_array[i];
		if (dyn == NULL)
		    continue;
		gaiaFreeDynamicLine (dyn);
	    }
	  free (dyn_array);
	  free (ext_array);
	  free (ln_array);
	  return NULL;
      }

    ok = 1;
    for (i = 0; i < lns; i++)
      {
	  /* testing interior/exterior relationships */
	  dyn = dyn_array[i];
	  if (dyn == NULL)
	      continue;
	  for (i2 = i + 1; i2 < lns; i2++)
	    {
		/* testing interior/exterior relationships */
		dyn2 = dyn_array[i2];
		if (dyn2 == NULL)
		    continue;
		test_interior_ring (dyn, dyn2, &contains, &within, &crosses);
		if (contains)
		    ext_array[i2] = dyn;
		if (within)
		    ext_array[i] = dyn2;
		if (crosses)
		    ok = 0;
	    }
      }
    if (ok == 0)
      {
	  /* invalid: quitting */
	  for (i = 0; i < lns; i++)
	    {
		dyn = dyn_array[i];
		if (dyn == NULL)
		    continue;
		gaiaFreeDynamicLine (dyn);
	    }
	  free (dyn_array);
	  free (ext_array);
	  free (ln_array);
	  return NULL;
      }

    if (geom->DimensionModel == GAIA_XY_Z_M)
	result = gaiaAllocGeomCollXYZM ();
    else if (geom->DimensionModel == GAIA_XY_Z)
	result = gaiaAllocGeomCollXYZ ();
    else if (geom->DimensionModel == GAIA_XY_M)
	result = gaiaAllocGeomCollXYM ();
    else
	result = gaiaAllocGeomColl ();
    result->Srid = geom->Srid;
    if (force_multi)
	result->DeclaredType = GAIA_MULTIPOLYGON;

    for (i = 0; i < lns; i++)
      {
	  /* creating Polygons */
	  dyn = dyn_array[i];
	  if (dyn == NULL)
	      continue;
	  if (ext_array[i] != NULL)
	    {
		/* skipping any INTERIOR RING */
		continue;
	    }
	  pts = 0;
	  pt = dyn->First;
	  while (pt)
	    {
		/* counting how many points are there */
		pts++;
		pt = pt->Next;
	    }
	  num_interiors = 0;
	  for (i2 = 0; i2 < lns; i2++)
	    {
		if (ext_array[i2] == dyn)
		    num_interiors++;
	    }
	  pg = gaiaAddPolygonToGeomColl (result, pts, num_interiors);
	  rng = pg->Exterior;
	  iv = 0;
	  pt = dyn->First;
	  while (pt)
	    {
		/* EXTERIOR RING */
		if (result->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaSetPointXYZM (rng->Coords, iv, pt->X, pt->Y, pt->Z,
					pt->M);
		  }
		else if (result->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaSetPointXYZ (rng->Coords, iv, pt->X, pt->Y, pt->Z);
		  }
		else if (result->DimensionModel == GAIA_XY_M)
		  {
		      gaiaSetPointXYM (rng->Coords, iv, pt->X, pt->Y, pt->M);
		  }
		else
		  {
		      gaiaSetPoint (rng->Coords, iv, pt->X, pt->Y);
		  }
		iv++;
		pt = pt->Next;
	    }
	  ib = 0;
	  for (i2 = 0; i2 < lns; i2++)
	    {
		/* inserting any INTERIOR RING */
		if (ext_array[i2] == dyn)
		  {
		      dyn2 = dyn_array[i2];
		      ok = 1;
		      pts = 0;
		      pt = dyn2->First;
		      while (pt)
			{
			    /* counting how many points are there */
			    pts++;
			    pt = pt->Next;
			}
		      rng = gaiaAddInteriorRing (pg, ib, pts);
		      ib++;
		      iv = 0;
		      pt = dyn2->First;
		      while (pt)
			{
			    /* INTERIOR RING */
			    if (result->DimensionModel == GAIA_XY_Z_M)
			      {
				  gaiaSetPointXYZM (rng->Coords, iv, pt->X,
						    pt->Y, pt->Z, pt->M);
			      }
			    else if (result->DimensionModel == GAIA_XY_Z)
			      {
				  gaiaSetPointXYZ (rng->Coords, iv, pt->X,
						   pt->Y, pt->Z);
			      }
			    else if (result->DimensionModel == GAIA_XY_M)
			      {
				  gaiaSetPointXYM (rng->Coords, iv, pt->X,
						   pt->Y, pt->M);
			      }
			    else
			      {
				  gaiaSetPoint (rng->Coords, iv, pt->X, pt->Y);
			      }
			    iv++;
			    pt = pt->Next;
			}
		  }
	    }
      }

/* memory cleanup */
    for (i = 0; i < lns; i++)
      {
	  dyn = dyn_array[i];
	  if (dyn == NULL)
	      continue;
	  gaiaFreeDynamicLine (dyn);
      }
    free (dyn_array);
    free (ext_array);
    free (ln_array);
    if (result->FirstPolygon == NULL)
      {
	  gaiaFreeGeomColl (result);
	  return NULL;
      }
    return result;
}

#endif /* 0 */

#endif /* end including GEOS */
