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
	case GAIA_XY_M:
		return point(pt->X, pt->Y, 0.0, pt->M);
	case GAIA_XY_Z:
		return point(pt->X, pt->Y, pt->Z);
    case GAIA_XY_Z_M:
    	return point(pt->X, pt->Y, pt->Z, pt->M);
    default:
    	return point(pt->X, pt->Y);
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
			  res.push_back(point(v[0], v[1], 0.0, v[2]));
			  v += 3;
			  break;
		  case GAIA_XY_Z_M:
			  res.push_back(point(v[0], v[1], v[2], v[3]));
			  v += 4;
			  break;
		  default:
			  res.push_back(point(v[0], v[1]));
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
			  res.push_back(point(v[0], v[1], 0.0, v[2]));
			  v += 3;
			  break;
		  case GAIA_XY_Z_M:
			  res.push_back(point(v[0], v[1], v[2], v[3]));
			  v += 4;
			  break;
		  default:
			  res.push_back(point(v[0], v[1]));
			  v += 2;
			  break;
		  };
		}
	}

	return res;
}

static polygon toBoostPolygon(const gaiaPolygonPtr pg) {
	int i;
	polygon res;
	linearring ring;

	ring = toBoostRing(pg->Exterior);
	res.outer().swap(ring);

	for(i=0; i<pg->NumInteriors; i++) {
		res.inners().push_back(toBoostRing(pg->Interiors + i));
	}

	return res;
}


static multi_point toBoostMultiPoint(gaiaPointPtr pt) {
	multi_point res;
	while(pt) {
		res.push_back(toBoostPoint(pt));
		pt = pt->Next;
	}
	return res;
}

static multi_linestring toBoostMultiLinestring(gaiaLinestringPtr ln) {
	multi_linestring res;
	while(ln) {
		res.push_back(toBoostLinestring(ln));
		ln = ln->Next;
	}
	return res;
}

static multi_polygon toBoostMultiPolygon(gaiaPolygonPtr pg) {
	multi_polygon res;
	while(pg) {
		res.push_back(toBoostPolygon(pg));
		pg = pg->Next;
	}
	return res;
}

geometry_collection
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

    geometry_collection res;

	pt = gaia->FirstPoint;
	while(pt) {
		res.points.push_back(toBoostPoint(pt));
		pt = pt->Next;
	}

	ln = gaia->FirstLinestring;
	while(ln) {
		res.linestrings.push_back(toBoostLinestring(ln));
		ln = ln->Next;
	}

	pg = gaia->FirstPolygon;
	while(pg) {
		res.polygons.push_back(toBoostPolygon(pg));
		pg = pg->Next;
	}

	res.Srid = gaia->Srid;
	res.DimensionModel = gaia->DimensionModel;
	res.DeclaredType = gaia->DeclaredType;

	return res;
}


static gaiaPointPtr fromBoostPoint(const point& p) {
	gaiaPointPtr res = gaiaAllocPointXYZM(p.x, p.y, p.z, p.m);
	res->DimensionModel = p.DimensionModel;
	return res;
}

static gaiaLinestringPtr fromBoostLinestring(const linestring& line) {
	gaiaLinestringPtr res;
	int i;
	switch(line.DimensionModel) {
		case GAIA_XY_M:
			res = gaiaAllocLinestringXYM(line.size());
			break;
		case GAIA_XY_Z:
			res = gaiaAllocLinestringXYZ(line.size());
			break;
		case GAIA_XY_Z_M:
			res = gaiaAllocLinestringXYZM(line.size());
			break;
		default:
			res = gaiaAllocLinestring(line.size());
			break;
	}

	for(i=0; i<line.size(); i++) {
		const point& p=line[i];
		gaiaLineSetPoint (res, i, p.x, p.y, p.z, p.m);
		if (p.x < res->MinX) res->MinX = p.x;
		if (p.y < res->MinY) res->MinY = p.y;
		if (p.x > res->MaxX) res->MaxX = p.x;
		if (p.y > res->MaxY) res->MaxY = p.y;
	}

	return res;
}

static gaiaRingPtr fromBoostRing(const linearring& line, gaiaRingPtr dest=NULL) {
	gaiaRingPtr res;
	int i;

	if(dest) {
		res = dest;
	} else {
		switch(line.DimensionModel) {
			case GAIA_XY_M:
				res = gaiaAllocRingXYM(line.size());
				break;
			case GAIA_XY_Z:
				res = gaiaAllocRingXYZ(line.size());
				break;
			case GAIA_XY_Z_M:
				res = gaiaAllocRingXYZM(line.size());
				break;
			default:
				res = gaiaAllocRing(line.size());
				break;
		}
	}

	for(i=0; i<line.size(); i++) {
		const point& p=line[i];
		gaiaRingSetPoint (res, i, p.x, p.y, p.z, p.m);
		if (p.x < res->MinX) res->MinX = p.x;
		if (p.y < res->MinY) res->MinY = p.y;
		if (p.x > res->MaxX) res->MaxX = p.x;
		if (p.y > res->MaxY) res->MaxY = p.y;
	}

	return res;
}

static gaiaPolygonPtr fromBoostPolygon(const polygon& p) {
	gaiaPolygonPtr res;
	int i;

	switch(p.DimensionModel) {
		case GAIA_XY_M:
			res = gaiaAllocPolygonXYM(p.outer().size(), 0);
			break;
		case GAIA_XY_Z:
			res = gaiaAllocPolygonXYZ(p.outer().size(), 0);
			break;
		case GAIA_XY_Z_M:
			res = gaiaAllocPolygonXYZM(p.outer().size(), 0);
			break;
		default:
			res = gaiaAllocPolygon(p.outer().size(), 0);
			break;
	}

	fromBoostRing(p.outer(), res->Exterior);

	for(i=0; i<p.inners().size(); i++) {
		gaiaAddRingToPolyg(res, fromBoostRing(p.inners()[i]));
	}

	return res;
}

/* this function should be added to gaia core */
GAIAGEO_DECLARE gaiaPointPtr
gaiaInsertPointInGeomColl (gaiaGeomCollPtr c, gaiaPointPtr p) {
	if (c->FirstPoint == NULL)
		c->FirstPoint = p;
	if (c->LastPoint != NULL)
		c->LastPoint->Next = p;
	    c->LastPoint = p;
}

/** this function should be added to gaia core and replace the old InsertPolygonColl */
GAIAGEO_DECLARE gaiaPolygonPtr
gaiaInsertPolygonInGeomColl2 (gaiaGeomCollPtr p, gaiaPolygonPtr pg)
{
    if (p->FirstPolygon == NULL)
    	p->FirstPolygon = pg;
    if (p->LastPolygon != NULL)
    	p->LastPolygon->Next = pg;
    p->LastPolygon = pg;
    return pg;
}

gaiaGeomCollPtr
fromBoostGeometry (const geometry_collection& geom)
{
/* converting a GEOS Geometry into a GAIA Geometry */
   int i;
   gaiaGeomCollPtr gaia;

    gaia = gaiaAllocGeomColl();
    gaia->DimensionModel = geom.DimensionModel;
    gaia->DeclaredType = geom.DeclaredType;
    gaia->Srid = geom.Srid;

    for(i=0;i<geom.points.size();i++) {
    	gaiaInsertPointInGeomColl(gaia, fromBoostPoint(geom.points[i]));
    }

    for(i=0;i<geom.linestrings.size();i++) {
		gaiaInsertLinestringInGeomColl(gaia, fromBoostLinestring(geom.linestrings[i]));
	}

    for(i=0;i<geom.polygons.size();i++) {
    	gaiaInsertPolygonInGeomColl2(gaia, fromBoostPolygon(geom.polygons[i]));
    }


    return gaia;
}

#endif /* end including BOOST */
