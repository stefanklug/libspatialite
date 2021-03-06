/*

 spatialite.c -- SQLite3 spatial extension

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
#include <math.h>
#include <float.h>
#include <locale.h>
#include <errno.h>

#if defined(_WIN32) || defined(WIN32)
#include <io.h>
#define isatty	_isatty
#else
#include <unistd.h>
#endif

#ifdef SPL_AMALGAMATION		/* spatialite-amalgamation */
#include <spatialite/sqlite3ext.h>
#else
#include <sqlite3ext.h>
#endif

#include <spatialite/gaiaaux.h>
#include <spatialite/gaiageo.h>
#include <spatialite/gaiaexif.h>
#include <spatialite/spatialite.h>
#include <spatialite.h>

#ifndef OMIT_GEOS		/* including GEOS */
#include <geos_c.h>
#endif

#ifndef OMIT_PROJ		/* including PROJ.4 */
#include <proj_api.h>
#endif

#ifdef _WIN32
#define strcasecmp	_stricmp
#endif /* not WIN32 */

#define GAIA_UNUSED() if (argc || argv) argc = argc;

#ifndef OMIT_GEOCALLBACKS	/* supporting RTree geometry callbacks */
struct gaia_rtree_mbr
{
/* a struct used by R*Tree GeometryCallback functions [MBR] */
    double minx;
    double miny;
    double maxx;
    double maxy;
};
#endif /* end RTree geometry callbacks */

static SQLITE_EXTENSION_INIT1 struct spatial_index_str
{
/* a struct to implement a linked list of spatial-indexes */
    char ValidRtree;
    char ValidCache;
    char *TableName;
    char *ColumnName;
    struct spatial_index_str *Next;
};

struct stddev_str
{
/* a struct to implement StandardVariation and Variance aggregate functions */
    int cleaned;
    double mean;
    double quot;
    double count;
};

struct fdo_table
{
/* a struct to implement a linked-list for FDO-ORG table names */
    char *table;
    struct fdo_table *next;
};

static void
fnct_spatialite_version (sqlite3_context * context, int argc,
			 sqlite3_value ** argv)
{
/* SQL function:
/ spatialite_version()
/
/ return a text string representing the current SpatiaLite version
*/
    int len;
    const char *p_result = spatialite_version ();
    GAIA_UNUSED ();
    len = strlen (p_result);
    sqlite3_result_text (context, p_result, len, SQLITE_TRANSIENT);
}

static void
fnct_geos_version (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ geos_version()
/
/ return a text string representing the current GEOS version
/ or NULL if GEOS is currently unsupported
*/

#ifndef OMIT_GEOS		/* GEOS version */
    int len;
    const char *p_result = GEOSversion ();
    GAIA_UNUSED ();
    len = strlen (p_result);
    sqlite3_result_text (context, p_result, len, SQLITE_TRANSIENT);
#else
    sqlite3_result_null (context);
#endif
}


static void
fnct_proj4_version (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ proj4_version()
/
/ return a text string representing the current PROJ.4 version
/ or NULL if PROJ.4 is currently unsupported
*/

#ifndef OMIT_PROJ		/* PROJ.4 version */
    int len;
    const char *p_result = pj_get_release ();
    GAIA_UNUSED ();
    len = strlen (p_result);
    sqlite3_result_text (context, p_result, len, SQLITE_TRANSIENT);
#else
    sqlite3_result_null (context);
#endif
}

static void
clean_sql_string (char *buf)
{
/* well-formatting a string to be used as an SQL string-value */
    char tmp[1024];
    char *in = tmp;
    char *out = buf;
    strcpy (tmp, buf);
    while (*in != '\0')
      {
	  if (*in == '\'')
	      *out++ = '\'';
	  *out++ = *in++;
      }
    *out = '\0';
}

static void
double_quoted_sql (char *buf)
{
/* well-formatting a string to be used as an SQL name */
    char tmp[1024];
    char *in = tmp;
    char *out = buf;
    strcpy (tmp, buf);
    *out++ = '"';
    while (*in != '\0')
      {
	  if (*in == '"')
	      *out++ = '"';
	  *out++ = *in++;
      }
    *out++ = '"';
    *out = '\0';
}

static void
fnct_GeometryConstraints (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ GeometryConstraints(BLOBencoded geometry, geometry-type, srid)
/ GeometryConstraints(BLOBencoded geometry, geometry-type, srid, dimensions)
/
/ checks geometry constraints, returning:
/
/ -1 - if some error occurred
/ 1 - if geometry constraints validation passes
/ 0 - if geometry constraints validation fails
/
*/
    int little_endian;
    int endian_arch = gaiaEndianArch ();
    unsigned char *p_blob = NULL;
    int n_bytes = 0;
    int srid;
    int geom_srid = -1;
    const unsigned char *type;
    int xtype;
    int geom_type = -1;
    int geom_normalized_type;
    const unsigned char *dimensions;
    int dims = GAIA_XY;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_BLOB
	|| sqlite3_value_type (argv[0]) == SQLITE_NULL)
	;
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_TEXT)
	type = sqlite3_value_text (argv[1]);
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	srid = sqlite3_value_int (argv[2]);
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (argc == 4)
      {
	  /* explicit dimensions - supporting XYZM */
	  dimensions = sqlite3_value_text (argv[3]);
	  if (strcasecmp ((char *) dimensions, "XYZ") == 0)
	      dims = GAIA_XY_Z;
	  else if (strcasecmp ((char *) dimensions, "XYM") == 0)
	      dims = GAIA_XY_M;
	  else if (strcasecmp ((char *) dimensions, "XYZM") == 0)
	      dims = GAIA_XY_Z_M;
	  else
	      dims = GAIA_XY;
      }
    if (sqlite3_value_type (argv[0]) == SQLITE_BLOB)
      {
	  p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
	  n_bytes = sqlite3_value_bytes (argv[0]);
      }
    if (p_blob)
      {
	  /* quick Geometry validation */
	  if (n_bytes < 45)
	      goto illegal_geometry;	/* cannot be an internal BLOB WKB geometry */
	  if (*(p_blob + 0) != GAIA_MARK_START)
	      goto illegal_geometry;	/* failed to recognize START signature */
	  if (*(p_blob + (n_bytes - 1)) != GAIA_MARK_END)
	      goto illegal_geometry;	/* failed to recognize END signature */
	  if (*(p_blob + 38) != GAIA_MARK_MBR)
	      goto illegal_geometry;	/* failed to recognize MBR signature */
	  if (*(p_blob + 1) == GAIA_LITTLE_ENDIAN)
	      little_endian = 1;
	  else if (*(p_blob + 1) == GAIA_BIG_ENDIAN)
	      little_endian = 0;
	  else
	      goto illegal_geometry;	/* unknown encoding; neither little-endian nor big-endian */
	  geom_type = gaiaImport32 (p_blob + 39, little_endian, endian_arch);
	  geom_srid = gaiaImport32 (p_blob + 2, little_endian, endian_arch);
	  goto valid_geometry;
	illegal_geometry:
	  sqlite3_result_int (context, -1);
	  return;
      }
  valid_geometry:
    xtype = GAIA_UNKNOWN;
    if (strcasecmp ((char *) type, "POINT") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_POINTZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_POINTM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_POINTZM;
		break;
	    default:
		xtype = GAIA_POINT;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "LINESTRING") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_LINESTRINGZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_LINESTRINGM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_LINESTRINGZM;
		break;
	    default:
		xtype = GAIA_LINESTRING;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "POLYGON") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_POLYGONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_POLYGONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_POLYGONZM;
		break;
	    default:
		xtype = GAIA_POLYGON;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "MULTIPOINT") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTIPOINTZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTIPOINTM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTIPOINTZM;
		break;
	    default:
		xtype = GAIA_MULTIPOINT;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "MULTILINESTRING") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTILINESTRINGZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTILINESTRINGM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTILINESTRINGZM;
		break;
	    default:
		xtype = GAIA_MULTILINESTRING;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "MULTIPOLYGON") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTIPOLYGONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTIPOLYGONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTIPOLYGONZM;
		break;
	    default:
		xtype = GAIA_MULTIPOLYGON;
		break;
	    };
      }
    if (strcasecmp ((char *) type, "GEOMETRYCOLLECTION") == 0)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_GEOMETRYCOLLECTIONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_GEOMETRYCOLLECTIONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_GEOMETRYCOLLECTIONZM;
		break;
	    default:
		xtype = GAIA_GEOMETRYCOLLECTION;
		break;
	    };
      }
    switch (geom_type)
      {
	  /* adjusting COMPRESSED Geometries */
      case GAIA_COMPRESSED_LINESTRING:
	  geom_normalized_type = GAIA_LINESTRING;
	  break;
      case GAIA_COMPRESSED_LINESTRINGZ:
	  geom_normalized_type = GAIA_LINESTRINGZ;
	  break;
      case GAIA_COMPRESSED_LINESTRINGM:
	  geom_normalized_type = GAIA_LINESTRINGM;
	  break;
      case GAIA_COMPRESSED_LINESTRINGZM:
	  geom_normalized_type = GAIA_LINESTRINGZM;
	  break;
      case GAIA_COMPRESSED_POLYGON:
	  geom_normalized_type = GAIA_POLYGON;
	  break;
      case GAIA_COMPRESSED_POLYGONZ:
	  geom_normalized_type = GAIA_POLYGONZ;
	  break;
      case GAIA_COMPRESSED_POLYGONM:
	  geom_normalized_type = GAIA_POLYGONM;
	  break;
      case GAIA_COMPRESSED_POLYGONZM:
	  geom_normalized_type = GAIA_POLYGONZM;
	  break;
      default:
	  geom_normalized_type = geom_type;
	  break;
      };
    if (strcasecmp ((char *) type, "GEOMETRY") == 0)
	xtype = -1;
    if (xtype == GAIA_UNKNOWN)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = 1;
	  if (p_blob)
	    {
		/* skipping NULL Geometry; this is assumed to be always good */
		if (geom_srid != srid)
		    ret = 0;
		if (xtype == -1)
		    ;
		else if (xtype != geom_normalized_type)
		    ret = 0;
	    }
	  sqlite3_result_int (context, ret);
      }
}

static void
fnct_RTreeAlign (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ RTreeAlign(RTree-table-name, PKID-value, BLOBencoded geometry)
/
/ attempts to update the associated R*Tree, returning:
/
/ -1 - if some invalid arg was passed
/ 1 - succesfull update
/ 0 - update failure
/
*/
    unsigned char *p_blob = NULL;
    int n_bytes = 0;
    sqlite3_int64 pkid;
    const unsigned char *rtree_table;
    gaiaGeomCollPtr geom = NULL;
    int ret;
    char table_name[1024];
    char sql[4192];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_TEXT)
	rtree_table = sqlite3_value_text (argv[0]);
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	pkid = sqlite3_value_int64 (argv[1]);
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_BLOB
	|| sqlite3_value_type (argv[2]) == SQLITE_NULL)
	;
    else
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_BLOB)
      {
	  p_blob = (unsigned char *) sqlite3_value_blob (argv[2]);
	  n_bytes = sqlite3_value_bytes (argv[2]);
	  geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
      }

    if (geom == NULL)
      {
	  /* NULL geometry: nothing to do */
	  sqlite3_result_int (context, 1);
      }
    else
      {
	  /* INSERTing into the R*Tree */
	  strcpy (table_name, rtree_table);
	  if (*(table_name + 0) == '"'
	      && *(table_name + strlen (table_name) - 1) == '"')
	      ;			/* earlier versions may pass an already quoted name */
	  else
	      double_quoted_sql (table_name);
#if defined(_WIN32) || defined(__MINGW32__)
/* CAVEAT: M$ runtime doesn't supports %lld for 64 bits */
	  sprintf (sql, "INSERT INTO %s (pkid, xmin, ymin, xmax, ymax) "
		   "VALUES (%I64d, %1.12f, %1.12f, %1.12f, %1.12f)",
		   table_name, pkid, geom->MinX, geom->MinY, geom->MaxX,
		   geom->MaxY);
#else
	  sprintf (sql, "INSERT INTO %s (pkid, xmin, ymin, xmax, ymax) "
		   "VALUES (%lld, %1.12f, %1.12f, %1.12f, %1.12f)",
		   table_name, pkid, geom->MinX, geom->MinY, geom->MaxX,
		   geom->MaxY);
#endif
	  gaiaFreeGeomColl (geom);
	  ret = sqlite3_exec (sqlite, sql, NULL, NULL, NULL);
	  if (ret != SQLITE_OK)
	      sqlite3_result_int (context, 0);
	  else
	      sqlite3_result_int (context, 1);
      }
}

static int
checkSpatialMetaData (sqlite3 * sqlite)
{
/* internal utility function:
/
/ for FDO-OGR interoperability:
/ tests the SpatialMetadata type, returning:
/
/ 0 - if no valid SpatialMetaData where found
/ 1 - if SpatiaLite-like SpatialMetadata where found
/ 2- if FDO-OGR-like SpatialMetadata where found
/
*/
    int spatialite_rs = 0;
    int fdo_rs = 0;
    int spatialite_gc = 0;
    int fdo_gc = 0;
    int rs_srid = 0;
    int auth_name = 0;
    int auth_srid = 0;
    int srtext = 0;
    int ref_sys_name = 0;
    int proj4text = 0;
    int f_table_name = 0;
    int f_geometry_column = 0;
    int geometry_type = 0;
    int coord_dimension = 0;
    int gc_srid = 0;
    int geometry_format = 0;
    int type = 0;
    int spatial_index_enabled = 0;
    char sql[1024];
    int ret;
    const char *name;
    int i;
    char **results;
    int rows;
    int columns;
/* checking the GEOMETRY_COLUMNS table */
    strcpy (sql, "PRAGMA table_info(geometry_columns)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	goto unknown;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "f_table_name") == 0)
		    f_table_name = 1;
		if (strcasecmp (name, "f_geometry_column") == 0)
		    f_geometry_column = 1;
		if (strcasecmp (name, "geometry_type") == 0)
		    geometry_type = 1;
		if (strcasecmp (name, "coord_dimension") == 0)
		    coord_dimension = 1;
		if (strcasecmp (name, "srid") == 0)
		    gc_srid = 1;
		if (strcasecmp (name, "geometry_format") == 0)
		    geometry_format = 1;
		if (strcasecmp (name, "type") == 0)
		    type = 1;
		if (strcasecmp (name, "spatial_index_enabled") == 0)
		    spatial_index_enabled = 1;
	    }
      }
    sqlite3_free_table (results);
    if (f_table_name
	&&
	f_geometry_column
	&& type && coord_dimension && gc_srid && spatial_index_enabled)
	spatialite_gc = 1;
    if (f_table_name
	&&
	f_geometry_column
	&& geometry_type && coord_dimension && gc_srid && geometry_format)
	fdo_gc = 1;
/* checking the SPATIAL_REF_SYS table */
    strcpy (sql, "PRAGMA table_info(spatial_ref_sys)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	goto unknown;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "srid") == 0)
		    rs_srid = 1;
		if (strcasecmp (name, "auth_name") == 0)
		    auth_name = 1;
		if (strcasecmp (name, "auth_srid") == 0)
		    auth_srid = 1;
		if (strcasecmp (name, "srtext") == 0)
		    srtext = 1;
		if (strcasecmp (name, "ref_sys_name") == 0)
		    ref_sys_name = 1;
		if (strcasecmp (name, "proj4text") == 0)
		    proj4text = 1;
	    }
      }
    sqlite3_free_table (results);
    if (rs_srid && auth_name && auth_srid && ref_sys_name && proj4text)
	spatialite_rs = 1;
    if (rs_srid && auth_name && auth_srid && srtext)
	fdo_rs = 1;
/* verifying the MetaData format */
    if (spatialite_gc && spatialite_rs)
	return 1;
    if (fdo_gc && fdo_rs)
	return 2;
  unknown:
    return 0;
}

static void
add_fdo_table (struct fdo_table **first, struct fdo_table **last,
	       const char *table, int len)
{
/* adds an FDO-OGR styled Geometry Table to corresponding linked list */
    struct fdo_table *p = malloc (sizeof (struct fdo_table));
    p->table = malloc (len + 1);
    strcpy (p->table, table);
    p->next = NULL;
    if (!(*first))
	(*first) = p;
    if ((*last))
	(*last)->next = p;
    (*last) = p;
}

static void
free_fdo_tables (struct fdo_table *first)
{
/* memory cleanup; destroying the FDO-OGR tables linked list */
    struct fdo_table *p;
    struct fdo_table *pn;
    p = first;
    while (p)
      {
	  pn = p->next;
	  if (p->table)
	      free (p->table);
	  free (p);
	  p = pn;
      }
}

static void
fnct_AutoFDOStart (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AutoFDOStart(void)
/
/ for FDO-OGR interoperability:
/ tests the SpatialMetadata type, then automatically
/ creating a VirtualFDO table for each FDO-OGR main table 
/ declared within FDO-styled SpatialMetadata
/
*/
    int ret;
    const char *name;
    int i;
    char **results;
    int rows;
    int columns;
    char sql[1024];
    int count = 0;
    struct fdo_table *first = NULL;
    struct fdo_table *last = NULL;
    struct fdo_table *p;
    int len;
    char xname[1024];
    char xtable[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (checkSpatialMetaData (sqlite) == 2)
      {
	  /* ok, creating VirtualFDO tables */
	  strcpy (sql, "SELECT DISTINCT f_table_name FROM geometry_columns");
	  ret =
	      sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (rows < 1)
	      ;
	  else
	    {
		for (i = 1; i <= rows; i++)
		  {
		      name = results[(i * columns) + 0];
		      if (name)
			{
			    len = strlen (name);
			    add_fdo_table (&first, &last, name, len);
			}
		  }
	    }
	  sqlite3_free_table (results);
	  p = first;
	  while (p)
	    {
		/* destroying the VirtualFDO table [if existing] */
		sprintf (xname, "fdo_%s", p->table);
		double_quoted_sql (xname);
		sprintf (sql, "DROP TABLE IF EXISTS %s", xname);
		ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
		if (ret != SQLITE_OK)
		    goto error;
		/* creating the VirtualFDO table */
		strcpy (xtable, p->table);
		double_quoted_sql (xtable);
		sprintf (sql, "CREATE VIRTUAL TABLE %s USING VirtualFDO(%s)",
			 xname, xtable);
		ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
		if (ret != SQLITE_OK)
		    goto error;
		count++;
		p = p->next;
	    }
	error:
	  free_fdo_tables (first);
	  sqlite3_result_int (context, count);
	  return;
      }
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_AutoFDOStop (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AutoFDOStop(void)
/
/ for FDO-OGR interoperability:
/ tests the SpatialMetadata type, then automatically
/ removes any VirtualFDO table 
/
*/
    int ret;
    const char *name;
    int i;
    char **results;
    int rows;
    int columns;
    char sql[1024];
    int count = 0;
    struct fdo_table *first = NULL;
    struct fdo_table *last = NULL;
    struct fdo_table *p;
    int len;
    char xname[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (checkSpatialMetaData (sqlite) == 2)
      {
	  /* ok, creating VirtualFDO tables */
	  strcpy (sql, "SELECT DISTINCT f_table_name FROM geometry_columns");
	  ret =
	      sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (rows < 1)
	      ;
	  else
	    {
		for (i = 1; i <= rows; i++)
		  {
		      name = results[(i * columns) + 0];
		      if (name)
			{
			    len = strlen (name);
			    add_fdo_table (&first, &last, name, len);
			}
		  }
	    }
	  sqlite3_free_table (results);
	  p = first;
	  while (p)
	    {
		/* destroying the VirtualFDO table [if existing] */
		sprintf (xname, "fdo_%s", p->table);
		double_quoted_sql (xname);
		sprintf (sql, "DROP TABLE IF EXISTS %s", xname);
		ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
		if (ret != SQLITE_OK)
		    goto error;
		count++;
		p = p->next;
	    }
	error:
	  free_fdo_tables (first);
	  sqlite3_result_int (context, count);
	  return;
      }
    sqlite3_result_int (context, 0);
    return;
}

static int
testSpatiaLiteHistory (sqlite3 * sqlite)
{
/* internal utility function:
/
/ checks if the SPATIALITE_HISTORY table already exists
/
*/
    int event_id = 0;
    int table_name = 0;
    int geometry_column = 0;
    int event = 0;
    int timestamp = 0;
    int ver_sqlite = 0;
    int ver_splite = 0;
    char sql[1024];
    int ret;
    const char *name;
    int i;
    char **results;
    int rows;
    int columns;
/* checking the SPATIALITE_HISTORY table */
    strcpy (sql, "PRAGMA table_info(spatialite_history)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return 0;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "event_id") == 0)
		    event_id = 1;
		if (strcasecmp (name, "table_name") == 0)
		    table_name = 1;
		if (strcasecmp (name, "geometry_column") == 0)
		    geometry_column = 1;
		if (strcasecmp (name, "event") == 0)
		    event = 1;
		if (strcasecmp (name, "timestamp") == 0)
		    timestamp = 1;
		if (strcasecmp (name, "ver_sqlite") == 0)
		    ver_sqlite = 1;
		if (strcasecmp (name, "ver_splite") == 0)
		    ver_splite = 1;
	    }
      }
    sqlite3_free_table (results);
    if (event_id && table_name && geometry_column && event && timestamp
	&& ver_sqlite && ver_splite)
	return 1;
    return 0;
}

static int
checkSpatiaLiteHistory (sqlite3 * sqlite)
{
/* internal utility function:
/
/ checks if the SPATIALITE_HISTORY table already exists
/ if not, such table will then be created
/
*/
    char sql[1024];
    char *errMsg = NULL;
    int ret;

    if (testSpatiaLiteHistory (sqlite))
	return 1;

/* creating the SPATIALITE_HISTORY table */
    strcpy (sql, "CREATE TABLE IF NOT EXISTS ");
    strcat (sql, "spatialite_history (\n");
    strcat (sql, "event_id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,\n");
    strcat (sql, "table_name TEXT NOT NULL,\n");
    strcat (sql, "geometry_column TEXT,\n");
    strcat (sql, "event TEXT NOT NULL,\n");
    strcat (sql, "timestamp TEXT NOT NULL,\n");
    strcat (sql, "ver_sqlite TEXT NOT NULL,\n");
    strcat (sql, "ver_splite TEXT NOT NULL)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;

    if (testSpatiaLiteHistory (sqlite))
	return 1;
    return 0;
}

static void
updateSpatiaLiteHistory (sqlite3 * sqlite, const char *table,
			 const char *geom, const char *operation)
{
/* inserting a row in SPATIALITE_HISTORY */
    char sql[2048];
    sqlite3_stmt *stmt = NULL;
    int ret;

    if (checkSpatiaLiteHistory (sqlite) == 0)
	return;

    strcpy (sql, "INSERT INTO spatialite_history ");
    strcat (sql, "(event_id, table_name, geometry_column, event, timestamp, ");
    strcat (sql, "ver_sqlite, ver_splite) ");
    strcat (sql,
	    "VALUES (NULL, ?, ?, ?, DateTime('now'), sqlite_version(), spatialite_version())");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "SQL error: %s\n%s\n", sql, sqlite3_errmsg (sqlite));
	  goto stop;
      }
    sqlite3_reset (stmt);
    sqlite3_clear_bindings (stmt);
    sqlite3_bind_text (stmt, 1, table, strlen (table), SQLITE_STATIC);
    if (!geom)
	sqlite3_bind_null (stmt, 2);
    else
	sqlite3_bind_text (stmt, 2, geom, strlen (geom), SQLITE_STATIC);
    sqlite3_bind_text (stmt, 3, operation, strlen (operation), SQLITE_STATIC);
    ret = sqlite3_step (stmt);
    if (ret == SQLITE_DONE || ret == SQLITE_ROW)
	goto stop;
    fprintf (stderr, "SQL error: %s\n", sqlite3_errmsg (sqlite));

  stop:
    if (stmt)
	sqlite3_finalize (stmt);
}

static int
createAdvancedMetaData (sqlite3 * sqlite)
{
/* creating the advanced MetaData tables */
    char sql[1024];
    char *errMsg = NULL;
    int ret;
/* creating the VIEWS_GEOMETRY_COLUMNS table */
    strcpy (sql, "CREATE TABLE IF NOT EXISTS ");
    strcat (sql, "views_geometry_columns (\n");
    strcat (sql, "view_name TEXT NOT NULL,\n");
    strcat (sql, "view_geometry TEXT NOT NULL,\n");
    strcat (sql, "view_rowid TEXT NOT NULL,\n");
    strcat (sql, "f_table_name VARCHAR(256) NOT NULL,\n");
    strcat (sql, "f_geometry_column VARCHAR(256) NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_geom_cols_views PRIMARY KEY ");
    strcat (sql, "(view_name, view_geometry),\n");
    strcat (sql, "CONSTRAINT fk_views_geom_cols FOREIGN KEY ");
    strcat (sql, "(f_table_name, f_geometry_column) ");
    strcat (sql, "REFERENCES geometry_columns ");
    strcat (sql, "(f_table_name, f_geometry_column) ");
    strcat (sql, "ON DELETE CASCADE)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
/* creating an INDEX supporting the GEOMETRY_COLUMNS FK */
    strcpy (sql, "CREATE INDEX IF NOT EXISTS ");
    strcat (sql, "idx_viewsjoin ON views_geometry_columns\n");
    strcat (sql, "(f_table_name, f_geometry_column)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
/* creating the VIRTS_GEOMETRY_COLUMNS table */
    strcpy (sql, "CREATE TABLE IF NOT EXISTS ");
    strcat (sql, "virts_geometry_columns (\n");
    strcat (sql, "virt_name TEXT NOT NULL,\n");
    strcat (sql, "virt_geometry TEXT NOT NULL,\n");
    strcat (sql, "type VARCHAR(30) NOT NULL,\n");
    strcat (sql, "srid INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_geom_cols_virts PRIMARY KEY ");
    strcat (sql, "(virt_name, virt_geometry),\n");
    strcat (sql, "CONSTRAINT fk_vgc_srid FOREIGN KEY ");
    strcat (sql, "(srid) REFERENCES spatial_ref_sys (srid))");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
/* creating an INDEX supporting the SPATIAL_REF_SYS FK */
    strcpy (sql, "CREATE INDEX IF NOT EXISTS ");
    strcat (sql, "idx_virtssrid ON virts_geometry_columns\n");
    strcat (sql, "(srid)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
/* creating the GEOMETRY_COLUMNS_AUTH table */
    strcpy (sql, "CREATE TABLE IF NOT EXISTS ");
    strcat (sql, "geometry_columns_auth (\n");
    strcat (sql, "f_table_name VARCHAR(256) NOT NULL,\n");
    strcat (sql, "f_geometry_column VARCHAR(256) NOT NULL,\n");
    strcat (sql, "read_only INTEGER NOT NULL,\n");
    strcat (sql, "hidden INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_gc_auth PRIMARY KEY ");
    strcat (sql, "(f_table_name, f_geometry_column),\n");
    strcat (sql, "CONSTRAINT fk_gc_auth FOREIGN KEY ");
    strcat (sql, "(f_table_name, f_geometry_column) ");
    strcat (sql, "REFERENCES geometry_columns ");
    strcat (sql, "(f_table_name, f_geometry_column) ");
    strcat (sql, "ON DELETE CASCADE)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
    return 1;
}

static void
fnct_CheckSpatialMetaData (sqlite3_context * context, int argc,
			   sqlite3_value ** argv)
{
/* SQL function:
/ CheckSpatialMetaData(void)
/
/ for FDO-OGR interoperability:
/ tests the SpatialMetadata type, returning:
/
/ 0 - if no valid SpatialMetaData where found
/ 1 - if SpatiaLite-like SpatialMetadata where found
/ 2- if FDO-OGR-like SpatialMetadata where found
/
*/
    sqlite3 *sqlite;
    int ret;
    GAIA_UNUSED ();
    sqlite = sqlite3_context_db_handle (context);
    ret = checkSpatialMetaData (sqlite);
    if (ret == 1)
      {
	  /* trying to create the advanced metadata tables */
	  createAdvancedMetaData (sqlite);
      }
    sqlite3_result_int (context, ret);
    return;
}

static void
fnct_InitSpatialMetaData (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ InitSpatialMetaData(void)
/
/ creates the SPATIAL_REF_SYS and GEOMETRY_COLUMNS tables
/ returns 1 on success
/ 0 on failure
*/
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
/* creating the SPATIAL_REF_SYS table */
    strcpy (sql, "CREATE TABLE spatial_ref_sys (\n");
    strcat (sql, "srid INTEGER NOT NULL PRIMARY KEY,\n");
    strcat (sql, "auth_name TEXT NOT NULL,\n");
    strcat (sql, "auth_srid INTEGER NOT NULL,\n");
    strcat (sql, "ref_sys_name TEXT,\n");
    strcat (sql, "proj4text TEXT NOT NULL,\n");
    strcat (sql, "srs_wkt TEXT) ");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    strcpy (sql, "CREATE UNIQUE INDEX idx_spatial_ref_sys \n");
    strcat (sql, "ON spatial_ref_sys (auth_srid, auth_name)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    updateSpatiaLiteHistory (sqlite, "spatial_ref_sys", NULL,
			     "table successfully created");
/* creating the GEOMETRY_COLUMN table */
    strcpy (sql, "CREATE TABLE geometry_columns (\n");
    strcat (sql, "f_table_name TEXT NOT NULL,\n");
    strcat (sql, "f_geometry_column TEXT NOT NULL,\n");
    strcat (sql, "type TEXT NOT NULL,\n");
    strcat (sql, "coord_dimension TEXT NOT NULL,\n");
    strcat (sql, "srid INTEGER NOT NULL,\n");
    strcat (sql, "spatial_index_enabled INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_geom_cols PRIMARY KEY ");
    strcat (sql, "(f_table_name, f_geometry_column),\n");
    strcat (sql, "CONSTRAINT fk_gc_srs FOREIGN KEY ");
    strcat (sql, "(srid) REFERENCES spatial_ref_sys (srid))");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    updateSpatiaLiteHistory (sqlite, "geometry_columns", NULL,
			     "table successfully created");
/* creating an INDEX corresponding to the SRID FK */
    strcpy (sql, "CREATE INDEX idx_srid_geocols ON geometry_columns\n");
    strcat (sql, "(srid) ");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/* creating the GEOM_COLS_REF_SYS view */
    strcpy (sql, "CREATE VIEW geom_cols_ref_sys AS\n");
    strcat (sql, "SELECT f_table_name, f_geometry_column, type,\n");
    strcat (sql, "coord_dimension, spatial_ref_sys.srid AS srid,\n");
    strcat (sql, "auth_name, auth_srid, ref_sys_name, proj4text\n");
    strcat (sql, "FROM geometry_columns, spatial_ref_sys\n");
    strcat (sql, "WHERE geometry_columns.srid = spatial_ref_sys.srid");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    if (!createAdvancedMetaData (sqlite))
	goto error;
/* creating the SpatialIndex VIRTUAL TABLE */
    strcpy (sql, "CREATE VIRTUAL TABLE SpatialIndex ");
    strcat (sql, "USING VirtualSpatialIndex()");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    if (spatial_ref_sys_init (sqlite, 0))
	updateSpatiaLiteHistory (sqlite, "spatial_ref_sys", NULL,
				 "table successfully populated");
    sqlite3_result_int (context, 1);
    return;
  error:
    fprintf (stderr, " InitSpatiaMetaData ()error:\"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static int
recoverGeomColumn (sqlite3 * sqlite, const unsigned char *table,
		   const unsigned char *column, int xtype, int dims, int srid)
{
/* checks if TABLE.COLUMN exists and has the required features */
    int ok = 1;
    char sql[1024];
    int type;
    sqlite3_stmt *stmt;
    gaiaGeomCollPtr geom;
    const void *blob_value;
    int len;
    int ret;
    int i_col;
    char xcolumn[1024];
    char xtable[1024];
    strcpy (xcolumn, (char *) column);
    double_quoted_sql (xcolumn);
    strcpy (xtable, (char *) table);
    double_quoted_sql (xtable);
    sprintf (sql, "SELECT %s FROM %s", xcolumn, xtable);
/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "recoverGeomColumn: error %d \"%s\"\n",
		   sqlite3_errcode (sqlite), sqlite3_errmsg (sqlite));
	  return 0;
      }
    while (1)
      {
	  /* scrolling the result set rows */
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;		/* end of result set */
	  if (ret == SQLITE_ROW)
	    {
		/* checking Geometry features */
		geom = NULL;
		for (i_col = 0; i_col < sqlite3_column_count (stmt); i_col++)
		  {
		      if (sqlite3_column_type (stmt, i_col) != SQLITE_BLOB)
			  ok = 0;
		      else
			{
			    blob_value = sqlite3_column_blob (stmt, i_col);
			    len = sqlite3_column_bytes (stmt, i_col);
			    geom = gaiaFromSpatiaLiteBlobWkb (blob_value, len);
			    if (!geom)
				ok = 0;
			    else
			      {
				  if (geom->DimensionModel != dims)
				      ok = 0;
				  if (geom->Srid != srid)
				      ok = 0;
				  type = gaiaGeometryType (geom);
				  if (xtype == -1)
				      ;	/* GEOMETRY */
				  else
				    {
					if (xtype == type)
					    ;
					else
					    ok = 0;
				    }
				  gaiaFreeGeomColl (geom);
			      }
			}
		  }
	    }
	  if (!ok)
	      break;
      }
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "recoverGeomColumn: error %d \"%s\"\n",
		   sqlite3_errcode (sqlite), sqlite3_errmsg (sqlite));
	  return 0;
      }
    return ok;
}

static void
buildSpatialIndex (sqlite3 * sqlite, const unsigned char *table,
		   const char *col_name)
{
/* loading a SpatialIndex [RTree] */
    char sql[2048];
    char sql2[1024];
    char *errMsg = NULL;
    int ret;
    char xname[1024];
    char xtable[1024];
    sprintf (xname, "idx_%s_%s", table, col_name);
    double_quoted_sql (xname);
    sprintf (sql, "INSERT INTO %s (pkid, xmin, xmax, ymin, ymax) ", xname);
    strcpy (xname, col_name);
    double_quoted_sql (xname);
    strcpy (xtable, (char *) table);
    double_quoted_sql (xtable);
    sprintf (sql2,
	     "SELECT ROWID, MbrMinX(%s), MbrMaxX(%s), MbrMinY(%s), MbrMaxY(%s) FROM %s",
	     xname, xname, xname, xname, xtable);
    strcat (sql, sql2);
    sprintf (sql2, " WHERE MbrMinX(%s) IS NOT NULL", xname);
    strcat (sql, sql2);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "buildSpatialIndex error: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
      }
}

static void
updateGeometryTriggers (sqlite3 * sqlite, const unsigned char *table,
			const unsigned char *column)
{
/* updates triggers for some Spatial Column */
    char sql[256];
    char trigger[4096];
    char **results;
    int ret;
    int rows;
    int columns;
    int i;
    char tblname[256];
    char colname[256];
    char col_type[32];
    char col_srid[32];
    char col_index[32];
    char col_dims[64];
    int index;
    int cached;
    int dims;
    char *txt_dims;
    int len;
    char *errMsg = NULL;
    char dummy[512];
    char sqltable[1024];
    char sqlcolumn[1024];
    char xname[1024];
    char xcolname[1024];
    char xtable[1024];
    char xindex[1024];
    struct spatial_index_str *first_idx = NULL;
    struct spatial_index_str *last_idx = NULL;
    struct spatial_index_str *curr_idx;
    struct spatial_index_str *next_idx;
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "SELECT f_table_name, f_geometry_column, type, srid, spatial_index_enabled, coord_dimension "
	     "FROM geometry_columns WHERE f_table_name LIKE '%s' AND f_geometry_column LIKE '%s'",
	     sqltable, sqlcolumn);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "updateTableTriggers: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    for (i = 1; i <= rows; i++)
      {
	  /* preparing the triggers */
	  strcpy (tblname, results[(i * columns)]);
	  strcpy (colname, results[(i * columns) + 1]);
	  strcpy (col_type, results[(i * columns) + 2]);
	  /* 
	     / Even Rouault - 3 Mar 2010 
	     / the OGR driver wrongly inserts a NULL SRID
	     / into GEOMETRY_COLUMNS, so we must check such
	     / an odd condition to avoid a crash
	   */
	  if (results[(i * columns) + 3] == NULL)
	      strcpy (col_srid, "-1");
	  else
	      strcpy (col_srid, results[(i * columns) + 3]);
	  strcpy (col_index, results[(i * columns) + 4]);
	  strcpy (col_dims, results[(i * columns) + 5]);
	  if (atoi (col_index) == 1)
	      index = 1;
	  else
	      index = 0;
	  if (atoi (col_index) == 2)
	      cached = 1;
	  else
	      cached = 0;
	  dims = GAIA_XY;
	  if (strcasecmp (col_dims, "XYZ") == 0)
	      dims = GAIA_XY_Z;
	  if (strcasecmp (col_dims, "XYM") == 0)
	      dims = GAIA_XY_M;
	  if (strcasecmp (col_dims, "XYZM") == 0)
	      dims = GAIA_XY_Z_M;
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		txt_dims = "XYZ";
		break;
	    case GAIA_XY_M:
		txt_dims = "XYM";
		break;
	    case GAIA_XY_Z_M:
		txt_dims = "XYZM";
		break;
	    default:
		txt_dims = "XY";
		break;
	    };

	  /* trying to delete old versions [v2.0, v2.2] triggers[if any] */
	  strcpy (sqltable, (char *) tblname);
	  clean_sql_string (sqltable);
	  strcpy (sqlcolumn, (char *) colname);
	  clean_sql_string (sqlcolumn);
	  sprintf (xname, "gti_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  sprintf (xname, "gtu_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  sprintf (xname, "gsi_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  sprintf (xname, "gsu_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  /* end deletion old versions [v2.0, v2.2] triggers[if any] */

	  /* deleting the old INSERT trigger TYPE [if any] */
	  sprintf (xname, "ggi_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  /* inserting the new INSERT trigger TYPE */
	  strcpy (xtable, tblname);
	  double_quoted_sql (xtable);
	  strcpy (xcolname, colname);
	  double_quoted_sql (xcolname);
	  sprintf (trigger, "CREATE TRIGGER %s BEFORE INSERT ON %s\n", xname,
		   xtable);
	  strcat (trigger, "FOR EACH ROW BEGIN\n");
	  sprintf (dummy,
		   "SELECT RAISE(ROLLBACK, '%s.%s violates Geometry constraint [geom-type or SRID not allowed]')\n",
		   sqltable, sqlcolumn);
	  strcat (trigger, dummy);
	  strcat (trigger, "WHERE (SELECT type FROM geometry_columns\n");
	  sprintf (dummy,
		   "WHERE f_table_name = '%s' AND f_geometry_column = '%s'\n",
		   sqltable, sqlcolumn);
	  strcat (trigger, dummy);
	  sprintf (dummy,
		   "AND GeometryConstraints(NEW.%s, type, srid, '%s') = 1) IS NULL;\n",
		   xcolname, txt_dims);
	  strcat (trigger, dummy);
	  strcat (trigger, "END;");
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  /* deleting the old UPDATE trigger TYPE [if any] */
	  sprintf (xname, "ggu_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  /* inserting the new UPDATE trigger TYPE */
	  sprintf (trigger, "CREATE TRIGGER %s BEFORE UPDATE ON %s\n", xname,
		   xtable);
	  strcat (trigger, "FOR EACH ROW BEGIN\n");
	  sprintf (dummy,
		   "SELECT RAISE(ROLLBACK, '%s.%s violates Geometry constraint [geom-type or SRID not allowed]')\n",
		   sqltable, sqlcolumn);
	  strcat (trigger, dummy);
	  strcat (trigger, "WHERE (SELECT type FROM geometry_columns\n");
	  sprintf (dummy,
		   "WHERE f_table_name = '%s' AND f_geometry_column = '%s'\n",
		   sqltable, sqlcolumn);
	  strcat (trigger, dummy);
	  sprintf (dummy,
		   "AND GeometryConstraints(NEW.%s, type, srid, '%s') = 1) IS NULL;\n",
		   xcolname, txt_dims);
	  strcat (trigger, dummy);
	  strcat (trigger, "END;");
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  /* inserting SpatialIndex information into the linked list */
	  curr_idx = malloc (sizeof (struct spatial_index_str));
	  len = strlen (tblname);
	  curr_idx->TableName = malloc (len + 1);
	  strcpy (curr_idx->TableName, tblname);
	  len = strlen ((char *) colname);
	  curr_idx->ColumnName = malloc (len + 1);
	  strcpy (curr_idx->ColumnName, (char *) colname);
	  curr_idx->ValidRtree = (char) index;
	  curr_idx->ValidCache = (char) cached;
	  curr_idx->Next = NULL;
	  if (!first_idx)
	      first_idx = curr_idx;
	  if (last_idx)
	      last_idx->Next = curr_idx;
	  last_idx = curr_idx;
	  /* deleting the old INSERT trigger SPATIAL_INDEX [if any] */
	  sprintf (xname, "gii_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (index)
	    {
		/* inserting the new INSERT trigger SRID */
		sprintf (xindex, "idx_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER INSERT ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy, "DELETE FROM %s WHERE pkid=NEW.ROWID;\n",
			 xindex);
		strcat (trigger, dummy);
		sprintf (xindex, "idx_%s_%s", tblname, colname);
		clean_sql_string (xindex);
		sprintf (dummy, "SELECT RTreeAlign('%s', NEW.ROWID, NEW.%s);",
			 xindex, xcolname);
		strcat (trigger, dummy);
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  /* deleting the old UPDATE trigger SPATIAL_INDEX [if any] */
	  sprintf (xname, "giu_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (index)
	    {
		/* inserting the new UPDATE trigger SRID */
		sprintf (xindex, "idx_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER UPDATE ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy, "DELETE FROM %s WHERE pkid=NEW.ROWID;\n",
			 xindex);
		strcat (trigger, dummy);
		sprintf (xindex, "idx_%s_%s", tblname, colname);
		clean_sql_string (xindex);
		sprintf (dummy, "SELECT RTreeAlign('%s', NEW.ROWID, NEW.%s);",
			 xindex, xcolname);
		strcat (trigger, dummy);
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  /* deleting the old UPDATE trigger SPATIAL_INDEX [if any] */
	  sprintf (xname, "gid_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (index)
	    {
		/* inserting the new DELETE trigger SRID */
		sprintf (xindex, "idx_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER DELETE ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy, "DELETE FROM %s WHERE pkid = OLD.ROWID;\n",
			 xindex);
		strcat (trigger, dummy);
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  /* deleting the old INSERT trigger MBR_CACHE [if any] */
	  sprintf (xname, "gci_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (cached)
	    {
		/* inserting the new INSERT trigger SRID */
		sprintf (xindex, "cache_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER INSERT ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy,
			 "INSERT INTO %s (rowid, mbr) VALUES (NEW.ROWID,\nBuildMbrFilter(",
			 xindex);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMinX(NEW.%s), ", xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMinY(NEW.%s), ", xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMaxX(NEW.%s), ", xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMaxY(NEW.%s)));\n", xcolname);
		strcat (trigger, dummy);
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  /* deleting the old UPDATE trigger MBR_CACHE [if any] */
	  sprintf (xname, "gcu_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (cached)
	    {
		/* inserting the new UPDATE trigger SRID */
		sprintf (xindex, "cache_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER UPDATE ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy, "UPDATE %s SET ", xindex);
		strcat (trigger, dummy);
		sprintf (dummy, "mbr = BuildMbrFilter(MbrMinX(NEW.%s), ",
			 xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMinY(NEW.%s), ", xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMaxX(NEW.%s), ", xcolname);
		strcat (trigger, dummy);
		sprintf (dummy, "MbrMaxY(NEW.%s))\n", xcolname);
		strcat (trigger, dummy);
		strcat (trigger, "WHERE rowid = NEW.ROWID;\n");
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  /* deleting the old UPDATE trigger MBR_CACHE [if any] */
	  sprintf (xname, "gcd_%s_%s", tblname, colname);
	  double_quoted_sql (xname);
	  sprintf (trigger, "DROP TRIGGER IF EXISTS %s", xname);
	  ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
	  if (ret != SQLITE_OK)
	      goto error;
	  if (cached)
	    {
		/* inserting the new DELETE trigger SRID */
		sprintf (xindex, "cache_%s_%s", tblname, colname);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE TRIGGER %s AFTER DELETE ON %s\n",
			 xname, xtable);
		strcat (trigger, "FOR EACH ROW BEGIN\n");
		sprintf (dummy, "DELETE FROM %s WHERE rowid = OLD.ROWID;\n",
			 xindex);
		strcat (trigger, dummy);
		strcat (trigger, "END;");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
      }
    sqlite3_free_table (results);
/* now we'll adjust any related SpatialIndex as required */
    curr_idx = first_idx;
    while (curr_idx)
      {
	  if (curr_idx->ValidRtree)
	    {
		/* building RTree SpatialIndex */
		sprintf (xindex, "idx_%s_%s", curr_idx->TableName,
			 curr_idx->ColumnName);
		double_quoted_sql (xindex);
		sprintf (trigger, "CREATE VIRTUAL TABLE %s USING rtree(\n",
			 xindex);
		strcat (trigger, "pkid, xmin, xmax, ymin, ymax)");
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
		buildSpatialIndex (sqlite,
				   (unsigned char *) (curr_idx->TableName),
				   curr_idx->ColumnName);
	    }
	  if (curr_idx->ValidCache)
	    {
		/* building MbrCache SpatialIndex */
		sprintf (xindex, "cache_%s_%s", curr_idx->TableName,
			 curr_idx->ColumnName);
		double_quoted_sql (xindex);
		strcpy (xtable, curr_idx->TableName);
		double_quoted_sql (xtable);
		strcpy (xcolname, curr_idx->ColumnName);
		double_quoted_sql (xcolname);
		sprintf (trigger,
			 "CREATE VIRTUAL TABLE %s USING MbrCache(%s, %s)\n",
			 xindex, xtable, xcolname);
		ret = sqlite3_exec (sqlite, trigger, NULL, NULL, &errMsg);
		if (ret != SQLITE_OK)
		    goto error;
	    }
	  curr_idx = curr_idx->Next;
      }
    goto index_cleanup;
  error:
    fprintf (stderr, "updateTableTriggers: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
  index_cleanup:
    curr_idx = first_idx;
    while (curr_idx)
      {
	  next_idx = curr_idx->Next;
	  if (curr_idx->TableName)
	      free (curr_idx->TableName);
	  if (curr_idx->ColumnName)
	      free (curr_idx->ColumnName);
	  free (curr_idx);
	  curr_idx = next_idx;
      }
}

static void
fnct_AddGeometryColumn (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ AddGeometryColumn(table, column, srid, type , dimension  [  , not-null ]  )
/
/ creates a new COLUMN of given TYPE into TABLE
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    const unsigned char *type;
    const unsigned char *txt_dims;
    int xtype;
    int srid = -1;
    int dimension = 2;
    int dims = -1;
    char dummy[32];
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    int i;
    char tblname[256];
    char xtable[1024];
    char xcolumn[1024];
    char sqltable[1024];
    char sqlcolumn[1024];
    int notNull = 0;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    if (sqlite3_value_type (argv[2]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 3 [SRID] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    srid = sqlite3_value_int (argv[2]);
    if (sqlite3_value_type (argv[3]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 4 [geometry_type] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    type = sqlite3_value_text (argv[3]);
    if (sqlite3_value_type (argv[4]) == SQLITE_INTEGER)
      {
	  dimension = sqlite3_value_int (argv[4]);
	  if (dimension == 2)
	      dims = GAIA_XY;
	  if (dimension == 3)
	      dims = GAIA_XY_Z;
      }
    else if (sqlite3_value_type (argv[4]) == SQLITE_TEXT)
      {
	  txt_dims = sqlite3_value_text (argv[4]);
	  if (strcasecmp ((char *) txt_dims, "XY") == 0)
	      dims = GAIA_XY;
	  if (strcasecmp ((char *) txt_dims, "XYZ") == 0)
	      dims = GAIA_XY_Z;
	  if (strcasecmp ((char *) txt_dims, "XYM") == 0)
	      dims = GAIA_XY_M;
	  if (strcasecmp ((char *) txt_dims, "XYZM") == 0)
	      dims = GAIA_XY_Z_M;
      }
    else
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 5 [dimension] is not of the Integer or Text type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (argc > 5)
      {
	  /* optional NOT NULL arg */
	  if (sqlite3_value_type (argv[5]) != SQLITE_INTEGER)
	    {
		fprintf (stderr,
			 "AddGeometryColumn() error: argument 6 [not null] is not of the Integer type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  notNull = sqlite3_value_int (argv[5]);
      }
    xtype = GAIA_UNKNOWN;
    if (strcasecmp ((char *) type, "POINT") == 0)
	xtype = GAIA_POINT;
    if (strcasecmp ((char *) type, "LINESTRING") == 0)
	xtype = GAIA_LINESTRING;
    if (strcasecmp ((char *) type, "POLYGON") == 0)
	xtype = GAIA_POLYGON;
    if (strcasecmp ((char *) type, "MULTIPOINT") == 0)
	xtype = GAIA_MULTIPOINT;
    if (strcasecmp ((char *) type, "MULTILINESTRING") == 0)
	xtype = GAIA_MULTILINESTRING;
    if (strcasecmp ((char *) type, "MULTIPOLYGON") == 0)
	xtype = GAIA_MULTIPOLYGON;
    if (strcasecmp ((char *) type, "GEOMETRYCOLLECTION") == 0)
	xtype = GAIA_GEOMETRYCOLLECTION;
    if (strcasecmp ((char *) type, "GEOMETRY") == 0)
	xtype = -1;
    if (xtype == GAIA_UNKNOWN)
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 4 [geometry_type] has an illegal value\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (dims == GAIA_XY || dims == GAIA_XY_Z || dims == GAIA_XY_M
	|| dims == GAIA_XY_Z_M)
	;
    else
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: argument 5 [dimension] ILLEGAL VALUE\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
/* checking if the table exists */
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE '%s'",
	     sqltable);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    *tblname = '\0';
    for (i = 1; i <= rows; i++)
	strcpy (tblname, results[(i * columns)]);
    sqlite3_free_table (results);
    if (*tblname == '\0')
      {
	  fprintf (stderr,
		   "AddGeometryColumn() error: table '%s' does not exist\n",
		   table);
	  sqlite3_result_int (context, 0);
	  return;
      }
/* trying to add the column */
    strcpy (xtable, (char *) table);
    double_quoted_sql (xtable);
    strcpy (xcolumn, (char *) column);
    double_quoted_sql (xcolumn);
    strcpy (sql, "ALTER TABLE ");
    strcat (sql, xtable);
    strcat (sql, " ADD COLUMN ");
    strcat (sql, xcolumn);
    strcat (sql, " ");
    switch (xtype)
      {
      case GAIA_POINT:
	  strcat (sql, "POINT");
	  break;
      case GAIA_LINESTRING:
	  strcat (sql, "LINESTRING");
	  break;
      case GAIA_POLYGON:
	  strcat (sql, "POLYGON");
	  break;
      case GAIA_MULTIPOINT:
	  strcat (sql, "MULTIPOINT");
	  break;
      case GAIA_MULTILINESTRING:
	  strcat (sql, "MULTILINESTRING");
	  break;
      case GAIA_MULTIPOLYGON:
	  strcat (sql, "MULTIPOLYGON");
	  break;
      case GAIA_GEOMETRYCOLLECTION:
	  strcat (sql, "GEOMETRYCOLLECTION");
	  break;
      case -1:
	  strcat (sql, "GEOMETRY");
	  break;
      };
    if (notNull)
      {
	  /* adding a NOT NULL clause */
	  strcat (sql, " NOT NULL DEFAULT ''");
      }
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/*ok, inserting into geometry_columns [Spatial Metadata] */
    strcpy (sql,
	    "INSERT INTO geometry_columns (f_table_name, f_geometry_column, type, ");
    strcat (sql, "coord_dimension, srid, spatial_index_enabled) VALUES (");
    strcat (sql, "'");
    strcat (sql, sqltable);
    strcat (sql, "', '");
    strcat (sql, sqlcolumn);
    strcat (sql, "', '");
    switch (xtype)
      {
      case GAIA_POINT:
	  strcat (sql, "POINT");
	  break;
      case GAIA_LINESTRING:
	  strcat (sql, "LINESTRING");
	  break;
      case GAIA_POLYGON:
	  strcat (sql, "POLYGON");
	  break;
      case GAIA_MULTIPOINT:
	  strcat (sql, "MULTIPOINT");
	  break;
      case GAIA_MULTILINESTRING:
	  strcat (sql, "MULTILINESTRING");
	  break;
      case GAIA_MULTIPOLYGON:
	  strcat (sql, "MULTIPOLYGON");
	  break;
      case GAIA_GEOMETRYCOLLECTION:
	  strcat (sql, "GEOMETRYCOLLECTION");
	  break;
      case -1:
	  strcat (sql, "GEOMETRY");
	  break;
      };
    strcat (sql, "', '");
    switch (dims)
      {
      case GAIA_XY:
	  strcat (sql, "XY");
	  break;
      case GAIA_XY_Z:
	  strcat (sql, "XYZ");
	  break;
      case GAIA_XY_M:
	  strcat (sql, "XYM");
	  break;
      case GAIA_XY_Z_M:
	  strcat (sql, "XYZM");
	  break;
      };
    strcat (sql, "', ");
    if (srid <= 0)
	strcat (sql, "-1");
    else
      {
	  sprintf (dummy, "%d", srid);
	  strcat (sql, dummy);
      }
    strcat (sql, ", 0)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "Geometry [");
    switch (xtype)
      {
      case GAIA_POINT:
	  strcat (sql, "POINT");
	  break;
      case GAIA_LINESTRING:
	  strcat (sql, "LINESTRING");
	  break;
      case GAIA_POLYGON:
	  strcat (sql, "POLYGON");
	  break;
      case GAIA_MULTIPOINT:
	  strcat (sql, "MULTIPOINT");
	  break;
      case GAIA_MULTILINESTRING:
	  strcat (sql, "MULTILINESTRING");
	  break;
      case GAIA_MULTIPOLYGON:
	  strcat (sql, "MULTIPOLYGON");
	  break;
      case GAIA_GEOMETRYCOLLECTION:
	  strcat (sql, "GEOMETRYCOLLECTION");
	  break;
      case -1:
	  strcat (sql, "GEOMETRY");
	  break;
      };
    strcat (sql, ",");
    switch (dims)
      {
      case GAIA_XY:
	  strcat (sql, "XY");
	  break;
      case GAIA_XY_Z:
	  strcat (sql, "XYZ");
	  break;
      case GAIA_XY_M:
	  strcat (sql, "XYM");
	  break;
      case GAIA_XY_Z_M:
	  strcat (sql, "XYZM");
	  break;
      };
    sprintf (sqlcolumn, ",SRID=%d", (srid <= 0) ? -1 : srid);
    strcat (sql, sqlcolumn);
    strcat (sql, "] successfully created");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "AddGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_RecoverGeometryColumn (sqlite3_context * context, int argc,
			    sqlite3_value ** argv)
{
/* SQL function:
/ RecoverGeometryColumn(table, column, srid, type , dimension )
/
/ checks if an existing TABLE.COLUMN satisfies the required geometric features
/ if yes adds it to SpatialMetaData and enabling triggers
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    const unsigned char *type;
    int xtype;
    int xxtype;
    int srid = -1;
    const unsigned char *txt_dims;
    int dimension = 2;
    int dims = -1;
    char dummy[32];
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    int i;
    char tblname[256];
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    if (sqlite3_value_type (argv[2]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 3 [SRID] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    srid = sqlite3_value_int (argv[2]);
    if (sqlite3_value_type (argv[3]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 4 [geometry_type] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    type = sqlite3_value_text (argv[3]);
    if (sqlite3_value_type (argv[4]) == SQLITE_INTEGER)
      {
	  dimension = sqlite3_value_int (argv[4]);
	  if (dimension == 2)
	      dims = GAIA_XY;
	  if (dimension == 3)
	      dims = GAIA_XY_Z;
      }
    else if (sqlite3_value_type (argv[4]) == SQLITE_TEXT)
      {
	  txt_dims = sqlite3_value_text (argv[4]);
	  if (strcasecmp ((char *) txt_dims, "XY") == 0)
	      dims = GAIA_XY;
	  if (strcasecmp ((char *) txt_dims, "XYZ") == 0)
	      dims = GAIA_XY_Z;
	  if (strcasecmp ((char *) txt_dims, "XYM") == 0)
	      dims = GAIA_XY_M;
	  if (strcasecmp ((char *) txt_dims, "XYZM") == 0)
	      dims = GAIA_XY_Z_M;
      }
    else
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 5 [dimension] is not of the Integer or Text type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    xtype = GAIA_UNKNOWN;
    if (strcasecmp ((char *) type, "POINT") == 0)
	xtype = GAIA_POINT;
    if (strcasecmp ((char *) type, "LINESTRING") == 0)
	xtype = GAIA_LINESTRING;
    if (strcasecmp ((char *) type, "POLYGON") == 0)
	xtype = GAIA_POLYGON;
    if (strcasecmp ((char *) type, "MULTIPOINT") == 0)
	xtype = GAIA_MULTIPOINT;
    if (strcasecmp ((char *) type, "MULTILINESTRING") == 0)
	xtype = GAIA_MULTILINESTRING;
    if (strcasecmp ((char *) type, "MULTIPOLYGON") == 0)
	xtype = GAIA_MULTIPOLYGON;
    if (strcasecmp ((char *) type, "GEOMETRYCOLLECTION") == 0)
	xtype = GAIA_GEOMETRYCOLLECTION;
    if (strcasecmp ((char *) type, "GEOMETRY") == 0)
	xtype = -1;
    if (xtype == GAIA_UNKNOWN)
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 4 [geometry_type] has an illegal value\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (dims == GAIA_XY || dims == GAIA_XY_Z || dims == GAIA_XY_M
	|| dims == GAIA_XY_Z_M)
	;
    else
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: argument 5 [dimension] ILLEGAL VALUE\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
/* checking if the table exists */
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE '%s'",
	     sqltable);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "RecoverGeometryColumn: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    *tblname = '\0';
    for (i = 1; i <= rows; i++)
      {
	  /* preparing the triggers */
	  strcpy (tblname, results[(i * columns)]);
      }
    sqlite3_free_table (results);
    if (*tblname == '\0')
      {
	  fprintf (stderr,
		   "RecoverGeometryColumn() error: table '%s' does not exist\n",
		   table);
	  sqlite3_result_int (context, 0);
	  return;
      }
/* adjusting the actual GeometryType */
    xxtype = xtype;
    xtype = GAIA_UNKNOWN;
    if (xxtype == GAIA_POINT)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_POINTZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_POINTM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_POINTZM;
		break;
	    default:
		xtype = GAIA_POINT;
		break;
	    };
      }
    if (xxtype == GAIA_LINESTRING)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_LINESTRINGZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_LINESTRINGM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_LINESTRINGZM;
		break;
	    default:
		xtype = GAIA_LINESTRING;
		break;
	    };
      }
    if (xxtype == GAIA_POLYGON)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_POLYGONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_POLYGONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_POLYGONZM;
		break;
	    default:
		xtype = GAIA_POLYGON;
		break;
	    };
      }
    if (xxtype == GAIA_MULTIPOINT)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTIPOINTZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTIPOINTM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTIPOINTZM;
		break;
	    default:
		xtype = GAIA_MULTIPOINT;
		break;
	    };
      }
    if (xxtype == GAIA_MULTILINESTRING)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTILINESTRINGZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTILINESTRINGM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTILINESTRINGZM;
		break;
	    default:
		xtype = GAIA_MULTILINESTRING;
		break;
	    };
      }
    if (xxtype == GAIA_MULTIPOLYGON)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_MULTIPOLYGONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_MULTIPOLYGONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_MULTIPOLYGONZM;
		break;
	    default:
		xtype = GAIA_MULTIPOLYGON;
		break;
	    };
      }
    if (xxtype == GAIA_GEOMETRYCOLLECTION)
      {
	  switch (dims)
	    {
	    case GAIA_XY_Z:
		xtype = GAIA_GEOMETRYCOLLECTIONZ;
		break;
	    case GAIA_XY_M:
		xtype = GAIA_GEOMETRYCOLLECTIONM;
		break;
	    case GAIA_XY_Z_M:
		xtype = GAIA_GEOMETRYCOLLECTIONZM;
		break;
	    default:
		xtype = GAIA_GEOMETRYCOLLECTION;
		break;
	    };
      }
    if (xxtype == -1)
	xtype = -1;		/* GEOMETRY */
    if (!recoverGeomColumn (sqlite, table, column, xtype, dims, srid))
      {
	  fprintf (stderr, "RecoverGeometryColumn(): validation failed\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    strcpy (sql,
	    "INSERT INTO geometry_columns (f_table_name, f_geometry_column, type, ");
    strcat (sql, "coord_dimension, srid, spatial_index_enabled) VALUES (");
    strcat (sql, "'");
    strcat (sql, sqltable);
    strcat (sql, "', '");
    strcat (sql, sqlcolumn);
    strcat (sql, "', '");
    switch (xtype)
      {
      case GAIA_POINT:
      case GAIA_POINTZ:
      case GAIA_POINTM:
      case GAIA_POINTZM:
	  strcat (sql, "POINT");
	  break;
      case GAIA_LINESTRING:
      case GAIA_LINESTRINGZ:
      case GAIA_LINESTRINGM:
      case GAIA_LINESTRINGZM:
	  strcat (sql, "LINESTRING");
	  break;
      case GAIA_POLYGON:
      case GAIA_POLYGONZ:
      case GAIA_POLYGONM:
      case GAIA_POLYGONZM:
	  strcat (sql, "POLYGON");
	  break;
      case GAIA_MULTIPOINT:
      case GAIA_MULTIPOINTZ:
      case GAIA_MULTIPOINTM:
      case GAIA_MULTIPOINTZM:
	  strcat (sql, "MULTIPOINT");
	  break;
      case GAIA_MULTILINESTRING:
      case GAIA_MULTILINESTRINGZ:
      case GAIA_MULTILINESTRINGM:
      case GAIA_MULTILINESTRINGZM:
	  strcat (sql, "MULTILINESTRING");
	  break;
      case GAIA_MULTIPOLYGON:
      case GAIA_MULTIPOLYGONZ:
      case GAIA_MULTIPOLYGONM:
      case GAIA_MULTIPOLYGONZM:
	  strcat (sql, "MULTIPOLYGON");
	  break;
      case GAIA_GEOMETRYCOLLECTION:
      case GAIA_GEOMETRYCOLLECTIONZ:
      case GAIA_GEOMETRYCOLLECTIONM:
      case GAIA_GEOMETRYCOLLECTIONZM:
	  strcat (sql, "GEOMETRYCOLLECTION");
	  break;
      case -1:
	  strcat (sql, "GEOMETRY");
	  break;
      };
    strcat (sql, "', '");
    switch (dims)
      {
      case GAIA_XY:
	  strcat (sql, "XY");
	  break;
      case GAIA_XY_Z:
	  strcat (sql, "XYZ");
	  break;
      case GAIA_XY_M:
	  strcat (sql, "XYM");
	  break;
      case GAIA_XY_Z_M:
	  strcat (sql, "XYZM");
	  break;
      };
    strcat (sql, "', ");
    if (srid <= 0)
	strcat (sql, "-1");
    else
      {
	  sprintf (dummy, "%d", srid);
	  strcat (sql, dummy);
      }
    strcat (sql, ", 0)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "Geometry [");
    switch (xtype)
      {
      case GAIA_POINT:
	  strcat (sql, "POINT");
	  break;
      case GAIA_LINESTRING:
	  strcat (sql, "LINESTRING");
	  break;
      case GAIA_POLYGON:
	  strcat (sql, "POLYGON");
	  break;
      case GAIA_MULTIPOINT:
	  strcat (sql, "MULTIPOINT");
	  break;
      case GAIA_MULTILINESTRING:
	  strcat (sql, "MULTILINESTRING");
	  break;
      case GAIA_MULTIPOLYGON:
	  strcat (sql, "MULTIPOLYGON");
	  break;
      case GAIA_GEOMETRYCOLLECTION:
	  strcat (sql, "GEOMETRYCOLLECTION");
	  break;
      case -1:
	  strcat (sql, "GEOMETRY");
	  break;
      };
    strcat (sql, ",");
    switch (dims)
      {
      case GAIA_XY:
	  strcat (sql, "XY");
	  break;
      case GAIA_XY_Z:
	  strcat (sql, "XYZ");
	  break;
      case GAIA_XY_M:
	  strcat (sql, "XYM");
	  break;
      case GAIA_XY_Z_M:
	  strcat (sql, "XYZM");
	  break;
      };
    sprintf (sqlcolumn, ",SRID=%d", (srid <= 0) ? -1 : srid);
    strcat (sql, sqlcolumn);
    strcat (sql, "] successfully recovered");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "RecoverGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_DiscardGeometryColumn (sqlite3_context * context, int argc,
			    sqlite3_value ** argv)
{
/* SQL function:
/ DiscardGeometryColumn(table, column)
/
/ removes TABLE.COLUMN from the Spatial MetaData [thus disabling triggers too]
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char xname[1024];
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DiscardGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DiscardGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "DELETE FROM geometry_columns WHERE f_table_name LIKE '%s' AND f_geometry_column LIKE '%s'",
	     sqltable, sqlcolumn);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/* removing triggers too */
    sprintf (xname, "ggi_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "ggu_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gii_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "giu_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gid_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gci_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gcu_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gcd_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;

    /* trying to delete old versions [v2.0, v2.2] triggers[if any] */
    sprintf (xname, "gti_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gtu_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gsi_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sprintf (xname, "gsu_%s_%s", (char *) table, (char *) column);
    double_quoted_sql (xname);
    sprintf (sql, "DROP TRIGGER IF EXISTS %s", xname);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    /* end deletion old versions [v2.0, v2.2] triggers[if any] */

    sqlite3_result_int (context, 1);
    strcpy (sql, "Geometry successfully discarded");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "DiscardGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_InitFDOSpatialMetaData (sqlite3_context * context, int argc,
			     sqlite3_value ** argv)
{
/* SQL function:
/ InitFDOSpatialMetaData(void)
/
/ creates the FDO-styled SPATIAL_REF_SYS and GEOMETRY_COLUMNS tables
/ returns 1 on success
/ 0 on failure
*/
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
/* creating the SPATIAL_REF_SYS tables */
    strcpy (sql, "CREATE TABLE spatial_ref_sys (\n");
    strcat (sql, "srid INTEGER PRIMARY KEY,\n");
    strcat (sql, "auth_name TEXT,\n");
    strcat (sql, "auth_srid INTEGER,\n");
    strcat (sql, "srtext TEXT)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/* creating the GEOMETRY_COLUMN tables */
    strcpy (sql, "CREATE TABLE geometry_columns (\n");
    strcat (sql, "f_table_name TEXT,\n");
    strcat (sql, "f_geometry_column TEXT,\n");
    strcat (sql, "geometry_type INTEGER,\n");
    strcat (sql, "coord_dimension INTEGER,\n");
    strcat (sql, "srid INTEGER,\n");
    strcat (sql, "geometry_format TEXT)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sqlite3_result_int (context, 1);
    return;
  error:fprintf (stderr, "InitFDOSpatiaMetaData() error: \"%s\"\n",
	     errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static int
recoverFDOGeomColumn (sqlite3 * sqlite, const unsigned char *table,
		      const unsigned char *column, int xtype, int srid)
{
/* checks if TABLE.COLUMN exists and has the required features */
    int ok = 1;
    char sql[1024];
    int type;
    sqlite3_stmt *stmt;
    gaiaGeomCollPtr geom;
    const void *blob_value;
    int len;
    int ret;
    int i_col;
    char xcolumn[1024];
    char xtable[1024];
    strcpy (xcolumn, (char *) column);
    double_quoted_sql (xcolumn);
    strcpy (xtable, (char *) table);
    double_quoted_sql (xtable);
    sprintf (sql, "SELECT %s FROM %s", xcolumn, xtable);
/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "recoverFDOGeomColumn: error %d \"%s\"\n",
		   sqlite3_errcode (sqlite), sqlite3_errmsg (sqlite));
	  return 0;
      }
    while (1)
      {
	  /* scrolling the result set rows */
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;		/* end of result set */
	  if (ret == SQLITE_ROW)
	    {
		/* checking Geometry features */
		geom = NULL;
		for (i_col = 0; i_col < sqlite3_column_count (stmt); i_col++)
		  {
		      if (sqlite3_column_type (stmt, i_col) != SQLITE_BLOB)
			  ok = 0;
		      else
			{
			    blob_value = sqlite3_column_blob (stmt, i_col);
			    len = sqlite3_column_bytes (stmt, i_col);
			    geom = gaiaFromSpatiaLiteBlobWkb (blob_value, len);
			    if (!geom)
				ok = 0;
			    else
			      {
				  if (geom->Srid != srid)
				      ok = 0;
				  type = gaiaGeometryType (geom);
				  if (xtype == type)
				      ;
				  else
				      ok = 0;
				  gaiaFreeGeomColl (geom);
			      }
			}
		  }
	    }
	  if (!ok)
	      break;
      }
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "recoverFDOGeomColumn: error %d \"%s\"\n",
		   sqlite3_errcode (sqlite), sqlite3_errmsg (sqlite));
	  return 0;
      }
    return ok;
}

static void
fnct_AddFDOGeometryColumn (sqlite3_context * context, int argc,
			   sqlite3_value ** argv)
{
/* SQL function:
/ AddFDOGeometryColumn(table, column, srid, geometry_type , dimension, geometry_format )
/
/ creates a new COLUMN of given TYPE into TABLE
/ returns 1 on success
/ 0 on failure
*/
    const char *table;
    const char *column;
    const char *format;
    char xformat[64];
    int type;
    int srid = -1;
    int dimension = 2;
    char dummy[32];
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    int i;
    char tblname[256];
    char xtable[1024];
    char xcolumn[1024];
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = (const char *) sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = (const char *) sqlite3_value_text (argv[1]);
    if (sqlite3_value_type (argv[2]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 3 [SRID] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    srid = sqlite3_value_int (argv[2]);
    if (sqlite3_value_type (argv[3]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 4 [geometry_type] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    type = sqlite3_value_int (argv[3]);
    if (sqlite3_value_type (argv[4]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 5 [dimension] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    dimension = sqlite3_value_int (argv[4]);
    if (sqlite3_value_type (argv[5]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 6 [geometry_format] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    format = (const char *) sqlite3_value_text (argv[5]);
    if (type ==
	GAIA_POINT
	|| type ==
	GAIA_LINESTRING
	|| type ==
	GAIA_POLYGON
	|| type ==
	GAIA_MULTIPOINT
	|| type ==
	GAIA_MULTILINESTRING
	|| type == GAIA_MULTIPOLYGON || type == GAIA_GEOMETRYCOLLECTION)
	;
    else
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 4 [geometry_type] has an illegal value\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (dimension < 2 || dimension > 4)
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 5 [dimension] current version only accepts dimension=2,3,4\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (strcasecmp (format, "WKT") == 0)
	strcpy (xformat, "WKT");
    else if (strcasecmp (format, "WKB") == 0)
	strcpy (xformat, "WKB");
    else if (strcasecmp (format, "FGF") == 0)
	strcpy (xformat, "FGF");
    else
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: argument 6 [geometry_format] has to be one of: WKT,WKB,FGF\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
/* checking if the table exists */
    strcpy (xtable, (char *) table);
    double_quoted_sql (xtable);
    strcpy (xcolumn, (char *) column);
    double_quoted_sql (xcolumn);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE '%s'",
	     sqltable);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddFDOGeometryColumn: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    *tblname = '\0';
    for (i = 1; i <= rows; i++)
      {
	  strcpy (tblname, results[(i * columns)]);
      }
    sqlite3_free_table (results);
    if (*tblname == '\0')
      {
	  fprintf (stderr,
		   "AddFDOGeometryColumn() error: table '%s' does not exist\n",
		   table);
	  sqlite3_result_int (context, 0);
	  return;
      }
/* trying to add the column */
    strcpy (sql, "ALTER TABLE ");
    strcat (sql, xtable);
    strcat (sql, " ADD COLUMN ");
    strcat (sql, xcolumn);
    strcat (sql, " BLOB");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/*ok, inserting into geometry_columns [FDO Spatial Metadata] */
    strcpy (sql,
	    "INSERT INTO geometry_columns (f_table_name, f_geometry_column, geometry_type, ");
    strcat (sql, "coord_dimension, srid, geometry_format) VALUES (");
    strcat (sql, "'");
    strcat (sql, sqltable);
    strcat (sql, "', '");
    strcat (sql, sqlcolumn);
    strcat (sql, "', ");
    sprintf (dummy, "%d, %d, ", type, dimension);
    strcat (sql, dummy);
    if (srid <= 0)
	strcat (sql, "-1");
    else
      {
	  sprintf (dummy, "%d", srid);
	  strcat (sql, dummy);
      }
    strcat (sql, ", '");
    strcat (sql, xformat);
    strcat (sql, "')");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sqlite3_result_int (context, 1);
    return;
  error:
    fprintf (stderr, "AddFDOGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_RecoverFDOGeometryColumn (sqlite3_context * context, int argc,
			       sqlite3_value ** argv)
{
/* SQL function:
/ RecoverFDOGeometryColumn(table, column, srid, geometry_type , dimension, geometry_format )
/
/ checks if an existing TABLE.COLUMN satisfies the required geometric features
/ if yes adds it to FDO-styled SpatialMetaData 
/ returns 1 on success
/ 0 on failure
*/
    const char *table;
    const char *column;
    const char *format;
    char xformat[64];
    int type;
    int srid = -1;
    int dimension = 2;
    char dummy[32];
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    int i;
    char tblname[256];
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = (const char *) sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = (const char *) sqlite3_value_text (argv[1]);
    if (sqlite3_value_type (argv[2]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 3 [SRID] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    srid = sqlite3_value_int (argv[2]);
    if (sqlite3_value_type (argv[3]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 4 [geometry_type] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    type = sqlite3_value_int (argv[3]);
    if (sqlite3_value_type (argv[4]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 5 [dimension] is not of the Integer type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    dimension = sqlite3_value_int (argv[4]);
    if (sqlite3_value_type (argv[5]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 6 [geometry_format] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    format = (const char *) sqlite3_value_text (argv[5]);
    if (type ==
	GAIA_POINT
	|| type ==
	GAIA_LINESTRING
	|| type ==
	GAIA_POLYGON
	|| type ==
	GAIA_MULTIPOINT
	|| type ==
	GAIA_MULTILINESTRING
	|| type == GAIA_MULTIPOLYGON || type == GAIA_GEOMETRYCOLLECTION)
	;
    else
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 4 [geometry_type] has an illegal value\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (dimension < 2 || dimension > 4)
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 5 [dimension] current version only accepts dimension=2,3,4\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (strcasecmp (format, "WKT") == 0)
	strcpy (xformat, "WKT");
    else if (strcasecmp (format, "WKB") == 0)
	strcpy (xformat, "WKB");
    else if (strcasecmp (format, "FGF") == 0)
	strcpy (xformat, "FGF");
    else
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: argument 6 [geometry_format] has to be one of: WKT,WKB,FGF\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
/* checking if the table exists */
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "SELECT name FROM sqlite_master WHERE type = 'table' AND name LIKE '%s'",
	     sqltable);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "RecoverFDOGeometryColumn: \"%s\"\n", errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    *tblname = '\0';
    for (i = 1; i <= rows; i++)
      {
	  strcpy (tblname, results[(i * columns)]);
      }
    sqlite3_free_table (results);
    if (*tblname == '\0')
      {
	  fprintf (stderr,
		   "RecoverFDOGeometryColumn() error: table '%s' does not exist\n",
		   table);
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (!recoverFDOGeomColumn
	(sqlite, (const unsigned char *) table,
	 (const unsigned char *) column, type, srid))
      {
	  fprintf (stderr, "RecoverFDOGeometryColumn(): validation failed\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    strcpy (sqltable, (char *) tblname);
    clean_sql_string (sqltable);
    strcpy (sql,
	    "INSERT INTO geometry_columns (f_table_name, f_geometry_column, geometry_type, ");
    strcat (sql, "coord_dimension, srid, geometry_format) VALUES (");
    strcat (sql, "'");
    strcat (sql, sqltable);
    strcat (sql, "', '");
    strcat (sql, sqlcolumn);
    strcat (sql, "', ");
    sprintf (dummy, "%d, %d, ", type, dimension);
    strcat (sql, dummy);
    if (srid <= 0)
	strcat (sql, "-1");
    else
      {
	  sprintf (dummy, "%d", srid);
	  strcat (sql, dummy);
      }
    strcat (sql, ", '");
    strcat (sql, xformat);
    strcat (sql, "')");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sqlite3_result_int (context, 1);
    return;
  error:
    fprintf (stderr, "RecoverFDOGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_DiscardFDOGeometryColumn (sqlite3_context * context, int argc,
			       sqlite3_value ** argv)
{
/* SQL function:
/ DiscardFDOGeometryColumn(table, column)
/
/ removes TABLE.COLUMN from the Spatial MetaData
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DiscardFDOGeometryColumn() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DiscardFDOGeometryColumn() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    sprintf (sql,
	     "DELETE FROM geometry_columns WHERE f_table_name LIKE '%s' AND f_geometry_column LIKE '%s'",
	     sqltable, sqlcolumn);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    sqlite3_result_int (context, 1);
    return;
  error:
    fprintf (stderr, "DiscardFDOGeometryColumn() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static int
eval_rtree_entry (int ok_geom, double geom_value, int ok_rtree,
		  double rtree_value)
{
/* evaluating geom-coord and rtree-coord */
    if (!ok_geom && !ok_rtree)
	return 1;
    if (ok_geom && ok_rtree)
      {
	  float g = (float) geom_value;
	  float r = (float) rtree_value;
	  if (g != r)
	      return 0;
	  return 1;
      }
    return 0;
}

static int
check_spatial_index (sqlite3 * sqlite, const unsigned char *table,
		     const unsigned char *geom)
{
/* attempting to check an R*Tree for consistency */
    char xtable[1024];
    char xgeom[1024];
    char idx_name[2048];
    char sql[8192];
    char sql2[2048];
    int ret;
    int is_defined = 0;
    sqlite3_stmt *stmt;
    sqlite3_int64 count_geom;
    sqlite3_int64 count_rtree;
    double g_xmin;
    double g_ymin;
    double g_xmax;
    double g_ymax;
    int ok_g_xmin;
    int ok_g_ymin;
    int ok_g_xmax;
    int ok_g_ymax;
    double i_xmin;
    double i_ymin;
    double i_xmax;
    double i_ymax;
    int ok_i_xmin;
    int ok_i_ymin;
    int ok_i_xmax;
    int ok_i_ymax;

/* checking if the R*Tree Spatial Index is defined */
    strcpy (xtable, (const char *) table);
    clean_sql_string (xtable);
    strcpy (xgeom, (const char *) geom);
    clean_sql_string (xgeom);
    strcpy (sql, "SELECT Count(*) FROM geometry_columns ");
    sprintf (sql2, "WHERE f_table_name LIKE '%s' ", xtable);
    strcat (sql, sql2);
    sprintf (sql2, "AND f_geometry_column LIKE '%s' ", xgeom);
    strcat (sql, sql2);
    strcat (sql, "AND spatial_index_enabled = 1");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	      is_defined = sqlite3_column_int (stmt, 0);
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    if (!is_defined)
	return -1;

    sprintf (xgeom, "%s", geom);
    double_quoted_sql (xgeom);
    strcpy (xtable, (const char *) table);
    double_quoted_sql (xtable);
    sprintf (idx_name, "idx_%s_%s", table, geom);
    double_quoted_sql (idx_name);

/* counting how many Geometries are set into the main-table */
    sprintf (sql, "SELECT Count(*) FROM %s ", xtable);
    sprintf (sql2, "WHERE ST_GeometryType(%s) IS NOT NULL", xgeom);
    strcat (sql, sql2);
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	      count_geom = sqlite3_column_int (stmt, 0);
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);

/* counting how many R*Tree entries are defined */
    sprintf (sql, "SELECT Count(*) FROM %s", idx_name);
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	      count_rtree = sqlite3_column_int (stmt, 0);
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    if (count_geom != count_rtree)
      {
	  /* unexpected count difference */
	  return 0;
      }

/* checking the geometry-table against the corresponding R*Tree */
    sprintf (sql, "SELECT ");
    sprintf (sql2, "MbrMinX(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMinY(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMaxX(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMaxY(g.%s), ", xgeom);
    strcat (sql, sql2);
    strcat (sql, "i.xmin, i.ymin, i.xmax, i.ymax\n");
    sprintf (sql2, "FROM %s AS g\n", xtable);
    strcat (sql, sql2);
    sprintf (sql2, "LEFT JOIN %s AS i ", idx_name);
    strcat (sql, sql2);
    strcat (sql, "ON (g.ROWID = i.pkid)");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	    {
		/* checking a row */
		ok_g_xmin = 1;
		ok_g_ymin = 1;
		ok_g_xmax = 1;
		ok_g_ymax = 1;
		ok_i_xmin = 1;
		ok_i_ymin = 1;
		ok_i_xmax = 1;
		ok_i_ymax = 1;
		if (sqlite3_column_type (stmt, 0) == SQLITE_NULL)
		    ok_g_xmin = 0;
		else
		    g_xmin = sqlite3_column_double (stmt, 0);
		if (sqlite3_column_type (stmt, 1) == SQLITE_NULL)
		    ok_g_ymin = 0;
		else
		    g_ymin = sqlite3_column_double (stmt, 1);
		if (sqlite3_column_type (stmt, 2) == SQLITE_NULL)
		    ok_g_xmax = 0;
		else
		    g_xmax = sqlite3_column_double (stmt, 2);
		if (sqlite3_column_type (stmt, 3) == SQLITE_NULL)
		    ok_g_ymax = 0;
		else
		    g_ymax = sqlite3_column_double (stmt, 3);
		if (sqlite3_column_type (stmt, 4) == SQLITE_NULL)
		    ok_i_xmin = 0;
		else
		    i_xmin = sqlite3_column_double (stmt, 4);
		if (sqlite3_column_type (stmt, 5) == SQLITE_NULL)
		    ok_i_ymin = 0;
		else
		    i_ymin = sqlite3_column_double (stmt, 5);
		if (sqlite3_column_type (stmt, 6) == SQLITE_NULL)
		    ok_i_xmax = 0;
		else
		    i_xmax = sqlite3_column_double (stmt, 6);
		if (sqlite3_column_type (stmt, 7) == SQLITE_NULL)
		    ok_i_ymax = 0;
		else
		    i_ymax = sqlite3_column_double (stmt, 7);
		if (eval_rtree_entry (ok_g_xmin, g_xmin, ok_i_xmin, i_xmin) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_ymin, g_ymin, ok_i_ymin, i_ymin) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_xmax, g_xmax, ok_i_xmax, i_xmax) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_ymax, g_ymax, ok_i_ymax, i_ymax) ==
		    0)
		    goto mismatching;
	    }
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
/* we have now to finalize the query [memory cleanup] */
    sqlite3_finalize (stmt);


/* now we'll check the R*Tree against the corresponding geometry-table */
    sprintf (sql, "SELECT ");
    sprintf (sql2, "MbrMinX(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMinY(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMaxX(g.%s), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "MbrMaxY(g.%s), ", xgeom);
    strcat (sql, sql2);
    strcat (sql, "i.xmin, i.ymin, i.xmax, i.ymax\n");
    sprintf (sql2, "FROM %s AS i\n", idx_name);
    strcat (sql, sql2);
    sprintf (sql2, "LEFT JOIN %s AS g ", xtable);
    strcat (sql, sql2);
    strcat (sql, "ON (g.ROWID = i.pkid)");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	    {
		/* checking a row */
		ok_g_xmin = 1;
		ok_g_ymin = 1;
		ok_g_xmax = 1;
		ok_g_ymax = 1;
		ok_i_xmin = 1;
		ok_i_ymin = 1;
		ok_i_xmax = 1;
		ok_i_ymax = 1;
		if (sqlite3_column_type (stmt, 0) == SQLITE_NULL)
		    ok_g_xmin = 0;
		else
		    g_xmin = sqlite3_column_double (stmt, 0);
		if (sqlite3_column_type (stmt, 1) == SQLITE_NULL)
		    ok_g_ymin = 0;
		else
		    g_ymin = sqlite3_column_double (stmt, 1);
		if (sqlite3_column_type (stmt, 2) == SQLITE_NULL)
		    ok_g_xmax = 0;
		else
		    g_xmax = sqlite3_column_double (stmt, 2);
		if (sqlite3_column_type (stmt, 3) == SQLITE_NULL)
		    ok_g_ymax = 0;
		else
		    g_ymax = sqlite3_column_double (stmt, 3);
		if (sqlite3_column_type (stmt, 4) == SQLITE_NULL)
		    ok_i_xmin = 0;
		else
		    i_xmin = sqlite3_column_double (stmt, 4);
		if (sqlite3_column_type (stmt, 5) == SQLITE_NULL)
		    ok_i_ymin = 0;
		else
		    i_ymin = sqlite3_column_double (stmt, 5);
		if (sqlite3_column_type (stmt, 6) == SQLITE_NULL)
		    ok_i_xmax = 0;
		else
		    i_xmax = sqlite3_column_double (stmt, 6);
		if (sqlite3_column_type (stmt, 7) == SQLITE_NULL)
		    ok_i_ymax = 0;
		else
		    i_ymax = sqlite3_column_double (stmt, 7);
		if (eval_rtree_entry (ok_g_xmin, g_xmin, ok_i_xmin, i_xmin) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_ymin, g_ymin, ok_i_ymin, i_ymin) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_xmax, g_xmax, ok_i_xmax, i_xmax) ==
		    0)
		    goto mismatching;
		if (eval_rtree_entry (ok_g_ymax, g_ymax, ok_i_ymax, i_ymax) ==
		    0)
		    goto mismatching;
	    }
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    strcpy (sql, "Check SpatialIndex: is valid");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) geom, sql);
    return 1;
  mismatching:
    sqlite3_finalize (stmt);
    strcpy (sql, "Check SpatialIndex: INCONSISTENCIES detected");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) geom, sql);
    return 0;
}

static int
check_any_spatial_index (sqlite3 * sqlite)
{
/* attempting to check any defined R*Tree for consistency */
    const unsigned char *table;
    const unsigned char *column;
    int status;
    char sql[1024];
    int ret;
    int invalid_rtree = 0;
    sqlite3_stmt *stmt;

/* retrieving any defined R*Tree */
    strcpy (sql,
	    "SELECT f_table_name, f_geometry_column FROM geometry_columns ");
    strcat (sql, "WHERE spatial_index_enabled = 1");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CheckSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	    {
		/* checking a single R*Tree */
		table = sqlite3_column_text (stmt, 0);
		column = sqlite3_column_text (stmt, 1);
		status = check_spatial_index (sqlite, table, column);
		if (status < 0)
		  {
		      sqlite3_finalize (stmt);
		      return -1;
		  }
		if (status == 0)
		    invalid_rtree = 1;
	    }
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    if (invalid_rtree)
	return 0;
    return 1;
}

static void
fnct_CheckSpatialIndex (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ CheckSpatialIndex()
/ CheckSpatialIndex(table, column)
/
/ checks a SpatialIndex for consistency, returning:
/ 1 - the R*Tree is fully consistent
/ 0 - the R*Tree is inconsistent
/ NULL on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    int status;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (argc == 0)
      {
	  /* no arguments: we must check any defined R*Tree */
	  status = check_any_spatial_index (sqlite);
	  if (status < 0)
	      sqlite3_result_null (context);
	  else if (status > 0)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }

    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CheckSpatialIndex() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_null (context);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CheckSpatialIndex() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_null (context);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    status = check_spatial_index (sqlite, table, column);
    if (status < 0)
	sqlite3_result_null (context);
    else if (status > 0)
	sqlite3_result_int (context, 1);
    else
	sqlite3_result_int (context, 0);
}

static int
recover_spatial_index (sqlite3 * sqlite, const unsigned char *table,
		       const unsigned char *geom)
{
/* attempting to rebuild an R*Tree */
    char sql[8192];
    char sql2[2048];
    char *errMsg = NULL;
    int ret;
    char xtable[1024];
    char xgeom[1024];
    char idx_name[2048];
    int is_defined = 0;
    sqlite3_stmt *stmt;

/* checking if the R*Tree Spatial Index is defined */
    strcpy (xtable, (const char *) table);
    clean_sql_string (xtable);
    strcpy (xgeom, (const char *) geom);
    clean_sql_string (xgeom);
    strcpy (sql, "SELECT Count(*) FROM geometry_columns ");
    sprintf (sql2, "WHERE f_table_name LIKE '%s' ", xtable);
    strcat (sql, sql2);
    sprintf (sql2, "AND f_geometry_column LIKE '%s' ", xgeom);
    strcat (sql, sql2);
    strcat (sql, "AND spatial_index_enabled = 1");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "RecoverSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	      is_defined = sqlite3_column_int (stmt, 0);
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    if (!is_defined)
	return -1;

/* erasing the R*Tree table */
    sprintf (idx_name, "idx_%s_%s", table, geom);
    double_quoted_sql (idx_name);
    sprintf (sql, "DELETE FROM %s", idx_name);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
/* populating the R*Tree table from scratch */
    buildSpatialIndex (sqlite, table, (const char *) geom);
    strcpy (sql, "SpatialIndex: successfully recovered");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) geom, sql);
    return 1;
  error:
    fprintf (stderr, "RecoverSpatialIndex() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    return 0;
}

static int
recover_any_spatial_index (sqlite3 * sqlite, int no_check)
{
/* attempting to rebuild any defined R*Tree */
    const unsigned char *table;
    const unsigned char *column;
    int status;
    char sql[1024];
    int ret;
    int to_be_fixed;
    sqlite3_stmt *stmt;

/* retrieving any defined R*Tree */
    strcpy (sql,
	    "SELECT f_table_name, f_geometry_column FROM geometry_columns ");
    strcat (sql, "WHERE spatial_index_enabled = 1");
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "RecoverSpatialIndex SQL error: %s\n",
		   sqlite3_errmsg (sqlite));
	  return -1;
      }
    while (1)
      {
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;
	  if (ret == SQLITE_ROW)
	    {
		/* checking a single R*Tree */
		table = sqlite3_column_text (stmt, 0);
		column = sqlite3_column_text (stmt, 1);
		to_be_fixed = 1;
		if (!no_check)
		  {
		      status = check_spatial_index (sqlite, table, column);
		      if (status < 0)
			{
			    /* some unexpected error occurred */
			    goto fatal_error;
			}
		      else if (status > 0)
			{
			    /* the Spatial Index is already valid */
			    to_be_fixed = 0;
			}
		  }
		if (to_be_fixed)
		  {
		      /* rebuilding the Spatial Index */
		      status = recover_spatial_index (sqlite, table, column);
		      if (status < 0)
			{
			    /* some unexpected error occurred */
			    goto fatal_error;
			}
		      else if (status == 0)
			  goto error;
		  }
	    }
	  else
	    {
		printf ("sqlite3_step() error: %s\n", sqlite3_errmsg (sqlite));
		sqlite3_finalize (stmt);
		return -1;
	    }
      }
    sqlite3_finalize (stmt);
    return 1;
  error:
    sqlite3_finalize (stmt);
    return 0;
  fatal_error:
    sqlite3_finalize (stmt);
    return -1;
}

static void
fnct_RecoverSpatialIndex (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ RecoverSpatialIndex()
/ RecoverSpatialIndex(no_check)
/ RecoverSpatialIndex(table, column)
/ RecoverSpatialIndex(table, column, no_check)
/
/ attempts to rebuild a SpatialIndex, returning:
/ 1 - on success
/ 0 - on failure
/ NULL if any syntax error is detected
*/
    const unsigned char *table;
    const unsigned char *column;
    int no_check = 0;
    int status;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (argc <= 1)
      {
	  /* no arguments: we must rebuild any defined R*Tree */
	  if (argc == 1)
	    {
		if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
		    no_check = sqlite3_value_int (argv[0]);
		else
		  {
		      fprintf (stderr,
			       "RecoverSpatialIndex() error: argument 1 [no_check] is not of the Integer type\n");
		      sqlite3_result_null (context);
		      return;
		  }
	    }
	  status = recover_any_spatial_index (sqlite, no_check);
	  if (status < 0)
	      sqlite3_result_null (context);
	  else if (status > 0)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }

    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverSpatialIndex() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_null (context);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RecoverSpatialIndex() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_null (context);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    if (argc == 3)
      {
	  if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	      no_check = sqlite3_value_int (argv[2]);
	  else
	    {
		fprintf (stderr,
			 "RecoverSpatialIndex() error: argument 2 [no_check] is not of the Integer type\n");
		sqlite3_result_null (context);
		return;
	    }
      }
    if (!no_check)
      {
	  /* checking the current SpatialIndex validity */
	  status = check_spatial_index (sqlite, table, column);
	  if (status < 0)
	    {
		/* some unexpected error occurred */
		sqlite3_result_null (context);
		return;
	    }
	  else if (status > 0)
	    {
		/* the Spatial Index is already valid */
		sqlite3_result_int (context, 1);
		return;
	    }
      }
/* rebuilding the Spatial Index */
    status = recover_spatial_index (sqlite, table, column);
    if (status < 0)
	sqlite3_result_null (context);
    else if (status > 0)
	sqlite3_result_int (context, 1);
    else
	sqlite3_result_int (context, 0);
}

static void
fnct_CreateSpatialIndex (sqlite3_context * context, int argc,
			 sqlite3_value ** argv)
{
/* SQL function:
/ CreateSpatialIndex(table, column )
/
/ creates a SpatialIndex based on Column and Table
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CreateSpatialIndex() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CreateSpatialIndex() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    strcpy (sql,
	    "UPDATE geometry_columns SET spatial_index_enabled = 1 WHERE f_table_name LIKE '");
    strcat (sql, sqltable);
    strcat (sql, "' AND f_geometry_column LIKE '");
    strcat (sql, sqlcolumn);
    strcat (sql, "' AND spatial_index_enabled = 0");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    if (sqlite3_changes (sqlite) == 0)
      {
	  fprintf (stderr,
		   "CreateSpatialIndex() error: either \"%s\".\"%s\" isn't a Geometry column or a SpatialIndex is already defined\n",
		   table, column);
	  sqlite3_result_int (context, 0);
	  return;
      }
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "R*Tree Spatial Index successfully created");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "CreateSpatialIndex() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_CreateMbrCache (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CreateMbrCache(table, column )
/
/ creates an MBR Cache based on Column and Table
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CreateMbrCache() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "CreateMbrCache() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    strcpy (sql,
	    "UPDATE geometry_columns SET spatial_index_enabled = 2 WHERE f_table_name LIKE '");
    strcat (sql, sqltable);
    strcat (sql, "' AND f_geometry_column LIKE '");
    strcat (sql, sqlcolumn);
    strcat (sql, "' AND spatial_index_enabled = 0");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    if (sqlite3_changes (sqlite) == 0)
      {
	  fprintf (stderr,
		   "CreateMbrCache() error: either \"%s\".\"%s\" isn't a Geometry column or a SpatialIndex is already defined\n",
		   table, column);
	  sqlite3_result_int (context, 0);
	  return;
      }
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "MbrCache successfully created");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "CreateMbrCache() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_DisableSpatialIndex (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ DisableSpatialIndex(table, column )
/
/ disables a SpatialIndex based on Column and Table
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DisableSpatialIndex() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "DisableSpatialIndex() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    strcpy (sql,
	    "UPDATE geometry_columns SET spatial_index_enabled = 0 WHERE f_table_name LIKE '");
    strcat (sql, sqltable);
    strcat (sql, "' AND f_geometry_column LIKE '");
    strcat (sql, sqlcolumn);
    strcat (sql, "' AND spatial_index_enabled <> 0");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &errMsg);
    if (ret != SQLITE_OK)
	goto error;
    if (sqlite3_changes (sqlite) == 0)
      {
	  fprintf (stderr,
		   "DisableSpatialIndex() error: either \"%s\".\"%s\" isn't a Geometry column or no SpatialIndex is defined\n",
		   table, column);
	  sqlite3_result_int (context, 0);
	  return;
      }
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "SpatialIndex successfully disabled");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "DisableSpatialIndex() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_RebuildGeometryTriggers (sqlite3_context * context, int argc,
			      sqlite3_value ** argv)
{
/* SQL function:
/ RebuildGeometryTriggers(table, column )
/
/ rebuilds Geometry Triggers (constraints)  based on Column and Table
/ returns 1 on success
/ 0 on failure
*/
    const unsigned char *table;
    const unsigned char *column;
    char sql[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    char sqltable[1024];
    char sqlcolumn[1024];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RebuildGeometryTriggers() error: argument 1 [table_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    table = sqlite3_value_text (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
      {
	  fprintf (stderr,
		   "RebuildGeometryTriggers() error: argument 2 [column_name] is not of the String type\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    column = sqlite3_value_text (argv[1]);
    strcpy (sqltable, (char *) table);
    clean_sql_string (sqltable);
    strcpy (sqlcolumn, (char *) column);
    clean_sql_string (sqlcolumn);
    strcpy (sql,
	    "SELECT f_table_name FROM geometry_columns WHERE f_table_name LIKE '");
    strcat (sql, sqltable);
    strcat (sql, "' AND f_geometry_column LIKE '");
    strcat (sql, sqlcolumn);
    strcat (sql, "'");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	goto error;
    sqlite3_free_table (results);
    if (rows <= 0)
      {
	  fprintf (stderr,
		   "RebuildGeometryTriggers() error: \"%s\".\"%s\" isn't a Geometry column\n",
		   table, column);
	  sqlite3_result_int (context, 0);
	  return;
      }
    updateGeometryTriggers (sqlite, table, column);
    sqlite3_result_int (context, 1);
    strcpy (sql, "Geometry Triggers successfully rebuilt");
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    fprintf (stderr, "RebuildGeometryTriggers() error: \"%s\"\n", errMsg);
    sqlite3_free (errMsg);
    sqlite3_result_int (context, 0);
    return;
}


static int
check_topo_table (sqlite3 * sqlite, const char *table, int is_view)
{
/* checking if some Topology-related table/view already exists */
    int exists = 0;
    char sql[2048];
    char sqltable[1024];
    char *errMsg = NULL;
    int ret;
    char **results;
    int rows;
    int columns;
    int i;
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT name FROM sqlite_master WHERE type = '%s' AND name LIKE '%s'",
	     (!is_view) ? "table" : "view", table);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
	return 0;
    for (i = 1; i <= rows; i++)
	exists = 1;
    sqlite3_free_table (results);
    return exists;
}

static int
create_topo_nodes (sqlite3 * sqlite, const char *table, int srid, int dims)
{
/* creating the topo_nodes table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "node_id INTEGER PRIMARY KEY AUTOINCREMENT,\n");
    strcat (sql, "gml_id TEXT)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT AddGeometryColumn('%s', 'Geometry', %d, 'POINT', '%s')",
	     sqltable, srid, (dims == GAIA_XY_Z) ? "XYZ" : "XY");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    sprintf (sql, "SELECT CreateSpatialIndex('%s', 'Geometry')", sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CreateSpatialIndex '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_gml", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (gml_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('gml_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_edges (sqlite3 * sqlite, const char *table, int srid, int dims)
{
/* creating the topo_edges table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "edge_id INTEGER PRIMARY KEY AUTOINCREMENT,\n");
    strcat (sql, "node_from_href TEXT,\n");
    strcat (sql, "node_to_href TEXT,\n");
    strcat (sql, "gml_id TEXT)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT AddGeometryColumn('%s', 'Geometry', %d, 'LINESTRING', '%s')",
	     sqltable, srid, (dims == GAIA_XY_Z) ? "XYZ" : "XY");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    sprintf (sql, "SELECT CreateSpatialIndex('%s', 'Geometry')", sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CreateSpatialIndex '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_gml_id", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (gml_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('gml_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_from", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (node_from_href)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('node_from_href') error: %s\n",
		   sqltable, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_to", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (node_to_href)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('node_to_href') error: %s\n",
		   sqltable, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_faces (sqlite3 * sqlite, const char *table, int srid, int dims)
{
/* creating the topo_faces table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "face_id INTEGER PRIMARY KEY AUTOINCREMENT,\n");
    strcat (sql, "gml_id TEXT)");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT AddGeometryColumn('%s', 'Geometry', %d, 'POLYGON', '%s')",
	     sqltable, srid, (dims == GAIA_XY_Z) ? "XYZ" : "XY");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    sprintf (sql, "SELECT CreateSpatialIndex('%s', 'Geometry')", sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CreateSpatialIndex '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_gml", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (gml_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('gml_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_faces_edges (sqlite3 * sqlite, const char *table,
			 const char *table2)
{
/* creating the topo_faces_edges table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "face_id INTEGER NOT NULL,\n");
    strcat (sql, "sub INTEGER NOT NULL,\n");
    strcat (sql, "gml_id TEXT,\n");
    strcat (sql, "orientation TEXT,\n");
    strcat (sql, "CONSTRAINT pk_faces_edges PRIMARY KEY ");
    strcat (sql, "(face_id, sub),\n");
    strcat (sql, "CONSTRAINT fk_faces_edges FOREIGN KEY ");
    strcat (sql, "(face_id) REFERENCES ");
    strcpy (sql2, table2);
    double_quoted_sql (sql2);
    strcat (sql, sql2);
    strcat (sql, " (face_id))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_edge", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (gml_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('gml_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_curves (sqlite3 * sqlite, const char *table, const char *table2,
		    int srid, int dims)
{
/* creating the topo_curves table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "curve_id INTEGER NOT NULL,\n");
    strcat (sql, "edge_id INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_curves PRIMARY KEY ");
    strcat (sql, "(curve_id, edge_id),\n");
    strcat (sql, "CONSTRAINT fk_curves FOREIGN KEY ");
    strcat (sql, "(edge_id) REFERENCES ");
    strcpy (sql2, table2);
    double_quoted_sql (sql2);
    strcat (sql, sql2);
    strcat (sql, " (edge_id))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT AddGeometryColumn('%s', 'Geometry', %d, 'MULTILINESTRING', '%s')",
	     sqltable, srid, (dims == GAIA_XY_Z) ? "XYZ" : "XY");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    sprintf (sql, "SELECT CreateSpatialIndex('%s', 'Geometry')", sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CreateSpatialIndex '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_edge", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (edge_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('edge_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_surfaces (sqlite3 * sqlite, const char *table, const char *table2,
		      int srid, int dims)
{
/* creating the topo_surfaces table */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE TABLE %s (\n", sqltable);
    strcat (sql, "surface_id INTEGER NOT NULL,\n");
    strcat (sql, "face_id INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT pk_surfaces PRIMARY KEY ");
    strcat (sql, "(surface_id, face_id),\n");
    strcat (sql, "CONSTRAINT fk_surfaces FOREIGN KEY ");
    strcat (sql, "(face_id) REFERENCES ");
    strcpy (sql2, table2);
    double_quoted_sql (sql2);
    strcat (sql, sql2);
    strcat (sql, " (face_id))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE '%s' error: %s\n", table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    clean_sql_string (sqltable);
    sprintf (sql,
	     "SELECT AddGeometryColumn('%s', 'Geometry', %d, 'MULTIPOLYGON', '%s')",
	     sqltable, srid, (dims == GAIA_XY_Z) ? "XYZ" : "XY");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "AddGeometryColumn '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    sprintf (sql, "SELECT CreateSpatialIndex('%s', 'Geometry')", sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CreateSpatialIndex '%s'.'Geometry' error: %s\n",
		   table, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    strcpy (sqltable, table);
    double_quoted_sql (sqltable);
    sprintf (sql2, "idx_%s_face", sqltable);
    double_quoted_sql (sql2);
    sprintf (sql, "CREATE INDEX %s ON %s (face_id)", sql2, sqltable);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "Create Index '%s'('face_id') error: %s\n", sqltable,
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_check_node_ids (sqlite3 * sqlite, const char *view,
		       const char *table_nodes)
{
/* creating the check node ids VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT gml_id AS gml_id, Count(node_id) AS count\n");
    strcpy (sqltable, table_nodes);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "GROUP BY gml_id\n");
    strcat (sql, "HAVING count > 1\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_check_node_geoms (sqlite3 * sqlite, const char *view,
			 const char *table_nodes)
{
/* creating the check node geoms VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT n1.node_id AS node1_id, n1.gml_id AS node1_gml_id, ");
    strcat (sql, "n2.node_id AS node2_id, n2.gml_id AS node2_gml_id\n");
    strcpy (sqltable, table_nodes);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s AS n1\n", sqltable);
    strcat (sql, sql2);
    sprintf (sql2, "JOIN %s AS n2 ON (\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "  n1.node_id <> n2.node_id AND\n");
    strcat (sql, "  ST_Equals(n1.Geometry, n2.Geometry) = 1 AND\n");
    strcat (sql, "  n2.node_id IN (\n");
    strcat (sql, "	SELECT ROWID FROM SpatialIndex\n");
    strcpy (sqltable, table_nodes);
    clean_sql_string (sqltable);
    sprintf (sql2, "	WHERE f_table_name = '%s' AND\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "	  search_frame = n1.Geometry))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_check_edge_ids (sqlite3 * sqlite, const char *view,
		       const char *table_edges)
{
/* creating the check edge ids VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT gml_id AS gml_id, Count(edge_id) AS count\n");
    strcpy (sqltable, table_edges);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "GROUP BY gml_id\n");
    strcat (sql, "HAVING count > 1\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_check_edge_geoms (sqlite3 * sqlite, const char *view,
			 const char *table_edges)
{
/* creating the check edge geoms VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT e1.edge_id AS edge1_id, e1.gml_id AS edge1_gml_id, ");
    strcat (sql, "e2.edge_id AS edge2_id, e2.gml_id AS edge2_gml_id\n");
    strcpy (sqltable, table_edges);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s AS e1\n", sqltable);
    strcat (sql, sql2);
    sprintf (sql2, "JOIN %s AS e2 ON (\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "  e1.edge_id <> e2.edge_id AND\n");
    strcat (sql, "NOT (e1.node_from_href = e2.node_from_href ");
    strcat (sql, "AND e1.node_to_href = e2.node_to_href) AND\n");
    strcat (sql, "  ST_Crosses(e1.Geometry, e2.Geometry) = 1 AND\n");
    strcat (sql, "  e2.edge_id IN (\n");
    strcat (sql, "	SELECT ROWID FROM SpatialIndex\n");
    strcpy (sqltable, table_edges);
    clean_sql_string (sqltable);
    sprintf (sql2, "	WHERE f_table_name = '%s' AND\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "	  search_frame = e1.Geometry))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_check_face_ids (sqlite3 * sqlite, const char *view,
		       const char *table_faces)
{
/* creating the check face ids VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT gml_id AS gml_id, Count(face_id) AS count\n");
    strcpy (sqltable, table_faces);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "GROUP BY gml_id\n");
    strcat (sql, "HAVING count > 1\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_faces_resolved (sqlite3 * sqlite, const char *view, const char *faces,
		       const char *faces_edges, const char *edges)
{
/* creating the Faces Resolved VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT f.face_id AS face_id, ");
    strcat (sql, "ST_Polygonize(e.Geometry) AS Geometry\n");
    strcpy (sqltable, faces);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s AS f\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "LEFT JOIN ");
    strcpy (sqltable, faces_edges);
    double_quoted_sql (sqltable);
    strcat (sql, sqltable);
    double_quoted_sql (sqltable);
    strcat (sql, " AS fe ON (fe.face_id = f.face_id)\n");
    strcat (sql, "LEFT JOIN ");
    strcpy (sqltable, edges);
    double_quoted_sql (sqltable);
    strcat (sql, sqltable);
    double_quoted_sql (sqltable);
    strcat (sql, " AS e ON (e.gml_id = fe.gml_id)\n");
    strcat (sql, "GROUP BY f.face_id\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_curves_resolved (sqlite3 * sqlite, const char *view, const char *curves,
			char *edges)
{
/* creating the Curves Resolved VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT c.curve_id AS curve_id, ");
    strcat (sql, "CastToMultiLinestring(Collect(e.Geometry)) AS Geometry\n");
    strcpy (sqltable, curves);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s AS c\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "LEFT JOIN ");
    strcpy (sqltable, edges);
    double_quoted_sql (sqltable);
    strcat (sql, sqltable);
    double_quoted_sql (sqltable);
    strcat (sql, " AS e ON (e.edge_id = c.edge_id)\n");
    strcat (sql, "GROUP BY c.curve_id\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_surfaces_resolved (sqlite3 * sqlite, const char *view,
			  const char *surfaces, const char *faces)
{
/* creating the Surfaces Resolved VIEW */
    char sql[2048];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;
    strcpy (sqltable, view);
    double_quoted_sql (sqltable);
    sprintf (sql, "CREATE VIEW %s AS\n", sqltable);
    strcat (sql, "SELECT s.surface_id AS surface_id, ");
    strcat (sql, "CastToMultipolygon(Collect(f.Geometry)) AS Geometry\n");
    strcpy (sqltable, surfaces);
    double_quoted_sql (sqltable);
    sprintf (sql2, "FROM %s AS s\n", sqltable);
    strcat (sql, sql2);
    strcat (sql, "LEFT JOIN ");
    strcpy (sqltable, faces);
    double_quoted_sql (sqltable);
    strcat (sql, sqltable);
    double_quoted_sql (sqltable);
    strcat (sql, " AS f ON (f.face_id = s.face_id)\n");
    strcat (sql, "GROUP BY s.surface_id\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE VIEW '%s' error: %s\n", view, err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
create_topo_master (sqlite3 * sqlite)
{
/* creating the topo_master table */
    char sql[4196];
    int ret;
    char *err_msg = NULL;

/* creating the table */
    strcpy (sql, "CREATE TABLE topology_master (\n");
    strcat (sql, "nodes TEXT NOT NULL,\n");
    strcat (sql, "edges TEXT NOT NULL,\n");
    strcat (sql, "faces TEXT NOT NULL,\n");
    strcat (sql, "faces_edges TEXT NOT NULL,\n");
    strcat (sql, "curves TEXT NOT NULL,\n");
    strcat (sql, "surfaces TEXT NOT NULL,\n");
    strcat (sql, "check_node_ids TEXT NOT NULL,\n");
    strcat (sql, "check_node_geoms TEXT NOT NULL,\n");
    strcat (sql, "check_edge_ids TEXT NOT NULL,\n");
    strcat (sql, "check_edge_geoms TEXT NOT NULL,\n");
    strcat (sql, "check_face_ids TEXT NOT NULL,\n");
    strcat (sql, "faces_resolved TEXT NOT NULL,\n");
    strcat (sql, "curves_resolved TEXT NOT NULL,\n");
    strcat (sql, "surfaces_resolved TEXT NOT NULL,\n");
    strcat (sql, "coord_dimension TEXT NOT NULL,\n");
    strcat (sql, "srid INTEGER NOT NULL,\n");
    strcat (sql, "CONSTRAINT fk_topo_master FOREIGN KEY \n");
    strcat (sql, "(srid) REFERENCES spatial_ref_sys (srid))\n");
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "CREATE TABLE 'topology_master' error: %s\n",
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static int
update_topo_master (sqlite3 * sqlite, const char *nodes, const char *edges,
		    const char *faces, const char *faces_edges,
		    const char *curves, const char *surfaces,
		    const char *check_nodes, const char *check_node_geoms,
		    const char *check_edges, const char *check_edge_geoms,
		    const char *check_faces, const char *faces_res,
		    const char *curves_res, const char *surfaces_res, int srid,
		    int dims)
{
/* updating the topo_master table */
    char sql[4196];
    char sql2[2048];
    char sqltable[1024];
    int ret;
    char *err_msg = NULL;

/* inserting Topology data into MASTER */
    strcpy (sql, "INSERT INTO topology_master ");
    strcat (sql, "(nodes, edges, faces, faces_edges, ");
    strcat (sql, "curves, surfaces, check_node_ids, ");
    strcat (sql, "check_node_geoms, check_edge_ids, ");
    strcat (sql, "check_edge_geoms, check_face_ids, ");
    strcat (sql, "faces_resolved, curves_resolved, ");
    strcat (sql, "surfaces_resolved, ");
    strcat (sql, "coord_dimension, srid) ");
    strcat (sql, "VALUES (");
    strcpy (sqltable, nodes);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, edges);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, faces);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, faces_edges);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, curves);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, surfaces);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, check_nodes);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, check_node_geoms);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, check_edges);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, check_edge_geoms);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, check_faces);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, faces_res);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, curves_res);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    strcpy (sqltable, surfaces_res);
    clean_sql_string (sqltable);
    sprintf (sql2, "'%s', ", sqltable);
    strcat (sql, sql2);
    sprintf (sql2, "'%s', %d)", (dims == GAIA_XY_Z) ? "XYZ" : "XY", srid);
    strcat (sql, sql2);
    ret = sqlite3_exec (sqlite, sql, NULL, NULL, &err_msg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "INSERT INTO 'topology_master' error: %s\n",
		   err_msg);
	  sqlite3_free (err_msg);
	  return 0;
      }
    return 1;
}

static void
fnct_CreateTopologyTables (sqlite3_context * context, int argc,
			   sqlite3_value ** argv)
{
/* SQL function:
/ CreateTopologyTables(srid, coord_dims)
/  or
/ CreateTopologyTables(prefix, srid, coord_dims)
/
/ creates any Topology related table 
/ returns 1 on success
/ 0 on failure
*/
    const char *prefix = "topo_";
    const unsigned char *txt_dims;
    int srid = -1;
    int dimension;
    int dims = -1;
    char table_curves[1024];
    char table_surfaces[1024];
    char table_nodes[1024];
    char table_edges[1024];
    char table_faces[1024];
    char table_faces_edges[1024];
    char view_check_node_ids[1024];
    char view_check_node_geoms[1024];
    char view_check_edge_ids[1024];
    char view_check_edge_geoms[1024];
    char view_check_face_ids[1024];
    char view_faces_resolved[1024];
    char view_curves_resolved[1024];
    char view_surfaces_resolved[1024];
    const char *tables[16];
    int views[16];
    int *p_view;
    const char **p_tbl;
    int ok_table;
    int create_master = 1;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (argc == 3)
      {
	  if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
	    {
		fprintf (stderr,
			 "CreateTopologyTables() error: argument 1 [table_prefix] is not of the String type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  prefix = (char *) sqlite3_value_text (argv[0]);
	  if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
	    {
		fprintf (stderr,
			 "CreateTopologyTables() error: argument 2 [SRID] is not of the Integer type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  srid = sqlite3_value_int (argv[1]);
	  if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	    {
		dimension = sqlite3_value_int (argv[2]);
		if (dimension == 2)
		    dims = GAIA_XY;
		if (dimension == 3)
		    dims = GAIA_XY_Z;
	    }
	  else if (sqlite3_value_type (argv[2]) == SQLITE_TEXT)
	    {
		txt_dims = sqlite3_value_text (argv[2]);
		if (strcasecmp ((char *) txt_dims, "XY") == 0)
		    dims = GAIA_XY;
		if (strcasecmp ((char *) txt_dims, "XYZ") == 0)
		    dims = GAIA_XY_Z;
	    }
	  else
	    {
		fprintf (stderr,
			 "CreateTopologyTables() error: argument 3 [dimension] is not of the Integer or Text type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
      }
    else
      {
	  if (sqlite3_value_type (argv[0]) != SQLITE_INTEGER)
	    {
		fprintf (stderr,
			 "CreateTopologyTables() error: argument 1 [SRID] is not of the Integer type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  srid = sqlite3_value_int (argv[0]);
	  if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	    {
		dimension = sqlite3_value_int (argv[1]);
		if (dimension == 2)
		    dims = GAIA_XY;
		if (dimension == 3)
		    dims = GAIA_XY_Z;
	    }
	  else if (sqlite3_value_type (argv[1]) == SQLITE_TEXT)
	    {
		txt_dims = sqlite3_value_text (argv[1]);
		if (strcasecmp ((char *) txt_dims, "XY") == 0)
		    dims = GAIA_XY;
		if (strcasecmp ((char *) txt_dims, "XYZ") == 0)
		    dims = GAIA_XY_Z;
	    }
	  else
	    {
		fprintf (stderr,
			 "CreateTopologyTables() error: argument 2 [dimension] is not of the Integer or Text type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
      }
    if (dims == GAIA_XY || dims == GAIA_XY_Z)
	;
    else
      {
	  fprintf (stderr,
		   "CreateTopologyTables() error: [dimension] ILLEGAL VALUE\n");
	  sqlite3_result_int (context, 0);
	  return;
      }
    if (srid <= 0)
      {
	  fprintf (stderr,
		   "CreateTopologyTables() error: [SRID] ILLEGAL VALUE\n");
	  sqlite3_result_int (context, 0);
	  return;
      }

/* checking Topology tables */
    tables[0] = "topology_master";
    views[0] = 0;
    sprintf (table_curves, "%scurves", prefix);
    tables[1] = table_curves;
    views[1] = 0;
    sprintf (table_surfaces, "%ssurfaces", prefix);
    tables[2] = table_surfaces;
    views[2] = 0;
    sprintf (table_nodes, "%snodes", prefix);
    tables[3] = table_nodes;
    views[3] = 0;
    sprintf (table_edges, "%sedges", prefix);
    tables[4] = table_edges;
    views[4] = 0;
    sprintf (table_faces, "%sfaces", prefix);
    tables[5] = table_faces;
    views[5] = 0;
    sprintf (table_faces_edges, "%sfaces_edges", prefix);
    tables[6] = table_faces_edges;
    views[6] = 0;
    sprintf (view_check_node_ids, "%snodes_check_dupl_ids", prefix);
    tables[7] = view_check_node_ids;
    views[7] = 1;
    sprintf (view_check_node_geoms, "%snodes_check_dupl_geoms", prefix);
    tables[8] = view_check_node_geoms;
    views[8] = 1;
    sprintf (view_check_edge_ids, "%sedges_check_dupl_ids", prefix);
    tables[9] = view_check_edge_ids;
    views[9] = 1;
    sprintf (view_check_edge_geoms, "%sedges_check_dupl_geoms", prefix);
    tables[10] = view_check_edge_geoms;
    views[10] = 1;
    sprintf (view_check_face_ids, "%sfaces_check_dupl_ids", prefix);
    tables[11] = view_check_face_ids;
    views[11] = 1;
    sprintf (view_faces_resolved, "%sfaces_resolved", prefix);
    tables[12] = view_faces_resolved;
    views[12] = 1;
    sprintf (view_curves_resolved, "%scurves_resolved", prefix);
    tables[13] = view_curves_resolved;
    views[13] = 1;
    sprintf (view_surfaces_resolved, "%ssurfaces_resolved", prefix);
    tables[14] = view_surfaces_resolved;
    views[14] = 1;
    tables[15] = NULL;
    p_view = views;
    p_tbl = tables;
    while (*p_tbl != NULL)
      {
	  ok_table = check_topo_table (sqlite, *p_tbl, *p_view);
	  if (ok_table)
	    {
		if (strcmp (*p_tbl, "topology_master") == 0)
		    create_master = 0;
		else
		  {
		      fprintf (stderr,
			       "CreateTopologyTables() error: table '%s' already exists\n",
			       *p_tbl);
		      sqlite3_result_int (context, 0);
		      return;
		  }
	    }
	  p_tbl++;
	  p_view++;
      }

/* creating Topology tables */
    if (create_master)
      {
	  if (!create_topo_master (sqlite))
	      goto error;
      }
    if (!create_topo_nodes (sqlite, table_nodes, srid, dims))
	goto error;
    if (!create_topo_edges (sqlite, table_edges, srid, dims))
	goto error;
    if (!create_topo_faces (sqlite, table_faces, srid, dims))
	goto error;
    if (!create_topo_faces_edges (sqlite, table_faces_edges, table_faces))
	goto error;
    if (!create_topo_curves (sqlite, table_curves, table_edges, srid, dims))
	goto error;
    if (!create_topo_surfaces (sqlite, table_surfaces, table_faces, srid, dims))
	goto error;
    if (!create_check_node_ids (sqlite, view_check_node_ids, table_nodes))
	goto error;
    if (!create_check_node_geoms (sqlite, view_check_node_geoms, table_nodes))
	goto error;
    if (!create_check_edge_ids (sqlite, view_check_edge_ids, table_edges))
	goto error;
    if (!create_check_edge_geoms (sqlite, view_check_edge_geoms, table_edges))
	goto error;
    if (!create_check_face_ids (sqlite, view_check_face_ids, table_faces))
	goto error;
    if (!create_faces_resolved
	(sqlite, view_faces_resolved, table_faces, table_faces_edges,
	 table_edges))
	goto error;
    if (!create_curves_resolved
	(sqlite, view_curves_resolved, table_curves, table_edges))
	goto error;
    if (!create_surfaces_resolved
	(sqlite, view_surfaces_resolved, table_surfaces, table_faces))
	goto error;
    if (!update_topo_master
	(sqlite, table_nodes, table_edges, table_faces, table_faces_edges,
	 table_curves, table_surfaces, view_check_node_ids,
	 view_check_node_geoms, view_check_edge_ids, view_check_edge_geoms,
	 view_check_face_ids, view_faces_resolved, view_curves_resolved,
	 view_surfaces_resolved, srid, dims))
	goto error;
    updateSpatiaLiteHistory (sqlite, "*** TOPOLOGY ***", NULL,
			     "Topology tables successfully created");
    sqlite3_result_int (context, 1);
    return;

  error:
    sqlite3_result_int (context, 0);
    return;
}

static void
fnct_UpdateLayerStatistics (sqlite3_context * context, int argc,
			    sqlite3_value ** argv)
{
/* SQL function:
/ UpdateLayerStatistics(table, column )
/
/ Updates LAYER_STATISTICS [based on Column and Table]
/ returns 1 on success
/ 0 on failure
*/
    const char *sql;
    const unsigned char *table = NULL;
    const unsigned char *column = NULL;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (argc >= 1)
      {
	  if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
	    {
		fprintf (stderr,
			 "UpdateLayerStatistics() error: argument 1 [table_name] is not of the String type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  table = sqlite3_value_text (argv[0]);
      }
    if (argc >= 2)
      {
	  if (sqlite3_value_type (argv[1]) != SQLITE_TEXT)
	    {
		fprintf (stderr,
			 "UpdateLayerStatistics() error: argument 2 [column_name] is not of the String type\n");
		sqlite3_result_int (context, 0);
		return;
	    }
	  column = sqlite3_value_text (argv[1]);
      }
    if (!update_layer_statistics (sqlite, table, column))
	goto error;
    sqlite3_result_int (context, 1);
    sql = "UpdateLayerStatistics";
    if (table == NULL)
	table = "ALL-TABLES";
    if (column == NULL)
	column = "ALL-GEOMETRY-COLUMNS";
    updateSpatiaLiteHistory (sqlite, (const char *) table,
			     (const char *) column, sql);
    return;
  error:
    sqlite3_result_int (context, 0);
    return;
}

static gaiaPointPtr
simplePoint (gaiaGeomCollPtr geo)
{
/* helper function
/ if this GEOMETRY contains only one POINT, and no other elementary geometry
/ the POINT address will be returned
/ otherwise NULL will be returned
*/
    int cnt = 0;
    gaiaPointPtr point;
    gaiaPointPtr this_point = NULL;
    if (!geo)
	return NULL;
    if (geo->FirstLinestring || geo->FirstPolygon)
	return NULL;
    point = geo->FirstPoint;
    while (point)
      {
	  /* counting how many POINTs are there */
	  cnt++;
	  this_point = point;
	  point = point->Next;
      }
    if (cnt == 1 && this_point)
	return this_point;
    return NULL;
}

static gaiaLinestringPtr
simpleLinestring (gaiaGeomCollPtr geo)
{
/* helper function
/ if this GEOMETRY contains only one LINESTRING, and no other elementary geometry
/ the LINESTRING address will be returned
/ otherwise NULL will be returned
*/
    int cnt = 0;
    gaiaLinestringPtr line;
    gaiaLinestringPtr this_line = NULL;
    if (!geo)
	return NULL;
    if (geo->FirstPoint || geo->FirstPolygon)
	return NULL;
    line = geo->FirstLinestring;
    while (line)
      {
	  /* counting how many LINESTRINGs are there */
	  cnt++;
	  this_line = line;
	  line = line->Next;
      }
    if (cnt == 1 && this_line)
	return this_line;
    return NULL;
}

static gaiaPolygonPtr
simplePolygon (gaiaGeomCollPtr geo)
{
/* helper function
/ if this GEOMETRY contains only one POLYGON, and no other elementary geometry
/ the POLYGON address will be returned
/ otherwise NULL will be returned
*/
    int cnt = 0;
    gaiaPolygonPtr polyg;
    gaiaPolygonPtr this_polyg = NULL;
    if (!geo)
	return NULL;
    if (geo->FirstPoint || geo->FirstLinestring)
	return NULL;
    polyg = geo->FirstPolygon;
    while (polyg)
      {
	  /* counting how many POLYGONs are there */
	  cnt++;
	  this_polyg = polyg;
	  polyg = polyg->Next;
      }
    if (cnt == 1 && this_polyg)
	return this_polyg;
    return NULL;
}

static void
fnct_AsText (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsText(BLOB encoded geometry)
/
/ returns the corresponding WKT encoded value
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    gaiaOutBufferInitialize (&out_buf);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaOutWkt (&out_buf, geo);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsWkt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsWkt(BLOB encoded geometry [, Integer precision])
/
/ returns the corresponding WKT encoded value
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int precision = 15;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (argc == 2)
      {
	  if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	      precision = sqlite3_value_int (argv[1]);
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    gaiaOutBufferInitialize (&out_buf);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaOutWktStrict (&out_buf, geo, precision);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

/*
/
/ AsSvg(geometry,[relative], [precision]) implementation
/
////////////////////////////////////////////////////////////
/
/ Author: Klaus Foerster klaus.foerster@svg.cc
/ version 0.9. 2008 September 21
 /
 */

static void
fnct_AsSvg (sqlite3_context * context, int argc, sqlite3_value ** argv,
	    int relative, int precision)
{
/* SQL function:
   AsSvg(BLOB encoded geometry, [int relative], [int precision])
   returns the corresponding SVG encoded value or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  /* make sure relative is 0 or 1 */
	  if (relative > 0)
	      relative = 1;
	  else
	      relative = 0;
	  /* make sure precision is between 0 and 15 - default to 6 if absent */
	  if (precision > GAIA_SVG_DEFAULT_MAX_PRECISION)
	      precision = GAIA_SVG_DEFAULT_MAX_PRECISION;
	  if (precision < 0)
	      precision = 0;
	  /* produce SVG-notation - actual work is done in gaiageo/gg_wkt.c */
	  gaiaOutBufferInitialize (&out_buf);
	  gaiaOutSvg (&out_buf, geo, relative, precision);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsSvg1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* called without additional arguments */
    fnct_AsSvg (context, argc, argv, GAIA_SVG_DEFAULT_RELATIVE,
		GAIA_SVG_DEFAULT_PRECISION);
}

static void
fnct_AsSvg2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* called with relative-switch */
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	fnct_AsSvg (context, argc, argv, sqlite3_value_int (argv[1]),
		    GAIA_SVG_DEFAULT_PRECISION);
    else
	sqlite3_result_null (context);
}

static void
fnct_AsSvg3 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* called with relative-switch and precision-argument */
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER
	&& sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	fnct_AsSvg (context, argc, argv, sqlite3_value_int (argv[1]),
		    sqlite3_value_int (argv[2]));
    else
	sqlite3_result_null (context);
}

/* END of Klaus Foerster AsSvg() implementation */

static void
proj_params (sqlite3 * sqlite, int srid, char *proj_params)
{
/* retrives the PROJ params from SPATIAL_SYS_REF table, if possible */
    char sql[256];
    char **results;
    int rows;
    int columns;
    int i;
    int ret;
    char *errMsg = NULL;
    *proj_params = '\0';
    sprintf (sql,
	     "SELECT proj4text FROM spatial_ref_sys WHERE srid = %d", srid);
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, &errMsg);
    if (ret != SQLITE_OK)
      {
	  fprintf (stderr, "unknown SRID: %d\t<%s>\n", srid, errMsg);
	  sqlite3_free (errMsg);
	  return;
      }
    for (i = 1; i <= rows; i++)
	strcpy (proj_params, results[(i * columns)]);
    if (*proj_params == '\0')
	fprintf (stderr, "unknown SRID: %d\n", srid);
    sqlite3_free_table (results);
}

#ifndef OMIT_PROJ		/* PROJ.4 is strictly required to support KML */
static void
fnct_AsKml1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsKml(BLOB encoded geometry [, Integer precision])
/
/ returns the corresponding 'bare geom' KML representation 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geo_wgs84;
    char proj_from[2048];
    char proj_to[2048];
    int precision = 15;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (argc == 2)
      {
	  if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	      precision = sqlite3_value_int (argv[1]);
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    gaiaOutBufferInitialize (&out_buf);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (geo->Srid == 4326)
	      ;			/* already WGS84 */
	  else if (geo->Srid == -1)
	    {
		/* unknown SRID: giving up */
		sqlite3_result_null (context);
		goto stop;
	    }
	  else
	    {
		/* attempting to reproject into WGS84 */
		proj_params (sqlite, geo->Srid, proj_from);
		proj_params (sqlite, 4326, proj_to);
		if (*proj_to == '\0' || *proj_from == '\0')
		  {
		      sqlite3_result_null (context);
		      goto stop;
		  }
		geo_wgs84 = gaiaTransform (geo, proj_from, proj_to);
		if (!geo_wgs84)
		  {
		      sqlite3_result_null (context);
		      goto stop;
		  }
		/* ok, reprojection was successful */
		gaiaFreeGeomColl (geo);
		geo = geo_wgs84;
	    }
	  /* produce KML-notation - actual work is done in gaiageo/gg_wkt.c */
	  gaiaOutBareKml (&out_buf, geo, precision);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
  stop:
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsKml3 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsKml(Anything name, Anything description, BLOB encoded geometry [, Integer precision])
/
/ returns the corresponding 'full' KML representation 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geo_wgs84;
    sqlite3_int64 int_value;
    double dbl_value;
    const char *name;
    const char *desc;
    char *name_malloc = NULL;
    char *desc_malloc = NULL;
    char dummy[128];
    char proj_from[2048];
    char proj_to[2048];
    int precision = 15;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    switch (sqlite3_value_type (argv[0]))
      {
      case SQLITE_TEXT:
	  name = (const char *) sqlite3_value_text (argv[0]);
	  len = strlen (name);
	  name_malloc = malloc (len + 1);
	  strcpy (name_malloc, name);
	  name = name_malloc;
	  break;
      case SQLITE_INTEGER:
	  int_value = sqlite3_value_int64 (argv[0]);
#if defined(_WIN32) || defined(__MINGW32__)
/* CAVEAT: M$ runtime doesn't supports %lld for 64 bits */
	  sprintf (dummy, "%I64d", int_value);
#else
	  sprintf (dummy, "%lld", int_value);
#endif
	  len = strlen (dummy);
	  name_malloc = malloc (len + 1);
	  strcpy (name_malloc, dummy);
	  name = name_malloc;
	  break;
      case SQLITE_FLOAT:
	  dbl_value = sqlite3_value_double (argv[0]);
	  sprintf (dummy, "%1.6f", dbl_value);
	  len = strlen (dummy);
	  name_malloc = malloc (len + 1);
	  strcpy (name_malloc, dummy);
	  name = name_malloc;
	  break;
      case SQLITE_BLOB:
	  name = "BLOB";
	  break;
      default:
	  name = "NULL";
	  break;
      };
    switch (sqlite3_value_type (argv[1]))
      {
      case SQLITE_TEXT:
	  desc = (const char *) sqlite3_value_text (argv[1]);
	  len = strlen (desc);
	  desc_malloc = malloc (len + 1);
	  strcpy (desc_malloc, desc);
	  desc = desc_malloc;
	  break;
      case SQLITE_INTEGER:
	  int_value = sqlite3_value_int64 (argv[1]);
#if defined(_WIN32) || defined(__MINGW32__)
/* CAVEAT: M$ runtime doesn't supports %lld for 64 bits */
	  sprintf (dummy, "%I64d", int_value);
#else
	  sprintf (dummy, "%lld", int_value);
#endif
	  len = strlen (dummy);
	  desc_malloc = malloc (len + 1);
	  strcpy (desc_malloc, dummy);
	  desc = desc_malloc;
	  break;
      case SQLITE_FLOAT:
	  dbl_value = sqlite3_value_double (argv[1]);
	  sprintf (dummy, "%1.6f", dbl_value);
	  len = strlen (dummy);
	  desc_malloc = malloc (len + 1);
	  strcpy (desc_malloc, dummy);
	  desc = desc_malloc;
	  break;
      case SQLITE_BLOB:
	  desc = "BLOB";
	  break;
      default:
	  desc = "NULL";
	  break;
      };
    gaiaOutBufferInitialize (&out_buf);
    if (sqlite3_value_type (argv[2]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  goto stop;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[2]);
    n_bytes = sqlite3_value_bytes (argv[2]);
    if (argc == 4)
      {
	  if (sqlite3_value_type (argv[3]) == SQLITE_INTEGER)
	      precision = sqlite3_value_int (argv[3]);
	  else
	    {
		sqlite3_result_null (context);
		goto stop;
	    }
      }
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (geo->Srid == 4326)
	      ;			/* already WGS84 */
	  else if (geo->Srid == -1)
	    {
		/* unknown SRID: giving up */
		sqlite3_result_null (context);
		goto stop;
	    }
	  else
	    {
		/* attempting to reproject into WGS84 */
		proj_params (sqlite, geo->Srid, proj_from);
		proj_params (sqlite, 4326, proj_to);
		if (*proj_to == '\0' || *proj_from == '\0')
		  {
		      sqlite3_result_null (context);
		      goto stop;
		  }
		geo_wgs84 = gaiaTransform (geo, proj_from, proj_to);
		if (!geo_wgs84)
		  {
		      sqlite3_result_null (context);
		      goto stop;
		  }
		/* ok, reprojection was successful */
		gaiaFreeGeomColl (geo);
		geo = geo_wgs84;
	    }
	  /* produce KML-notation - actual work is done in gaiageo/gg_wkt.c */
	  gaiaOutFullKml (&out_buf, name, desc, geo, precision);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
  stop:
    gaiaFreeGeomColl (geo);
    if (name_malloc)
	free (name_malloc);
    if (desc_malloc)
	free (desc_malloc);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsKml (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsKml(Anything name, Anything description, BLOB encoded geometry)
/     or
/ AsKml(BLOB encoded geometry)
/
/ returns the corresponding KML representation 
/ or NULL if any error is encountered
*/
    if (argc == 3 || argc == 4)
	fnct_AsKml3 (context, argc, argv);
    else
	fnct_AsKml1 (context, argc, argv);
}
#endif /* end including PROJ.4 */

static void
fnct_AsGml (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsGml(BLOB encoded geometry)
/    or
/ AsGml(integer version, BLOB encoded geometry)
/    or
/ AsGml(integer version, BLOB encoded geometry, integer precision)
/
/ *version* may be 2 (GML 2.1.2) or 3 (GML 3.1.1)
/ default *version*: 2
/
/ *precision* is the number of output decimal digits
/ default *precision*: 15
/
/ returns the corresponding GML representation 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int version = 2;
    int precision = 15;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (argc == 3)
      {
	  if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
	      version = sqlite3_value_int (argv[0]);
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
	  if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
	    {
		sqlite3_result_null (context);
		return;
	    }
	  p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
	  n_bytes = sqlite3_value_bytes (argv[1]);
	  if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	      precision = sqlite3_value_int (argv[2]);
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    else if (argc == 2)
      {
	  if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER
	      && sqlite3_value_type (argv[1]) == SQLITE_BLOB)
	    {
		version = sqlite3_value_int (argv[0]);
		p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
		n_bytes = sqlite3_value_bytes (argv[1]);
	    }
	  else if (sqlite3_value_type (argv[0]) == SQLITE_BLOB
		   && sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	    {
		p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
		n_bytes = sqlite3_value_bytes (argv[0]);
		precision = sqlite3_value_int (argv[1]);
	    }
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    else
      {
	  if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
	    {
		sqlite3_result_null (context);
		return;
	    }
	  p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
	  n_bytes = sqlite3_value_bytes (argv[0]);
      }
    gaiaOutBufferInitialize (&out_buf);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (geo->Srid == -1)
	      sqlite3_result_null (context);	/* unknown SRID: giving up */
	  else
	    {
		/* produce GML-notation - actual work is done in gaiageo/gg_wkt.c */
		gaiaOutGml (&out_buf, version, precision, geo);
		if (out_buf.Error || out_buf.Buffer == NULL)
		    sqlite3_result_null (context);
		else
		  {
		      len = out_buf.WriteOffset;
		      sqlite3_result_text (context, out_buf.Buffer, len, free);
		      out_buf.Buffer = NULL;
		  }
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsGeoJSON (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsGeoJSON(BLOB encoded geometry)
/    or
/ AsGeoJSON(BLOB encoded geometry, integer precision)
/    or
/ AsGeoJSON(BLOB encoded geometry, integer precision, integer options)
/
/ *precision* is the number of output decimal digits
/ default *precision*: 15
/
/ *options* may be one of the followings:
/   0 = no options [default]
/   1 = GeoJSON MBR
/   2 = GeoJSON Short CRS (e.g EPSG:4326) 
/   3 = 1 + 2 (Mbr + shortCrs)
/   4 = GeoJSON Long CRS (e.g urn:ogc:def:crs:EPSG::4326)
/   5 = 1 + 4 (Mbr + longCrs)
/
/ returns the corresponding GML representation 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int precision = 15;
    int options = 0;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (argc == 3)
      {
	  if (sqlite3_value_type (argv[0]) == SQLITE_BLOB
	      && sqlite3_value_type (argv[1]) == SQLITE_INTEGER
	      && sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	    {
		p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
		n_bytes = sqlite3_value_bytes (argv[0]);
		precision = sqlite3_value_int (argv[1]);
		options = sqlite3_value_int (argv[2]);
		if (options >= 1 && options <= 5)
		    ;
		else
		    options = 0;
	    }
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    else if (argc == 2)
      {
	  if (sqlite3_value_type (argv[0]) == SQLITE_BLOB
	      && sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	    {
		p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
		n_bytes = sqlite3_value_bytes (argv[0]);
		precision = sqlite3_value_int (argv[1]);
	    }
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    else
      {
	  if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
	    {
		sqlite3_result_null (context);
		return;
	    }
	  p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
	  n_bytes = sqlite3_value_bytes (argv[0]);
      }
    gaiaOutBufferInitialize (&out_buf);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (geo->Srid == -1)
	      sqlite3_result_null (context);	/* unknown SRID: giving up */
	  else
	    {
		/* produce GeoJSON-notation - actual work is done in gaiageo/gg_wkt.c */
		gaiaOutGeoJSON (&out_buf, geo, precision, options);
		if (out_buf.Error || out_buf.Buffer == NULL)
		    sqlite3_result_null (context);
		else
		  {
		      len = out_buf.WriteOffset;
		      sqlite3_result_text (context, out_buf.Buffer, len, free);
		      out_buf.Buffer = NULL;
		  }
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_AsBinary (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsBinary(BLOB encoded geometry)
/
/ returns the corresponding WKB encoded value
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaToWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_AsFGF (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsFGF(BLOB encoded geometry)
/
/ returns the corresponding FGF encoded value
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    int coord_dims;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  fprintf (stderr,
		   "AsFGF() error: argument 2 [geom_coords] is not of the Integer type\n");
	  sqlite3_result_null (context);
	  return;
      }
    coord_dims = sqlite3_value_int (argv[1]);
    if (coord_dims
	== 0 || coord_dims == 1 || coord_dims == 2 || coord_dims == 3)
	;
    else
      {
	  fprintf (stderr,
		   "AsFGF() error: argument 2 [geom_coords] out of range [0,1,2,3]\n");
	  sqlite3_result_null (context);
	  return;
      }
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaToFgf (geo, &p_result, &len, coord_dims);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_MakePoint1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MakePoint(double X, double Y)
/
/ builds a POINT 
/ or NULL if any error is encountered
*/
    int len;
    int int_value;
    unsigned char *p_result = NULL;
    double x;
    double y;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaMakePoint (x, y, -1, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_MakePoint2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MakePoint(double X, double Y, int SRID)
/
/ builds a POINT 
/ or NULL if any error is encountered
*/
    int len;
    int int_value;
    unsigned char *p_result = NULL;
    double x;
    double y;
    int srid;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	srid = sqlite3_value_int (argv[2]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaMakePoint (x, y, srid, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
addGeomPointToDynamicLine (gaiaDynamicLinePtr dyn, gaiaGeomCollPtr geom)
{
/* appending a simple-Point Geometry to a Dynamic Line */
    int pts;
    int lns;
    int pgs;
    gaiaPointPtr pt;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;

    if (dyn == NULL)
	return;
    if (dyn->Error)
	return;
/* checking if GEOM simply is a POINT */
    if (geom == NULL)
      {
	  dyn->Error = 1;
	  return;
      }
    pts = 0;
    lns = 0;
    pgs = 0;
    pt = geom->FirstPoint;
    while (pt)
      {
	  pts++;
	  pt = pt->Next;
      }
    ln = geom->FirstLinestring;
    while (ln)
      {
	  lns++;
	  ln = ln->Next;
      }
    pg = geom->FirstPolygon;
    while (pg)
      {
	  pgs++;
	  pg = pg->Next;
      }
    if (pts == 1 && lns == 0 && pgs == 0)
	;
    else
      {
	  /* failure: not a simple POINT */
	  dyn->Error = 1;
	  return;
      }

    if (dyn->Srid != geom->Srid)
      {
	  /* failure: SRID mismatch */
	  dyn->Error = 1;
	  return;
      }

    switch (geom->FirstPoint->DimensionModel)
      {
      case GAIA_XY_Z_M:
	  gaiaAppendPointZMToDynamicLine (dyn, geom->FirstPoint->X,
					  geom->FirstPoint->Y,
					  geom->FirstPoint->Z,
					  geom->FirstPoint->M);
	  break;
      case GAIA_XY_Z:
	  gaiaAppendPointZToDynamicLine (dyn, geom->FirstPoint->X,
					 geom->FirstPoint->Y,
					 geom->FirstPoint->Z);
	  break;
      case GAIA_XY_M:
	  gaiaAppendPointMToDynamicLine (dyn, geom->FirstPoint->X,
					 geom->FirstPoint->Y,
					 geom->FirstPoint->M);
	  break;
      default:
	  gaiaAppendPointToDynamicLine (dyn, geom->FirstPoint->X,
					geom->FirstPoint->Y);
	  break;
      }
}

static void
fnct_MakeLine_step (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MakeLine(BLOBencoded geom)
/
/ aggregate function - STEP
/
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    gaiaDynamicLinePtr *p;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geom)
	return;
    p = sqlite3_aggregate_context (context, sizeof (gaiaDynamicLinePtr));
    if (!(*p))
      {
	  /* this is the first row */
	  *p = gaiaAllocDynamicLine ();
	  (*p)->Srid = geom->Srid;
	  addGeomPointToDynamicLine (*p, geom);
	  gaiaFreeGeomColl (geom);
      }
    else
      {
	  /* subsequent rows */
	  addGeomPointToDynamicLine (*p, geom);
	  gaiaFreeGeomColl (geom);
      }
}

static gaiaGeomCollPtr
geomFromDynamicLine (gaiaDynamicLinePtr dyn)
{
/* attempting to build a Geometry from a Dynamic Line */
    gaiaGeomCollPtr geom = NULL;
    gaiaLinestringPtr ln = NULL;
    gaiaPointPtr pt;
    int iv;
    int count = 0;
    int dims = GAIA_XY;

    if (dyn == NULL)
	return NULL;
    if (dyn->Error)
	return NULL;

    pt = dyn->First;
    while (pt)
      {
	  /* counting points and checking dims */
	  count++;
	  if (dims == GAIA_XY && pt->DimensionModel != GAIA_XY)
	      dims = pt->DimensionModel;
	  if (dims == GAIA_XY_Z
	      && (pt->DimensionModel == GAIA_XY_M
		  || pt->DimensionModel == GAIA_XY_Z_M))
	      dims = GAIA_XY_Z_M;
	  if (dims == GAIA_XY_M
	      && (pt->DimensionModel == GAIA_XY_Z
		  || pt->DimensionModel == GAIA_XY_Z_M))
	      dims = GAIA_XY_Z_M;
	  pt = pt->Next;
      }
    if (count == 0)
	return NULL;

    switch (dims)
      {
      case GAIA_XY_Z_M:
	  geom = gaiaAllocGeomCollXYZM ();
	  ln = gaiaAllocLinestringXYZM (count);
	  break;
      case GAIA_XY_Z:
	  geom = gaiaAllocGeomCollXYZ ();
	  ln = gaiaAllocLinestringXYZ (count);
	  break;
      case GAIA_XY_M:
	  geom = gaiaAllocGeomCollXYM ();
	  ln = gaiaAllocLinestringXYM (count);
	  break;
      default:
	  geom = gaiaAllocGeomColl ();
	  ln = gaiaAllocLinestring (count);
	  break;
      };

    if (geom != NULL && ln != NULL)
      {
	  gaiaInsertLinestringInGeomColl (geom, ln);
	  geom->Srid = dyn->Srid;
      }
    else
      {
	  if (geom)
	      gaiaFreeGeomColl (geom);
	  if (ln)
	      gaiaFreeLinestring (ln);
	  return NULL;
      }

    iv = 0;
    pt = dyn->First;
    while (pt)
      {
	  /* setting linestring points */
	  if (dims == GAIA_XY_Z_M)
	    {
		gaiaSetPointXYZM (ln->Coords, iv, pt->X, pt->Y, pt->Z, pt->M);
	    }
	  else if (dims == GAIA_XY_Z)
	    {
		gaiaSetPointXYZ (ln->Coords, iv, pt->X, pt->Y, pt->Z);
	    }
	  else if (dims == GAIA_XY_M)
	    {
		gaiaSetPointXYM (ln->Coords, iv, pt->X, pt->Y, pt->M);
	    }
	  else
	    {
		gaiaSetPoint (ln->Coords, iv, pt->X, pt->Y);
	    }
	  iv++;
	  pt = pt->Next;
      }
    return geom;
}

static void
fnct_MakeLine_final (sqlite3_context * context)
{
/* SQL function:
/ MakeLine(BLOBencoded geom)
/
/ aggregate function - FINAL
/
*/
    gaiaGeomCollPtr result;
    gaiaDynamicLinePtr *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    result = geomFromDynamicLine (*p);
    gaiaFreeDynamicLine (*p);
    if (!result)
	sqlite3_result_null (context);
    else
      {
	  /* builds the BLOB geometry to be returned */
	  int len;
	  unsigned char *p_result = NULL;
	  gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
	  gaiaFreeGeomColl (result);
      }
}

static void
fnct_MakeLine (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MakeLine(point-geometry geom1, point-geometry geom2)
/
/ builds a SEGMENT joining two POINTs 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_blob;
    int n_bytes;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  goto stop;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1)
      {
	  sqlite3_result_null (context);
	  goto stop;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  goto stop;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo2)
      {
	  sqlite3_result_null (context);
	  goto stop;
      }
    gaiaMakeLine (geo1, geo2, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
  stop:
    if (geo1)
	gaiaFreeGeomColl (geo1);
    if (geo2)
	gaiaFreeGeomColl (geo2);
}

static void
fnct_Collect_step (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Collect(BLOBencoded geom)
/
/ aggregate function - STEP
/
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr *p;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geom)
	return;
    p = sqlite3_aggregate_context (context, sizeof (gaiaGeomCollPtr));
    if (!(*p))
      {
	  /* this is the first row */
	  *p = geom;
      }
    else
      {
	  /* subsequent rows */
	  result = gaiaMergeGeometries (*p, geom);
	  gaiaFreeGeomColl (*p);
	  *p = result;
	  gaiaFreeGeomColl (geom);
      }
}

static void
fnct_Collect_final (sqlite3_context * context)
{
/* SQL function:
/ Collect(BLOBencoded geom)
/
/ aggregate function - FINAL
/
*/
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    result = *p;
    if (!result)
	sqlite3_result_null (context);
    else if (gaiaIsEmpty (result))
      {
	  gaiaFreeGeomColl (result);
	  sqlite3_result_null (context);
      }
    else
      {
	  /* builds the BLOB geometry to be returned */
	  int len;
	  unsigned char *p_result = NULL;
	  gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
	  gaiaFreeGeomColl (result);
      }
}

static void
fnct_Collect (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Collect(geometry geom1, geometry geom2)
/
/ merges two generic GEOMETRIES into a single one 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaMergeGeometries (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else if (gaiaIsEmpty (result))
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
geom_from_text1 (sqlite3_context * context, int argc, sqlite3_value ** argv,
		 short type)
{
/* SQL function:
/ GeomFromText(WKT encoded geometry)
/
/ returns the current geometry by parsing WKT encoded string 
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, type);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
geom_from_text2 (sqlite3_context * context, int argc, sqlite3_value ** argv,
		 short type)
{
/* SQL function:
/ GeomFromText(WKT encoded geometry, SRID)
/
/ returns the current geometry by parsing WKT encoded string 
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, type);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static int
check_wkb (const unsigned char *wkb, int size, short type)
{
/* checking type coherency for WKB encoded GEOMETRY */
    int little_endian;
    int wkb_type;
    int endian_arch = gaiaEndianArch ();
    if (size < 5)
	return 0;		/* too short to be a WKB */
    if (*(wkb + 0) == 0x01)
	little_endian = GAIA_LITTLE_ENDIAN;
    else if (*(wkb + 0) == 0x00)
	little_endian = GAIA_BIG_ENDIAN;
    else
	return 0;		/* illegal byte ordering; neither BIG-ENDIAN nor LITTLE-ENDIAN */
    wkb_type = gaiaImport32 (wkb + 1, little_endian, endian_arch);
    if (wkb_type == GAIA_POINT || wkb_type == GAIA_LINESTRING
	|| wkb_type == GAIA_POLYGON || wkb_type == GAIA_MULTIPOINT
	|| wkb_type == GAIA_MULTILINESTRING || wkb_type == GAIA_MULTIPOLYGON
	|| wkb_type == GAIA_GEOMETRYCOLLECTION || wkb_type == GAIA_POINTZ
	|| wkb_type == GAIA_LINESTRINGZ || wkb_type == GAIA_POLYGONZ
	|| wkb_type == GAIA_MULTIPOINTZ || wkb_type == GAIA_MULTILINESTRINGZ
	|| wkb_type == GAIA_MULTIPOLYGONZ
	|| wkb_type == GAIA_GEOMETRYCOLLECTIONZ || wkb_type == GAIA_POINTM
	|| wkb_type == GAIA_LINESTRINGM || wkb_type == GAIA_POLYGONM
	|| wkb_type == GAIA_MULTIPOINTM || wkb_type == GAIA_MULTILINESTRINGM
	|| wkb_type == GAIA_MULTIPOLYGONM
	|| wkb_type == GAIA_GEOMETRYCOLLECTIONM || wkb_type == GAIA_POINTZM
	|| wkb_type == GAIA_LINESTRINGZM || wkb_type == GAIA_POLYGONZM
	|| wkb_type == GAIA_MULTIPOINTZM || wkb_type == GAIA_MULTILINESTRINGZM
	|| wkb_type == GAIA_MULTIPOLYGONZM
	|| wkb_type == GAIA_GEOMETRYCOLLECTIONZM)
	;
    else
	return 0;		/* illegal GEOMETRY CLASS */
    if (type < 0)
	;			/* no restrinction about GEOMETRY CLASS TYPE */
    else
      {
	  if (wkb_type != type)
	      return 0;		/* invalid CLASS TYPE for request */
      }
    return 1;
}

static void
geom_from_wkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv,
		short type)
{
/* SQL function:
/ GeomFromWKB(WKB encoded geometry)
/
/ returns the current geometry by parsing a WKB encoded blob 
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    int n_bytes;
    unsigned char *p_result = NULL;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, type))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
geom_from_wkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv,
		short type)
{
/* SQL function:
/ GeomFromWKB(WKB encoded geometry, SRID)
/
/ returns the current geometry by parsing a WKB encoded blob
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    int n_bytes;
    unsigned char *p_result = NULL;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, type))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_GeometryFromFGF1 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromFGF(FGF encoded geometry)
/
/ returns the current geometry by parsing an FGF encoded blob 
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    int n_bytes;
    unsigned char *p_result = NULL;
    const unsigned char *fgf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    fgf = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromFgf (fgf, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_GeometryFromFGF2 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromFGF(FGF encoded geometry, SRID)
/
/ returns the current geometry by parsing an FGF encoded string 
/ or NULL if any error is encountered
/
/ if *type* is a negative value can accept any GEOMETRY CLASS
/ otherwise only requests conforming with required CLASS are valid
*/
    int len;
    int n_bytes;
    unsigned char *p_result = NULL;
    const unsigned char *fgf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    fgf = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromFgf (fgf, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

/*
/ the following functions simply readdress the request to geom_from_text?()
/ setting the appropriate GEOMETRY CLASS TYPE
*/

static void
fnct_GeomFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) -1);
}

static void
fnct_GeomFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) -1);
}

static void
fnct_GeomCollFromText1 (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_GEOMETRYCOLLECTION);
}

static void
fnct_GeomCollFromText2 (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_GEOMETRYCOLLECTION);
}

static void
fnct_LineFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_LINESTRING);
}

static void
fnct_LineFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_LINESTRING);
}

static void
fnct_PointFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_POINT);
}

static void
fnct_PointFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_POINT);
}

static void
fnct_PolyFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_POLYGON);
}

static void
fnct_PolyFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_POLYGON);
}

static void
fnct_MLineFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_MULTILINESTRING);
}

static void
fnct_MLineFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_MULTILINESTRING);
}

static void
fnct_MPointFromText1 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_MULTIPOINT);
}

static void
fnct_MPointFromText2 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_MULTIPOINT);
}

static void
fnct_MPolyFromText1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text1 (context, argc, argv, (short) GAIA_MULTIPOLYGON);
}

static void
fnct_MPolyFromText2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_text2 (context, argc, argv, (short) GAIA_MULTIPOLYGON);
}

/*
/ the following functions simply readdress the request to geom_from_wkb?()
/ setting the appropriate GEOMETRY CLASS TYPE
*/

static void
fnct_GeomFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) -1);
}

static void
fnct_GeomFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) -1);
}

static void
fnct_GeomCollFromWkb1 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_GEOMETRYCOLLECTION);
}

static void
fnct_GeomCollFromWkb2 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_GEOMETRYCOLLECTION);
}

static void
fnct_LineFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_LINESTRING);
}

static void
fnct_LineFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_LINESTRING);
}

static void
fnct_PointFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_POINT);
}

static void
fnct_PointFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_POINT);
}

static void
fnct_PolyFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_POLYGON);
}

static void
fnct_PolyFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_POLYGON);
}

static void
fnct_MLineFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_MULTILINESTRING);
}

static void
fnct_MLineFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_MULTILINESTRING);
}

static void
fnct_MPointFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_MULTIPOINT);
}

static void
fnct_MPointFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_MULTIPOINT);
}

static void
fnct_MPolyFromWkb1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb1 (context, argc, argv, (short) GAIA_MULTIPOLYGON);
}

static void
fnct_MPolyFromWkb2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    geom_from_wkb2 (context, argc, argv, (short) GAIA_MULTIPOLYGON);
}

static void
fnct_CompressGeometry (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ CompressGeometry(BLOB encoded geometry)
/
/ returns a COMPRESSED geometry [if a valid Geometry was supplied]
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaToCompressedBlobWkb (geo, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_UncompressGeometry (sqlite3_context * context, int argc,
			 sqlite3_value ** argv)
{
/* SQL function:
/ UncompressGeometry(BLOB encoded geometry)
/
/ returns an UNCOMPRESSED geometry [if a valid Geometry was supplied] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_SanitizeGeometry (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ SanitizeGeometry(BLOB encoded geometry)
/
/ returns a SANITIZED geometry [if a valid Geometry was supplied]
/ or NULL in any other case
/
/ Sanitizing includes:
/ - repeated vertices suppression
/ - enforcing ring closure
/
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr sanitized = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  sanitized = gaiaSanitize (geo);
	  gaiaToSpatiaLiteBlobWkb (sanitized, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
    gaiaFreeGeomColl (sanitized);
}

static void
cast_count (gaiaGeomCollPtr geom, int *pts, int *lns, int *pgs)
{
/* counting elementary geometries */
    int n_pts = 0;
    int n_lns = 0;
    int n_pgs = 0;
    gaiaPointPtr pt;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;
    if (geom)
      {
	  pt = geom->FirstPoint;
	  while (pt)
	    {
		n_pts++;
		pt = pt->Next;
	    }
	  ln = geom->FirstLinestring;
	  while (ln)
	    {
		n_lns++;
		ln = ln->Next;
	    }
	  pg = geom->FirstPolygon;
	  while (pg)
	    {
		n_pgs++;
		pg = pg->Next;
	    }
      }
    *pts = n_pts;
    *lns = n_lns;
    *pgs = n_pgs;
}

static void
fnct_CastToPoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToPoint(BLOB encoded geometry)
/
/ returns a POINT-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts == 1 && lns == 0 && pgs == 0)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_POINT;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToLinestring (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ CastToLinestring(BLOB encoded geometry)
/
/ returns a LINESTRING-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts == 0 && lns == 1 && pgs == 0)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_LINESTRING;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToPolygon (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToPolygon(BLOB encoded geometry)
/
/ returns a POLYGON-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts == 0 && lns == 0 && pgs == 1)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_POLYGON;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToMultiPoint (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ CastToMultiPoint(BLOB encoded geometry)
/
/ returns a MULTIPOINT-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts >= 1 && lns == 0 && pgs == 0)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTIPOINT;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToMultiLinestring (sqlite3_context * context, int argc,
			    sqlite3_value ** argv)
{
/* SQL function:
/ CastToMultiLinestring(BLOB encoded geometry)
/
/ returns a MULTILINESTRING-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts == 0 && lns >= 1 && pgs == 0)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTILINESTRING;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToMultiPolygon (sqlite3_context * context, int argc,
			 sqlite3_value ** argv)
{
/* SQL function:
/ CastToMultiPolygon(BLOB encoded geometry)
/
/ returns a MULTIPOLYGON-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts == 0 && lns == 0 && pgs >= 1)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTIPOLYGON;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToGeometryCollection (sqlite3_context * context, int argc,
			       sqlite3_value ** argv)
{
/* SQL function:
/ CastToGeometryCollection(BLOB encoded geometry)
/
/ returns a GEOMETRYCOLLECTION-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts >= 1 || lns >= 1 || pgs >= 1)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_GEOMETRYCOLLECTION;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToMulti (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToMulti(BLOB encoded geometry)
/
/ returns a MULTIPOINT, MULTILINESTRING, MULTIPOLYGON or
/ GEOMETRYCOLLECTION-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts >= 1 || lns >= 1 || pgs >= 1)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		if (pts >= 1 && lns == 0 && pgs == 0)
		    geom2->DeclaredType = GAIA_MULTIPOINT;
		else if (pts == 0 && lns >= 1 && pgs == 0)
		    geom2->DeclaredType = GAIA_MULTILINESTRING;
		else if (pts == 0 && lns == 0 && pgs >= 1)
		    geom2->DeclaredType = GAIA_MULTIPOLYGON;
		else
		    geom2->DeclaredType = GAIA_GEOMETRYCOLLECTION;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToSingle (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToSingle(BLOB encoded geometry)
/
/ returns a POINT, LINESTRING or POLYGON-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    int ok;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  ok = 0;
	  if (pts == 1 && lns == 0 && pgs == 0)
	      ok = 1;
	  if (pts == 0 && lns == 1 && pgs == 0)
	      ok = 1;
	  if (pts == 0 && lns == 0 && pgs == 1)
	      ok = 1;
	  if (ok)
	    {
		geom2 = gaiaCloneGeomColl (geo);
		geom2->Srid = geo->Srid;
		if (pts == 1)
		    geom2->DeclaredType = GAIA_POINT;
		else if (lns == 1)
		    geom2->DeclaredType = GAIA_LINESTRING;
		else
		    geom2->DeclaredType = GAIA_POLYGON;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToXY (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToXY(BLOB encoded geometry)
/
/ returns an XY-dimension Geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  geom2 = gaiaCastGeomCollToXY (geo);
	  if (geom2)
	    {
		geom2->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToXYZ (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToXY(BLOB encoded geometry)
/
/ returns an XY-dimension Geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  geom2 = gaiaCastGeomCollToXYZ (geo);
	  if (geom2)
	    {
		geom2->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToXYM (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToXY(BLOB encoded geometry)
/
/ returns an XYM-dimension Geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  geom2 = gaiaCastGeomCollToXYM (geo);
	  if (geom2)
	    {
		geom2->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CastToXYZM (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CastToXY(BLOB encoded geometry)
/
/ returns an XYZM-dimension Geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  geom2 = gaiaCastGeomCollToXYZM (geo);
	  if (geom2)
	    {
		geom2->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ExtractMultiPoint (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ ExtractMultiPoint(BLOB encoded geometry)
/
/ returns a MULTIPOINT-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pts >= 1)
	    {
		geom2 = gaiaCloneGeomCollPoints (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTIPOINT;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ExtractMultiLinestring (sqlite3_context * context, int argc,
			     sqlite3_value ** argv)
{
/* SQL function:
/ ExtractMultiLinestring(BLOB encoded geometry)
/
/ returns a MULTILINESTRING-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (lns >= 1)
	    {
		geom2 = gaiaCloneGeomCollLinestrings (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTILINESTRING;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ExtractMultiPolygon (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ ExtractMultiPolygon(BLOB encoded geometry)
/
/ returns a MULTIPOLYGON-type geometry [if conversion is possible] 
/ or NULL in any other case
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    int pts;
    int lns;
    int pgs;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  cast_count (geo, &pts, &lns, &pgs);
	  if (pgs >= 1)
	    {
		geom2 = gaiaCloneGeomCollPolygons (geo);
		geom2->Srid = geo->Srid;
		geom2->DeclaredType = GAIA_MULTIPOLYGON;
		gaiaToSpatiaLiteBlobWkb (geom2, &p_result, &len);
		gaiaFreeGeomColl (geom2);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Dimension (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Dimension(BLOB encoded geometry)
/
/ returns:
/ 0 if geometry is a POINT or MULTIPOINT
/ 1 if geometry is a LINESTRING or MULTILINESTRING
/ 2 if geometry is a POLYGON or MULTIPOLYGON
/ 0, 1, 2, for GEOMETRYCOLLECTIONS according to geometries contained inside
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int dim;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  dim = gaiaDimension (geo);
	  sqlite3_result_int (context, dim);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CoordDimension (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CoordDimension(BLOB encoded geometry)
/
/ returns:
/ 'XY', 'XYM', 'XYZ', 'XYZM'
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    char *p_dim = NULL;
    char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (geo->DimensionModel == GAIA_XY)
	      p_dim = "XY";
	  else if (geo->DimensionModel == GAIA_XY_Z)
	      p_dim = "XYZ";
	  else if (geo->DimensionModel == GAIA_XY_M)
	      p_dim = "XYM";
	  else if (geo->DimensionModel == GAIA_XY_Z_M)
	      p_dim = "XYZM";
	  if (p_dim)
	    {
		len = strlen (p_dim);
		p_result = malloc (len + 1);
		strcpy (p_result, p_dim);
	    }
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	    {
		len = strlen (p_result);
		sqlite3_result_text (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_GeometryType (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeometryType(BLOB encoded geometry)
/
/ returns the class for current geometry:
/ 'POINT' or 'MULTIPOINT' [Z, M, ZM]
/ 'LINESTRING' or 'MULTILINESTRING' [Z, M, ZM]
/ 'POLYGON' or 'MULTIPOLYGON' [Z, M, ZM]
/ 'GEOMETRYCOLLECTION'  [Z, M, ZM]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int type;
    char *p_type = NULL;
    char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  type = gaiaGeometryType (geo);
	  switch (type)
	    {
	    case GAIA_POINT:
		p_type = "POINT";
		break;
	    case GAIA_POINTZ:
		p_type = "POINT Z";
		break;
	    case GAIA_POINTM:
		p_type = "POINT M";
		break;
	    case GAIA_POINTZM:
		p_type = "POINT ZM";
		break;
	    case GAIA_MULTIPOINT:
		p_type = "MULTIPOINT";
		break;
	    case GAIA_MULTIPOINTZ:
		p_type = "MULTIPOINT Z";
		break;
	    case GAIA_MULTIPOINTM:
		p_type = "MULTIPOINT M";
		break;
	    case GAIA_MULTIPOINTZM:
		p_type = "MULTIPOINT ZM";
		break;
	    case GAIA_LINESTRING:
	    case GAIA_COMPRESSED_LINESTRING:
		p_type = "LINESTRING";
		break;
	    case GAIA_LINESTRINGZ:
	    case GAIA_COMPRESSED_LINESTRINGZ:
		p_type = "LINESTRING Z";
		break;
	    case GAIA_LINESTRINGM:
	    case GAIA_COMPRESSED_LINESTRINGM:
		p_type = "LINESTRING M";
		break;
	    case GAIA_LINESTRINGZM:
	    case GAIA_COMPRESSED_LINESTRINGZM:
		p_type = "LINESTRING ZM";
		break;
	    case GAIA_MULTILINESTRING:
		p_type = "MULTILINESTRING";
		break;
	    case GAIA_MULTILINESTRINGZ:
		p_type = "MULTILINESTRING Z";
		break;
	    case GAIA_MULTILINESTRINGM:
		p_type = "MULTILINESTRING M";
		break;
	    case GAIA_MULTILINESTRINGZM:
		p_type = "MULTILINESTRING ZM";
		break;
	    case GAIA_POLYGON:
	    case GAIA_COMPRESSED_POLYGON:
		p_type = "POLYGON";
		break;
	    case GAIA_POLYGONZ:
	    case GAIA_COMPRESSED_POLYGONZ:
		p_type = "POLYGON Z";
		break;
	    case GAIA_POLYGONM:
	    case GAIA_COMPRESSED_POLYGONM:
		p_type = "POLYGON M";
		break;
	    case GAIA_POLYGONZM:
	    case GAIA_COMPRESSED_POLYGONZM:
		p_type = "POLYGON ZM";
		break;
	    case GAIA_MULTIPOLYGON:
		p_type = "MULTIPOLYGON";
		break;
	    case GAIA_MULTIPOLYGONZ:
		p_type = "MULTIPOLYGON Z";
		break;
	    case GAIA_MULTIPOLYGONM:
		p_type = "MULTIPOLYGON M";
		break;
	    case GAIA_MULTIPOLYGONZM:
		p_type = "MULTIPOLYGON ZM";
		break;
	    case GAIA_GEOMETRYCOLLECTION:
		p_type = "GEOMETRYCOLLECTION";
		break;
	    case GAIA_GEOMETRYCOLLECTIONZ:
		p_type = "GEOMETRYCOLLECTION Z";
		break;
	    case GAIA_GEOMETRYCOLLECTIONM:
		p_type = "GEOMETRYCOLLECTION M";
		break;
	    case GAIA_GEOMETRYCOLLECTIONZM:
		p_type = "GEOMETRYCOLLECTION ZM";
		break;
	    };
	  if (p_type)
	    {
		len = strlen (p_type);
		p_result = malloc (len + 1);
		strcpy (p_result, p_type);
	    }
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	    {
		len = strlen (p_result);
		sqlite3_result_text (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_GeometryAliasType (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ GeometryAliasType(BLOB encoded geometry)
/
/ returns the alias-class for current geometry:
/ 'POINT'
/ 'LINESTRING'
/ 'POLYGON'
/ 'MULTIPOINT'
/ 'MULTILINESTRING'
/ 'MULTIPOLYGON'
/ 'GEOMETRYCOLLECTION' 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int type;
    char *p_type = NULL;
    char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  type = gaiaGeometryAliasType (geo);
	  switch (type)
	    {
	    case GAIA_POINT:
		p_type = "POINT";
		break;
	    case GAIA_MULTIPOINT:
		p_type = "MULTIPOINT";
		break;
	    case GAIA_LINESTRING:
		p_type = "LINESTRING";
		break;
	    case GAIA_MULTILINESTRING:
		p_type = "MULTILINESTRING";
		break;
	    case GAIA_POLYGON:
		p_type = "POLYGON";
		break;
	    case GAIA_MULTIPOLYGON:
		p_type = "MULTIPOLYGON";
		break;
	    case GAIA_GEOMETRYCOLLECTION:
		p_type = "GEOMETRYCOLLECTION";
		break;
	    };
	  if (p_type)
	    {
		len = strlen (p_type);
		p_result = malloc (len + 1);
		strcpy (p_result, p_type);
	    }
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	    {
		len = strlen (p_result);
		sqlite3_result_text (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_SridFromAuthCRS (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ SridFromAuthCRS(auth_name, auth_srid)
/
/ returns the SRID
/ or NULL if any error is encountered
*/
    const unsigned char *auth_name;
    int auth_srid;
    int srid = -1;
    char sql[1024];
    char sql2[1024];
    char **results;
    int n_rows;
    int n_columns;
    char *err_msg = NULL;
    int ret;
    int i;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    auth_name = sqlite3_value_text (argv[0]);
    auth_srid = sqlite3_value_int (argv[1]);

    sprintf (sql, "SELECT srid FROM spatial_ref_sys ");
    sprintf (sql2, "WHERE auth_name LIKE '%s' AND auth_srid = %d", auth_name,
	     auth_srid);
    strcat (sql, sql2);
    ret =
	sqlite3_get_table (sqlite, sql, &results, &n_rows, &n_columns,
			   &err_msg);
    if (ret != SQLITE_OK)
	goto done;
    if (n_rows >= 1)
      {
	  for (i = 1; i <= n_rows; i++)
	      srid = atoi (results[(i * n_columns) + 0]);
      }
    sqlite3_free_table (results);
  done:
    sqlite3_result_int (context, srid);
}

static void
fnct_SRID (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Srid(BLOB encoded geometry)
/
/ returns the SRID
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
	sqlite3_result_int (context, geo->Srid);
    gaiaFreeGeomColl (geo);
}

static void
fnct_SetSRID (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ SetSrid(BLOBencoded geometry, srid)
/
/ returns a new geometry that is the original one received, but with the new SRID [no coordinates translation is applied]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    int srid;
    unsigned char *p_result = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	srid = sqlite3_value_int (argv[1]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  geo->Srid = srid;
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &n_bytes);
	  sqlite3_result_blob (context, p_result, n_bytes, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_IsEmpty (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ IsEmpty(BLOB encoded geometry)
/
/ returns:
/ 1 if this geometry contains no elementary geometries
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_int (context, 1);
    else
	sqlite3_result_int (context, gaiaIsEmpty (geo));
    gaiaFreeGeomColl (geo);
}

static void
fnct_Envelope (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Envelope(BLOB encoded geometry)
/
/ returns the MBR for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr bbox;
    gaiaPolygonPtr polyg;
    gaiaRingPtr rect;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (gaiaIsEmpty (geo))
	      sqlite3_result_null (context);
	  else
	    {
		gaiaMbrGeometry (geo);
		bbox = gaiaAllocGeomColl ();
		bbox->Srid = geo->Srid;
		polyg = gaiaAddPolygonToGeomColl (bbox, 5, 0);
		rect = polyg->Exterior;
		gaiaSetPoint (rect->Coords, 0, geo->MinX, geo->MinY);	/* vertex # 1 */
		gaiaSetPoint (rect->Coords, 1, geo->MaxX, geo->MinY);	/* vertex # 2 */
		gaiaSetPoint (rect->Coords, 2, geo->MaxX, geo->MaxY);	/* vertex # 3 */
		gaiaSetPoint (rect->Coords, 3, geo->MinX, geo->MaxY);	/* vertex # 4 */
		gaiaSetPoint (rect->Coords, 4, geo->MinX, geo->MinY);	/* vertex # 5 [same as vertex # 1 to close the polygon] */
		gaiaToSpatiaLiteBlobWkb (bbox, &p_result, &len);
		gaiaFreeGeomColl (bbox);
		sqlite3_result_blob (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
build_filter_mbr (sqlite3_context * context, int argc,
		  sqlite3_value ** argv, int mode)
{
/* SQL functions:
/ BuildMbrFilter(double X1, double Y1, double X2, double Y2)
/ FilterMBRWithin(double X1, double Y1, double X2, double Y2)
/ FilterMBRContain(double X1, double Y1, double X2, double Y2)
/ FilterMBRIntersects(double X1, double Y1, double X2, double Y2)
/
/ builds a generic filter for MBR from two points (identifying a rectangle's diagonal) 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    double x1;
    double y1;
    double x2;
    double y2;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x1 = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y1 = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	x2 = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  x2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[3]) == SQLITE_FLOAT)
	y2 = sqlite3_value_double (argv[3]);
    else if (sqlite3_value_type (argv[3]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[3]);
	  y2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaBuildFilterMbr (x1, y1, x2, y2, mode, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

/*
/ the following functions simply readdress the request to build_filter_mbr()
/ setting the appropriate MODe
*/

static void
fnct_BuildMbrFilter (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    build_filter_mbr (context, argc, argv, GAIA_FILTER_MBR_DECLARE);
}

static void
fnct_FilterMbrWithin (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
    build_filter_mbr (context, argc, argv, GAIA_FILTER_MBR_WITHIN);
}

static void
fnct_FilterMbrContains (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
    build_filter_mbr (context, argc, argv, GAIA_FILTER_MBR_CONTAINS);
}

static void
fnct_FilterMbrIntersects (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
    build_filter_mbr (context, argc, argv, GAIA_FILTER_MBR_INTERSECTS);
}

static void
fnct_BuildMbr1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ BuildMBR(double X1, double Y1, double X2, double Y2)
/
/ builds an MBR from two points (identifying a rectangle's diagonal) 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    double x1;
    double y1;
    double x2;
    double y2;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x1 = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y1 = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	x2 = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  x2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[3]) == SQLITE_FLOAT)
	y2 = sqlite3_value_double (argv[3]);
    else if (sqlite3_value_type (argv[3]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[3]);
	  y2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaBuildMbr (x1, y1, x2, y2, -1, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_BuildMbr2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ BuildMBR(double X1, double Y1, double X2, double Y2, int SRID)
/
/ builds an MBR from two points (identifying a rectangle's diagonal) 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    double x1;
    double y1;
    double x2;
    double y2;
    int int_value;
    int srid;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x1 = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y1 = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	x2 = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  x2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[3]) == SQLITE_FLOAT)
	y2 = sqlite3_value_double (argv[3]);
    else if (sqlite3_value_type (argv[3]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[3]);
	  y2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[4]) == SQLITE_INTEGER)
	srid = sqlite3_value_int (argv[4]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaBuildMbr (x1, y1, x2, y2, srid, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_BuildCircleMbr1 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BuildCircleMBR(double X, double Y, double radius)
/
/ builds an MBR from two points (identifying a rectangle's diagonal) 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    double x;
    double y;
    double radius;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	radius = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  radius = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaBuildCircleMbr (x, y, radius, -1, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_BuildCircleMbr2 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BuildCircleMBR(double X, double Y, double radius, int SRID)
/
/ builds an MBR from two points (identifying a rectangle's diagonal) 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    double x;
    double y;
    double radius;
    int int_value;
    int srid;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	radius = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  radius = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[3]) == SQLITE_INTEGER)
	srid = sqlite3_value_int (argv[3]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaBuildCircleMbr (x, y, radius, srid, &p_result, &len);
    if (!p_result)
	sqlite3_result_null (context);
    else
	sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_Extent_step (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Extent(BLOBencoded geom)
/
/ aggregate function - STEP
/
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    double **p;
    double *max_min;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geom)
	return;
    gaiaMbrGeometry (geom);
    p = sqlite3_aggregate_context (context, sizeof (double **));
    if (!(*p))
      {
	  /* this is the first row */
	  max_min = malloc (sizeof (double) * 4);
	  *(max_min + 0) = geom->MinX;
	  *(max_min + 1) = geom->MinY;
	  *(max_min + 2) = geom->MaxX;
	  *(max_min + 3) = geom->MaxY;
	  *p = max_min;
      }
    else
      {
	  /* subsequent rows */
	  max_min = *p;
	  if (geom->MinX < *(max_min + 0))
	      *(max_min + 0) = geom->MinX;
	  if (geom->MinY < *(max_min + 1))
	      *(max_min + 1) = geom->MinY;
	  if (geom->MaxX > *(max_min + 2))
	      *(max_min + 2) = geom->MaxX;
	  if (geom->MaxY > *(max_min + 3))
	      *(max_min + 3) = geom->MaxY;
      }
    gaiaFreeGeomColl (geom);
}

static void
fnct_Extent_final (sqlite3_context * context)
{
/* SQL function:
/ Extent(BLOBencoded geom)
/
/ aggregate function - FINAL
/
*/
    gaiaGeomCollPtr result;
    gaiaPolygonPtr polyg;
    gaiaRingPtr rect;
    double *max_min;
    double minx;
    double miny;
    double maxx;
    double maxy;
    double **p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    max_min = *p;
    if (!max_min)
      {
	  sqlite3_result_null (context);
	  return;
      }
    result = gaiaAllocGeomColl ();
    if (!result)
	sqlite3_result_null (context);
    else
      {
	  /* builds the BLOB geometry to be returned */
	  int len;
	  unsigned char *p_result = NULL;
	  polyg = gaiaAddPolygonToGeomColl (result, 5, 0);
	  rect = polyg->Exterior;
	  minx = *(max_min + 0);
	  miny = *(max_min + 1);
	  maxx = *(max_min + 2);
	  maxy = *(max_min + 3);
	  gaiaSetPoint (rect->Coords, 0, minx, miny);	/* vertex # 1 */
	  gaiaSetPoint (rect->Coords, 1, maxx, miny);	/* vertex # 2 */
	  gaiaSetPoint (rect->Coords, 2, maxx, maxy);	/* vertex # 3 */
	  gaiaSetPoint (rect->Coords, 3, minx, maxy);	/* vertex # 4 */
	  gaiaSetPoint (rect->Coords, 4, minx, miny);	/* vertex # 5 [same as vertex # 1 to close the polygon] */
	  gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
	  gaiaFreeGeomColl (result);
      }
    free (max_min);
}

static void
fnct_MbrMinX (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MbrMinX(BLOB encoded GEMETRY)
/
/ returns the MinX coordinate for current geometry's MBR 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double coord;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!gaiaGetMbrMinX (p_blob, n_bytes, &coord))
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, coord);
}

static void
fnct_MbrMaxX (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MbrMaxX(BLOB encoded GEMETRY)
/
/ returns the MaxX coordinate for current geometry's MBR 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double coord;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!gaiaGetMbrMaxX (p_blob, n_bytes, &coord))
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, coord);
}

static void
fnct_MbrMinY (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MbrMinY(BLOB encoded GEMETRY)
/
/ returns the MinY coordinate for current geometry's MBR 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double coord;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!gaiaGetMbrMinY (p_blob, n_bytes, &coord))
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, coord);
}

static void
fnct_MbrMaxY (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ MbrMaxY(BLOB encoded GEMETRY)
/
/ returns the MaxY coordinate for current geometry's MBR 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double coord;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!gaiaGetMbrMaxY (p_blob, n_bytes, &coord))
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, coord);
}

#ifndef OMIT_GEOCALLBACKS	/* supporting RTree geometry callbacks */
static void
gaia_mbr_del (void *p)
{
/* freeing data used by R*Tree Geometry Callback */
    sqlite3_free (p);
}

static int
fnct_RTreeIntersects (sqlite3_rtree_geometry * p, int nCoord, double *aCoord,
		      int *pRes)
{
/* R*Tree Geometry callback function:
/ ... MATCH RTreeIntersects(double x1, double y1, double x2, double y2)
*/
    struct gaia_rtree_mbr *mbr;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    float fminx;
    float fminy;
    float fmaxx;
    float fmaxy;
    double tic;
    double tic2;

    if (p->pUser == 0)
      {
	  /* first call: we must check args and then initialize the MBR struct */
	  if (nCoord != 4)
	      return SQLITE_ERROR;
	  if (p->nParam != 4)
	      return SQLITE_ERROR;
	  mbr = (struct gaia_rtree_mbr *) (p->pUser =
					   sqlite3_malloc (sizeof
							   (struct
							    gaia_rtree_mbr)));
	  if (!mbr)
	      return SQLITE_NOMEM;
	  p->xDelUser = gaia_mbr_del;
	  xmin = p->aParam[0];
	  ymin = p->aParam[1];
	  xmax = p->aParam[2];
	  ymax = p->aParam[3];
	  if (xmin > xmax)
	    {
		xmin = p->aParam[2];
		xmax = p->aParam[0];
	    }
	  if (ymin > ymax)
	    {
		ymin = p->aParam[3];
		ymax = p->aParam[1];
	    }

	  /* adjusting the MBR so to compensate for DOUBLE/FLOAT truncations */
	  fminx = (float) xmin;
	  fminy = (float) ymin;
	  fmaxx = (float) xmax;
	  fmaxy = (float) ymax;
	  tic = fabs (xmin - fminx);
	  tic2 = fabs (ymin - fminy);
	  if (tic2 > tic)
	      tic = tic2;
	  tic2 = fabs (xmax - fmaxx);
	  if (tic2 > tic)
	      tic = tic2;
	  tic2 = fabs (ymax - fmaxy);
	  if (tic2 > tic)
	      tic = tic2;
	  tic *= 2.0;

	  mbr->minx = xmin - tic;
	  mbr->miny = ymin - tic;
	  mbr->maxx = xmax + tic;
	  mbr->maxy = ymax + tic;
      }

    mbr = (struct gaia_rtree_mbr *) (p->pUser);
    xmin = aCoord[0];
    xmax = aCoord[1];
    ymin = aCoord[2];
    ymax = aCoord[3];
    *pRes = 1;
/* evaluating Intersects relationship */
    if (xmin > mbr->maxx)
	*pRes = 0;
    if (xmax < mbr->minx)
	*pRes = 0;
    if (ymin > mbr->maxy)
	*pRes = 0;
    if (ymax < mbr->miny)
	*pRes = 0;
    return SQLITE_OK;
}

static int
fnct_RTreeDistWithin (sqlite3_rtree_geometry * p, int nCoord, double *aCoord,
		      int *pRes)
{
/* R*Tree Geometry callback function:
/ ... MATCH RTreeDistWithin(double x, double y, double radius)
*/
    struct gaia_rtree_mbr *mbr;
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    if (p->pUser == 0)
      {
	  /* first call: we must check args and then initialize the MBR struct */
	  if (nCoord != 4)
	      return SQLITE_ERROR;
	  if (p->nParam != 3)
	      return SQLITE_ERROR;
	  mbr = (struct gaia_rtree_mbr *) (p->pUser =
					   sqlite3_malloc (sizeof
							   (struct
							    gaia_rtree_mbr)));
	  if (!mbr)
	      return SQLITE_NOMEM;
	  p->xDelUser = gaia_mbr_del;
	  mbr->minx = p->aParam[0] - p->aParam[2];
	  mbr->miny = p->aParam[1] - p->aParam[2];
	  mbr->maxx = p->aParam[0] + p->aParam[2];
	  mbr->maxy = p->aParam[1] + p->aParam[2];
      }

    mbr = (struct gaia_rtree_mbr *) (p->pUser);
    xmin = aCoord[0];
    xmax = aCoord[1];
    ymin = aCoord[2];
    ymax = aCoord[3];
    *pRes = 1;
/* evaluating Intersects relationship */
    if (xmin > mbr->maxx)
	*pRes = 0;
    if (xmax < mbr->minx)
	*pRes = 0;
    if (ymin > mbr->maxy)
	*pRes = 0;
    if (ymax < mbr->miny)
	*pRes = 0;
    return SQLITE_OK;
}
#endif /* end RTree geometry callbacks */

static void
fnct_X (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ X(BLOB encoded POINT)
/
/ returns the X coordinate for current POINT geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaPointPtr point;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = simplePoint (geo);
	  if (!point)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, point->X);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Y (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Y(BLOB encoded POINT)
/
/ returns the Y coordinate for current POINT geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaPointPtr point;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = simplePoint (geo);
	  if (!point)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, point->Y);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Z (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Z(BLOB encoded POINT)
/
/ returns the Z coordinate for current POINT geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaPointPtr point;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = simplePoint (geo);
	  if (!point)
	      sqlite3_result_null (context);
	  else
	    {
		if (point->DimensionModel == GAIA_XY_Z
		    || point->DimensionModel == GAIA_XY_Z_M)
		    sqlite3_result_double (context, point->Z);
		else
		    sqlite3_result_null (context);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_M (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ M(BLOB encoded POINT)
/
/ returns the M coordinate for current POINT geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaPointPtr point;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = simplePoint (geo);
	  if (!point)
	      sqlite3_result_null (context);
	  else
	    {
		if (point->DimensionModel == GAIA_XY_M
		    || point->DimensionModel == GAIA_XY_Z_M)
		    sqlite3_result_double (context, point->M);
		else
		    sqlite3_result_null (context);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_NumPoints (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ NumPoints(BLOB encoded LINESTRING)
/
/ returns the number of vertices for current LINESTRING geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  line = simpleLinestring (geo);
	  if (!line)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_int (context, line->Points);
      }
    gaiaFreeGeomColl (geo);
}

static void
point_n (sqlite3_context * context, int argc, sqlite3_value ** argv,
	 int request)
{
/* SQL functions:
/ StartPoint(BLOB encoded LINESTRING geometry)
/ EndPoint(BLOB encoded LINESTRING geometry)
/ PointN(BLOB encoded LINESTRING geometry, integer point_no)
/
/ returns the Nth POINT for current LINESTRING geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int vertex;
    int len;
    double x;
    double y;
    double z;
    double m;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (request == GAIA_POINTN)
      {
	  /* PointN() requires point index to be defined as an SQL function argument */
	  if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
	    {
		sqlite3_result_null (context);
		return;
	    }
	  vertex = sqlite3_value_int (argv[1]);
      }
    else if (request == GAIA_END_POINT)
	vertex = -1;		/* EndPoint() specifies a negative point index */
    else
	vertex = 1;		/* StartPoint() */
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  line = simpleLinestring (geo);
	  if (!line)
	      sqlite3_result_null (context);
	  else
	    {
		if (vertex < 0)
		    vertex = line->Points - 1;
		else
		    vertex -= 1;	/* decreasing the point index by 1, because PointN counts starting at index 1 */
		if (vertex >= 0 && vertex < line->Points)
		  {
		      if (line->DimensionModel == GAIA_XY_Z)
			{
			    gaiaGetPointXYZ (line->Coords, vertex, &x, &y, &z);
			    result = gaiaAllocGeomCollXYZ ();
			    result->Srid = geo->Srid;
			    gaiaAddPointToGeomCollXYZ (result, x, y, z);
			}
		      else if (line->DimensionModel == GAIA_XY_M)
			{
			    gaiaGetPointXYM (line->Coords, vertex, &x, &y, &m);
			    result = gaiaAllocGeomCollXYM ();
			    result->Srid = geo->Srid;
			    gaiaAddPointToGeomCollXYM (result, x, y, m);
			}
		      else if (line->DimensionModel == GAIA_XY_Z_M)
			{
			    gaiaGetPointXYZM (line->Coords, vertex, &x, &y, &z,
					      &m);
			    result = gaiaAllocGeomCollXYZM ();
			    result->Srid = geo->Srid;
			    gaiaAddPointToGeomCollXYZM (result, x, y, z, m);
			}
		      else
			{
			    gaiaGetPoint (line->Coords, vertex, &x, &y);
			    result = gaiaAllocGeomColl ();
			    result->Srid = geo->Srid;
			    gaiaAddPointToGeomColl (result, x, y);
			}
		  }
		else
		    result = NULL;
		if (!result)
		    sqlite3_result_null (context);
		else
		  {
		      gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		      gaiaFreeGeomColl (result);
		      sqlite3_result_blob (context, p_result, len, free);
		  }
	    }
      }
    gaiaFreeGeomColl (geo);
}

/*
/ the following functions simply readdress the request to point_n()
/ setting the appropriate request mode
*/

static void
fnct_StartPoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    point_n (context, argc, argv, GAIA_START_POINT);
}

static void
fnct_EndPoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    point_n (context, argc, argv, GAIA_END_POINT);
}

static void
fnct_PointN (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    point_n (context, argc, argv, GAIA_POINTN);
}

static void
fnct_ExteriorRing (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL functions:
/ ExteriorRing(BLOB encoded POLYGON geometry)
/
/ returns the EXTERIOR RING for current POLYGON geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int iv;
    double x;
    double y;
    double z;
    double m;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    gaiaPolygonPtr polyg;
    gaiaRingPtr ring;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  polyg = simplePolygon (geo);
	  if (!polyg)
	      sqlite3_result_null (context);
	  else
	    {
		ring = polyg->Exterior;
		if (ring->DimensionModel == GAIA_XY_Z)
		    result = gaiaAllocGeomCollXYZ ();
		else if (ring->DimensionModel == GAIA_XY_M)
		    result = gaiaAllocGeomCollXYM ();
		else if (ring->DimensionModel == GAIA_XY_Z_M)
		    result = gaiaAllocGeomCollXYZM ();
		else
		    result = gaiaAllocGeomColl ();
		result->Srid = geo->Srid;
		line = gaiaAddLinestringToGeomColl (result, ring->Points);
		for (iv = 0; iv < line->Points; iv++)
		  {
		      if (ring->DimensionModel == GAIA_XY_Z)
			{
			    gaiaGetPointXYZ (ring->Coords, iv, &x, &y, &z);
			    gaiaSetPointXYZ (line->Coords, iv, x, y, z);
			}
		      else if (ring->DimensionModel == GAIA_XY_M)
			{
			    gaiaGetPointXYM (ring->Coords, iv, &x, &y, &m);
			    gaiaSetPointXYM (line->Coords, iv, x, y, m);
			}
		      else if (ring->DimensionModel == GAIA_XY_Z_M)
			{
			    gaiaGetPointXYZM (ring->Coords, iv, &x, &y, &z, &m);
			    gaiaSetPointXYZM (line->Coords, iv, x, y, z, m);
			}
		      else
			{
			    gaiaGetPoint (ring->Coords, iv, &x, &y);
			    gaiaSetPoint (line->Coords, iv, x, y);
			}
		  }
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		gaiaFreeGeomColl (result);
		sqlite3_result_blob (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_NumInteriorRings (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ NumInteriorRings(BLOB encoded POLYGON)
/
/ returns the number of INTERIOR RINGS for current POLYGON geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaPolygonPtr polyg;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  polyg = simplePolygon (geo);
	  if (!polyg)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_int (context, polyg->NumInteriors);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_InteriorRingN (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL functions:
/ InteriorRingN(BLOB encoded POLYGON geometry)
/
/ returns the Nth INTERIOR RING for current POLYGON geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int border;
    int iv;
    double x;
    double y;
    double z;
    double m;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    gaiaPolygonPtr polyg;
    gaiaRingPtr ring;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    border = sqlite3_value_int (argv[1]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  polyg = simplePolygon (geo);
	  if (!polyg)
	      sqlite3_result_null (context);
	  else
	    {
		if (border >= 1 && border <= polyg->NumInteriors)
		  {
		      ring = polyg->Interiors + (border - 1);
		      if (ring->DimensionModel == GAIA_XY_Z)
			  result = gaiaAllocGeomCollXYZ ();
		      else if (ring->DimensionModel == GAIA_XY_M)
			  result = gaiaAllocGeomCollXYM ();
		      else if (ring->DimensionModel == GAIA_XY_Z_M)
			  result = gaiaAllocGeomCollXYZM ();
		      else
			  result = gaiaAllocGeomColl ();
		      result->Srid = geo->Srid;
		      line = gaiaAddLinestringToGeomColl (result, ring->Points);
		      for (iv = 0; iv < line->Points; iv++)
			{
			    if (ring->DimensionModel == GAIA_XY_Z)
			      {
				  gaiaGetPointXYZ (ring->Coords, iv, &x, &y,
						   &z);
				  gaiaSetPointXYZ (line->Coords, iv, x, y, z);
			      }
			    else if (ring->DimensionModel == GAIA_XY_M)
			      {
				  gaiaGetPointXYM (ring->Coords, iv, &x, &y,
						   &m);
				  gaiaSetPointXYM (line->Coords, iv, x, y, m);
			      }
			    else if (ring->DimensionModel == GAIA_XY_Z_M)
			      {
				  gaiaGetPointXYZM (ring->Coords, iv, &x, &y,
						    &z, &m);
				  gaiaSetPointXYZM (line->Coords, iv, x, y, z,
						    m);
			      }
			    else
			      {
				  gaiaGetPoint (ring->Coords, iv, &x, &y);
				  gaiaSetPoint (line->Coords, iv, x, y);
			      }
			}
		      gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		      gaiaFreeGeomColl (result);
		      sqlite3_result_blob (context, p_result, len, free);
		  }
		else
		    sqlite3_result_null (context);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_NumGeometries (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ NumGeometries(BLOB encoded GEOMETRYCOLLECTION)
/
/ returns the number of elementary geometries for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int cnt = 0;
    gaiaPointPtr point;
    gaiaLinestringPtr line;
    gaiaPolygonPtr polyg;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = geo->FirstPoint;
	  while (point)
	    {
		/* counts how many points are there */
		cnt++;
		point = point->Next;
	    }
	  line = geo->FirstLinestring;
	  while (line)
	    {
		/* counts how many linestrings are there */
		cnt++;
		line = line->Next;
	    }
	  polyg = geo->FirstPolygon;
	  while (polyg)
	    {
		/* counts how many polygons are there */
		cnt++;
		polyg = polyg->Next;
	    }
	  sqlite3_result_int (context, cnt);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_GeometryN (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeometryN(BLOB encoded GEOMETRYCOLLECTION geometry)
/
/ returns the Nth geometry for current GEOMETRYCOLLECTION or MULTIxxxx geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int entity;
    int len;
    int cnt = 0;
    int iv;
    int ib;
    double x;
    double y;
    double z;
    double m;
    gaiaPointPtr point;
    gaiaLinestringPtr line;
    gaiaLinestringPtr line2;
    gaiaPolygonPtr polyg;
    gaiaPolygonPtr polyg2;
    gaiaRingPtr ring_in;
    gaiaRingPtr ring_out;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    entity = sqlite3_value_int (argv[1]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  point = geo->FirstPoint;
	  while (point)
	    {
		/* counts how many points are there */
		cnt++;
		if (cnt == entity)
		  {
		      /* ok, required elementary geometry is this POINT */
		      if (point->DimensionModel == GAIA_XY_Z)
			  result = gaiaAllocGeomCollXYZ ();
		      else if (point->DimensionModel == GAIA_XY_M)
			  result = gaiaAllocGeomCollXYM ();
		      else if (point->DimensionModel == GAIA_XY_Z_M)
			  result = gaiaAllocGeomCollXYZM ();
		      else
			  result = gaiaAllocGeomColl ();
		      result->Srid = geo->Srid;
		      if (point->DimensionModel == GAIA_XY_Z)
			  gaiaAddPointToGeomCollXYZ (result, point->X,
						     point->Y, point->Z);
		      else if (point->DimensionModel == GAIA_XY_M)
			  gaiaAddPointToGeomCollXYM (result, point->X,
						     point->Y, point->M);
		      else if (point->DimensionModel == GAIA_XY_Z_M)
			  gaiaAddPointToGeomCollXYZM (result, point->X,
						      point->Y, point->Z,
						      point->M);
		      else
			  gaiaAddPointToGeomColl (result, point->X, point->Y);
		      goto skip;
		  }
		point = point->Next;
	    }
	  line = geo->FirstLinestring;
	  while (line)
	    {
		/* counts how many linestrings are there */
		cnt++;
		if (cnt == entity)
		  {
		      /* ok, required elementary geometry is this LINESTRING */
		      if (line->DimensionModel == GAIA_XY_Z)
			  result = gaiaAllocGeomCollXYZ ();
		      else if (line->DimensionModel == GAIA_XY_M)
			  result = gaiaAllocGeomCollXYM ();
		      else if (line->DimensionModel == GAIA_XY_Z_M)
			  result = gaiaAllocGeomCollXYZM ();
		      else
			  result = gaiaAllocGeomColl ();
		      result->Srid = geo->Srid;
		      line2 =
			  gaiaAddLinestringToGeomColl (result, line->Points);
		      for (iv = 0; iv < line2->Points; iv++)
			{
			    if (line->DimensionModel == GAIA_XY_Z)
			      {
				  gaiaGetPointXYZ (line->Coords, iv, &x, &y,
						   &z);
				  gaiaSetPointXYZ (line2->Coords, iv, x, y, z);
			      }
			    else if (line->DimensionModel == GAIA_XY_M)
			      {
				  gaiaGetPointXYM (line->Coords, iv, &x, &y,
						   &m);
				  gaiaSetPointXYM (line2->Coords, iv, x, y, m);
			      }
			    else if (line->DimensionModel == GAIA_XY_Z_M)
			      {
				  gaiaGetPointXYZM (line->Coords, iv, &x, &y,
						    &z, &m);
				  gaiaSetPointXYZM (line2->Coords, iv, x, y,
						    z, m);
			      }
			    else
			      {
				  gaiaGetPoint (line->Coords, iv, &x, &y);
				  gaiaSetPoint (line2->Coords, iv, x, y);
			      }
			}
		      goto skip;
		  }
		line = line->Next;
	    }
	  polyg = geo->FirstPolygon;
	  while (polyg)
	    {
		/* counts how many polygons are there */
		cnt++;
		if (cnt == entity)
		  {
		      /* ok, required elementary geometry is this POLYGON */
		      if (polyg->DimensionModel == GAIA_XY_Z)
			  result = gaiaAllocGeomCollXYZ ();
		      else if (polyg->DimensionModel == GAIA_XY_M)
			  result = gaiaAllocGeomCollXYM ();
		      else if (polyg->DimensionModel == GAIA_XY_Z_M)
			  result = gaiaAllocGeomCollXYZM ();
		      else
			  result = gaiaAllocGeomColl ();
		      result->Srid = geo->Srid;
		      ring_in = polyg->Exterior;
		      polyg2 =
			  gaiaAddPolygonToGeomColl (result, ring_in->Points,
						    polyg->NumInteriors);
		      ring_out = polyg2->Exterior;
		      for (iv = 0; iv < ring_out->Points; iv++)
			{
			    /* copying the exterior ring POINTs */
			    gaiaGetPoint (ring_in->Coords, iv, &x, &y);
			    gaiaSetPoint (ring_out->Coords, iv, x, y);
			}
		      for (ib = 0; ib < polyg2->NumInteriors; ib++)
			{
			    /* processing the interior rings */
			    ring_in = polyg->Interiors + ib;
			    ring_out =
				gaiaAddInteriorRing (polyg2, ib,
						     ring_in->Points);
			    for (iv = 0; iv < ring_out->Points; iv++)
			      {
				  if (ring_in->DimensionModel == GAIA_XY_Z)
				    {
					gaiaGetPointXYZ (ring_in->Coords, iv,
							 &x, &y, &z);
					gaiaSetPointXYZ (ring_out->Coords, iv,
							 x, y, z);
				    }
				  else if (ring_in->DimensionModel == GAIA_XY_M)
				    {
					gaiaGetPointXYM (ring_in->Coords, iv,
							 &x, &y, &m);
					gaiaSetPointXYM (ring_out->Coords, iv,
							 x, y, m);
				    }
				  else if (ring_in->DimensionModel ==
					   GAIA_XY_Z_M)
				    {
					gaiaGetPointXYZM (ring_in->Coords, iv,
							  &x, &y, &z, &m);
					gaiaSetPointXYZM (ring_out->Coords,
							  iv, x, y, z, m);
				    }
				  else
				    {
					gaiaGetPoint (ring_in->Coords, iv, &x,
						      &y);
					gaiaSetPoint (ring_out->Coords, iv, x,
						      y);
				    }
			      }
			}
		      goto skip;
		  }
		polyg = polyg->Next;
	    }
	skip:
	  if (result)
	    {
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		gaiaFreeGeomColl (result);
		sqlite3_result_blob (context, p_result, len, free);
	    }
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo);
}

static void
mbrs_eval (sqlite3_context * context, int argc, sqlite3_value ** argv,
	   int request)
{
/* SQL function:
/ MBRsomething(BLOB encoded GEOMETRY-1, BLOB encoded GEOMETRY-2)
/
/ returns:
/ 1 if the required spatial relationship between the two MBRs is TRUE
/ 0 otherwise
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int ret;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobMbr (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobMbr (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  ret = 0;
	  gaiaMbrGeometry (geo1);
	  gaiaMbrGeometry (geo2);
	  switch (request)
	    {
	    case GAIA_MBR_CONTAINS:
		ret = gaiaMbrsContains (geo1, geo2);
		break;
	    case GAIA_MBR_DISJOINT:
		ret = gaiaMbrsDisjoint (geo1, geo2);
		break;
	    case GAIA_MBR_EQUAL:
		ret = gaiaMbrsEqual (geo1, geo2);
		break;
	    case GAIA_MBR_INTERSECTS:
		ret = gaiaMbrsIntersects (geo1, geo2);
		break;
	    case GAIA_MBR_OVERLAPS:
		ret = gaiaMbrsOverlaps (geo1, geo2);
		break;
	    case GAIA_MBR_TOUCHES:
		ret = gaiaMbrsTouches (geo1, geo2);
		break;
	    case GAIA_MBR_WITHIN:
		ret = gaiaMbrsWithin (geo1, geo2);
		break;
	    }
	  if (ret < 0)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

/*
/ the following functions simply readdress the mbr_eval()
/ setting the appropriate request mode
*/

static void
fnct_MbrContains (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_CONTAINS);
}

static void
fnct_MbrDisjoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_DISJOINT);
}

static void
fnct_MbrEqual (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_EQUAL);
}

static void
fnct_MbrIntersects (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_INTERSECTS);
}

static void
fnct_MbrOverlaps (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_OVERLAPS);
}

static void
fnct_MbrTouches (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_TOUCHES);
}

static void
fnct_MbrWithin (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    mbrs_eval (context, argc, argv, GAIA_MBR_WITHIN);
}

static void
fnct_ShiftCoords (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ShiftCoords(BLOBencoded geometry, shiftX, shiftY)
/
/ returns a new geometry that is the original one received, but with shifted coordinates
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    double shift_x;
    double shift_y;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	shift_x = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  shift_x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	shift_y = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  shift_y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaShiftCoords (geo, shift_x, shift_y);
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ScaleCoords (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ScaleCoords(BLOBencoded geometry, scale_factor_x [, scale_factor_y])
/
/ returns a new geometry that is the original one received, but with scaled coordinates
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    double scale_x;
    double scale_y;
    int int_value;
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	scale_x = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  scale_x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (argc == 2)
	scale_y = scale_x;	/* this one is an isotropic scaling request */
    else
      {
	  /* an anisotropic scaling is requested */
	  if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	      scale_y = sqlite3_value_double (argv[2]);
	  else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	    {
		int_value = sqlite3_value_int (argv[2]);
		scale_y = int_value;
	    }
	  else
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaScaleCoords (geo, scale_x, scale_y);
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_RotateCoords (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ RotateCoords(BLOBencoded geometry, angle)
/
/ returns a new geometry that is the original one received, but with rotated coordinates
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    double angle;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	angle = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  angle = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaRotateCoords (geo, angle);
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ReflectCoords (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ReflectCoords(BLOBencoded geometry, x_axis,  y_axis)
/
/ returns a new geometry that is the original one received, but with mirrored coordinates
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    int x_axis;
    int y_axis;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	x_axis = sqlite3_value_int (argv[1]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	y_axis = sqlite3_value_int (argv[2]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaReflectCoords (geo, x_axis, y_axis);
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_SwapCoords (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ SwapCoords(BLOBencoded geometry)
/
/ returns a new geometry that is the original one received, but with swapped x- and y-coordinate
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaSwapCoords (geo);
	  gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
	  if (!p_result)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_blob (context, p_result, len, free);
      }
    gaiaFreeGeomColl (geo);
}

static int
get_ellipse_params (sqlite3 * sqlite, int srid, double *a, double *b,
		    double *rf)
{
/* 
/ retrieves the PROJ +ellps=xx [+a=xx +b=xx] params 
/from SPATIAL_SYS_REF table, if possible 
*/
    char proj4text[2048];
    char *p_proj;
    char *p_ellps;
    char *p_a;
    char *p_b;
    char *p_end;
    proj_params (sqlite, srid, proj4text);
    if (*proj4text == '\0')
	return 0;
/* parsing the proj4text geodesic string */
    p_proj = strstr (proj4text, "+proj=");
    p_ellps = strstr (proj4text, "+ellps=");
    p_a = strstr (proj4text, "+a=");
    p_b = strstr (proj4text, "+b=");
/* checking if +proj=longlat is true */
    if (!p_proj)
	return 0;
    p_end = strchr (p_proj, ' ');
    if (p_end)
	*p_end = '\0';
    if (strcmp (p_proj + 6, "longlat") != 0)
	return 0;
    if (p_ellps)
      {
	  /* trying to retrieve the ellipsoid params by name */
	  p_end = strchr (p_ellps, ' ');
	  if (p_end)
	      *p_end = '\0';
	  if (gaiaEllipseParams (p_ellps + 7, a, b, rf))
	      return 1;
      }
    if (p_a && p_b)
      {
	  /* trying to retrieve the +a=xx and +b=xx args */
	  p_end = strchr (p_a, ' ');
	  if (p_end)
	      *p_end = '\0';
	  p_end = strchr (p_b, ' ');
	  if (p_end)
	      *p_end = '\0';
	  *a = atof (p_a + 3);
	  *b = atof (p_b + 3);
	  *rf = 1.0 / ((*a - *b) / *a);
	  return 1;
      }
    return 0;
}

static void
fnct_FromEWKB (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromEWKB(EWKB encoded geometry)
/
/ returns the current geometry by parsing Geos/PostGis EWKB encoded string 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaFromEWKB (text);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_ToEWKB (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsEWKB(BLOB encoded geometry)
/
/ returns a text string corresponding to Geos/PostGIS EWKB notation 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
      {
	  sqlite3_result_null (context);
	  return;
      }
    else
      {
	  gaiaOutBufferInitialize (&out_buf);
	  gaiaToEWKB (&out_buf, geo);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_ToEWKT (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ AsEWKT(BLOB encoded geometry)
/
/ returns the corresponding PostGIS EWKT encoded value
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    gaiaOutBuffer out_buf;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    gaiaOutBufferInitialize (&out_buf);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  gaiaToEWKT (&out_buf, geo);
	  if (out_buf.Error || out_buf.Buffer == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		len = out_buf.WriteOffset;
		sqlite3_result_text (context, out_buf.Buffer, len, free);
		out_buf.Buffer = NULL;
	    }
      }
    gaiaFreeGeomColl (geo);
    gaiaOutBufferReset (&out_buf);
}

static void
fnct_FromEWKT (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromEWKT(EWKT encoded geometry)
/
/ returns the current geometry by parsing EWKT  (PostGIS) encoded string 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseEWKT (text);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_FromGeoJSON (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromGeoJSON(GeoJSON encoded geometry)
/
/ returns the current geometry by parsing GeoJSON encoded string 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseGeoJSON (text);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_FromKml (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromKml(KML encoded geometry)
/
/ returns the current geometry by parsing KML encoded string 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseKml (text);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_FromGml (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromGml(GML encoded geometry)
/
/ returns the current geometry by parsing GML encoded string 
/ or NULL if any error is encountered
*/
    int len;
    unsigned char *p_result = NULL;
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseGml (text, sqlite);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geo, &p_result, &len);
    gaiaFreeGeomColl (geo);
    sqlite3_result_blob (context, p_result, len, free);
}

static void
fnct_LinesFromRings (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ LinesFromRings(BLOBencoded geometry, BOOL multi_linestring)
/
/ returns a new geometry [LINESTRING or MULTILINESTRING] representing 
/ the linearization for current (MULTI)POLYGON geometry
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr geom_new = NULL;
    int len;
    int multi_linestring;
    unsigned char *p_result = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (argc == 2)
      {
	  if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	      multi_linestring = sqlite3_value_int (argv[1]);
      }
    geom_new = gaiaLinearize (geo, multi_linestring);
    if (!geom_new)
	goto invalid;
    gaiaFreeGeomColl (geo);
    gaiaToSpatiaLiteBlobWkb (geom_new, &p_result, &len);
    gaiaFreeGeomColl (geom_new);
    sqlite3_result_blob (context, p_result, len, free);
    return;
  invalid:
    if (geo)
	gaiaFreeGeomColl (geo);
    sqlite3_result_null (context);
}

#ifndef OMIT_GEOS		/* including GEOS */

static void
fnct_BuildArea (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ BuildArea(BLOBencoded geometry)
/
/ Assuming that Geometry represents a set of sparse Linestrings,
/ this function will attempt to reassemble a single Polygon
/ (or a set of Polygons)
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaPolygonize (geo, 0);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}


static void
fnct_Polygonize_step (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ Polygonize(BLOBencoded geom)
/
/ aggregate function - STEP
/
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr *p;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geom)
	return;
    p = sqlite3_aggregate_context (context, sizeof (gaiaGeomCollPtr));
    if (!(*p))
      {
	  /* this is the first row */
	  *p = geom;
      }
    else
      {
	  /* subsequent rows */
	  result = gaiaMergeGeometries (*p, geom);
	  gaiaFreeGeomColl (*p);
	  *p = result;
	  gaiaFreeGeomColl (geom);
      }
}

static void
fnct_Polygonize_final (sqlite3_context * context)
{
/* SQL function:
/ Polygonize(BLOBencoded geom)
/
/ aggregate function - FINAL
/
*/
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr geom;
    gaiaGeomCollPtr *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    result = *p;
    if (!result)
	sqlite3_result_null (context);
    else
      {
	  geom = gaiaPolygonize (result, 0);
	  if (geom == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		geom->Srid = result->Srid;
		gaiaToSpatiaLiteBlobWkb (geom, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (geom);
	    }
	  gaiaFreeGeomColl (result);
      }
}

#endif /* end including GEOS */

static void
fnct_DissolveSegments (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ DissolveSegments(BLOBencoded geometry)
/
/ Dissolves any LINESTRING or RING into elementary segments
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaDissolveSegments (geo);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_DissolvePoints (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ DissolvePoints(BLOBencoded geometry)
/
/ Dissolves any LINESTRING or RING into elementary Vertices
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaDissolvePoints (geo);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_CollectionExtract (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ CollectionExtract(BLOBencoded geometry, Integer type)
/
/ Extracts from a GEOMETRYCOLLECTION any item of the required TYPE
/ 1=Point - 2=Linestring - 3=Polygon
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    int type;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	type = sqlite3_value_int (argv[1]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (type == 1 || type == 2 || type == 3)
	;
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  switch (type)
	    {
	    case 1:
		result = gaiaExtractPointsFromGeomColl (geo);
		break;
	    case 2:
		result = gaiaExtractLinestringsFromGeomColl (geo);
		break;
	    case 3:
		result = gaiaExtractPolygonsFromGeomColl (geo);
		break;
	    };
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

#ifndef OMIT_PROJ		/* including PROJ.4 */

static void
fnct_Transform (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Transform(BLOBencoded geometry, srid)
/
/ returns a new geometry that is the original one received, but with the new SRID [no coordinates translation is applied]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    int srid_from;
    int srid_to;
    char proj_from[2048];
    char proj_to[2048];
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
	srid_to = sqlite3_value_int (argv[1]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  *proj_from = '\0';
	  *proj_to = '\0';
	  srid_from = geo->Srid;
	  proj_params (sqlite, srid_from, proj_from);
	  proj_params (sqlite, srid_to, proj_to);
	  if (*proj_to == '\0' || *proj_from == '\0')
	    {
		gaiaFreeGeomColl (geo);
		sqlite3_result_null (context);
		return;
	    }
	  result = gaiaTransform (geo, proj_from, proj_to);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = srid_to;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

#endif /* end including PROJ.4 */

#ifndef OMIT_GEOS		/* including GEOS */

static void
fnct_Boundary (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Boundary(BLOB encoded geometry)
/
/ returns the combinatorial boundary for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr boundary;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (gaiaIsEmpty (geo))
	      sqlite3_result_null (context);
	  else
	    {
		boundary = gaiaBoundary (geo);
		if (!boundary)
		    sqlite3_result_null (context);
		else
		  {
		      gaiaToSpatiaLiteBlobWkb (boundary, &p_result, &len);
		      gaiaFreeGeomColl (boundary);
		      sqlite3_result_blob (context, p_result, len, free);
		  }
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_IsClosed (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ IsClosed(BLOB encoded LINESTRING or MULTILINESTRING geometry)
/
/ returns:
/ 1 if this LINESTRING is closed [or if this is a MULTILINESTRING and every LINESTRINGs are closed] 
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_int (context, -1);
    else
      {
	  line = simpleLinestring (geo);
	  if (!line < 0)
	      sqlite3_result_int (context, -1);
	  else
	      sqlite3_result_int (context, gaiaIsClosed (line));
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_IsSimple (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ IsSimple(BLOB encoded GEOMETRY)
/
/ returns:
/ 1 if this GEOMETRY is simple
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int ret;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaIsSimple (geo);
	  if (ret < 0)
	      sqlite3_result_int (context, -1);
	  else
	      sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_IsRing (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ IsRing(BLOB encoded LINESTRING geometry)
/
/ returns:
/ 1 if this LINESTRING is a valid RING
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int ret;
    gaiaGeomCollPtr geo = NULL;
    gaiaLinestringPtr line;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_int (context, -1);
    else
      {
	  line = simpleLinestring (geo);
	  if (!line < 0)
	      sqlite3_result_int (context, -1);
	  else
	    {
		ret = gaiaIsRing (line);
		sqlite3_result_int (context, ret);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_IsValid (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ IsValid(BLOB encoded GEOMETRY)
/
/ returns:
/ 1 if this GEOMETRY is a valid one
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int ret;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaIsValid (geo);
	  if (ret < 0)
	      sqlite3_result_int (context, -1);
	  else
	      sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Length (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GLength(BLOB encoded GEOMETRYCOLLECTION)
/
/ returns  the total length for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double length = 0.0;
    int ret;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  ret = gaiaGeomCollLength (geo, &length);
	  if (!ret)
	      sqlite3_result_null (context);
	  sqlite3_result_double (context, length);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Area (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Area(BLOB encoded GEOMETRYCOLLECTION)
/
/ returns the total area for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double area = 0.0;
    int ret;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  ret = gaiaGeomCollArea (geo, &area);
	  if (!ret)
	      sqlite3_result_null (context);
	  sqlite3_result_double (context, area);
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Centroid (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Centroid(BLOBencoded POLYGON or MULTIPOLYGON geometry)
/
/ returns a POINT representing the centroid for current POLYGON / MULTIPOLYGON geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    int ret;
    double x;
    double y;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (gaiaIsEmpty (geo))
	      sqlite3_result_null (context);
	  else
	    {
		ret = gaiaGeomCollCentroid (geo, &x, &y);
		if (!ret)
		    sqlite3_result_null (context);
		else
		  {
		      result = gaiaAllocGeomColl ();
		      result->Srid = geo->Srid;
		      gaiaAddPointToGeomColl (result, x, y);
		      gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		      gaiaFreeGeomColl (result);
		      sqlite3_result_blob (context, p_result, len, free);
		  }
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_PointOnSurface (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ PointOnSurface(BLOBencoded POLYGON or MULTIPOLYGON geometry)
/
/ returns a POINT guaranteed to lie on the Surface
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    double x;
    double y;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (!gaiaGetPointOnSurface (geo, &x, &y))
	      sqlite3_result_null (context);
	  else
	    {
		result = gaiaAllocGeomColl ();
		gaiaAddPointToGeomColl (result, x, y);
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		gaiaFreeGeomColl (result);
		sqlite3_result_blob (context, p_result, len, free);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Simplify (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Simplify(BLOBencoded geometry, tolerance)
/
/ returns a new geometry that is a caricature of the original one received, but simplified using the Douglas-Peuker algorihtm
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    int int_value;
    double tolerance;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	tolerance = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  tolerance = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeomCollSimplify (geo, tolerance);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_SimplifyPreserveTopology (sqlite3_context * context, int argc,
			       sqlite3_value ** argv)
{
/* SQL function:
/ SimplifyPreserveTopology(BLOBencoded geometry, tolerance)
/
/ returns a new geometry that is a caricature of the original one received, but simplified using the Douglas-Peuker algorihtm [preserving topology]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    int int_value;
    double tolerance;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	tolerance = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  tolerance = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeomCollSimplifyPreserveTopology (geo, tolerance);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ConvexHull (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ConvexHull(BLOBencoded geometry)
/
/ returns a new geometry representing the CONVEX HULL for current geometry
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int len;
    unsigned char *p_result = NULL;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaConvexHull (geo);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Buffer (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Buffer(BLOBencoded geometry, radius)
/
/ returns a new geometry representing the BUFFER for current geometry
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    double radius;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	radius = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  radius = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeomCollBuffer (geo, radius, 30);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_Intersection (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Intersection(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns a new geometry representing the INTERSECTION of both geometries
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeometryIntersection (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else if (gaiaIsEmpty (result))
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Union_step (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Union(BLOBencoded geom)
/
/ aggregate function - STEP
/
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr *p;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geom = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geom)
	return;
    p = sqlite3_aggregate_context (context, sizeof (gaiaGeomCollPtr));
    if (!(*p))
      {
	  /* this is the first row */
	  *p = geom;
      }
    else
      {
	  /* subsequent rows */
	  result = gaiaGeometryUnion (*p, geom);
	  gaiaFreeGeomColl (*p);
	  *p = result;
	  gaiaFreeGeomColl (geom);
      }
}

static void
fnct_Union_final (sqlite3_context * context)
{
/* SQL function:
/ Union(BLOBencoded geom)
/
/ aggregate function - FINAL
/
*/
    gaiaGeomCollPtr result;
    gaiaGeomCollPtr *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    result = *p;
    if (!result)
	sqlite3_result_null (context);
    else if (gaiaIsEmpty (result))
      {
	  gaiaFreeGeomColl (result);
	  sqlite3_result_null (context);
      }
    else
      {
	  /* builds the BLOB geometry to be returned */
	  int len;
	  unsigned char *p_result = NULL;
	  gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
	  sqlite3_result_blob (context, p_result, len, free);
	  gaiaFreeGeomColl (result);
      }
}

static void
fnct_Union (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Union(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns a new geometry representing the UNION of both geometries
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeometryUnion (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else if (gaiaIsEmpty (result))
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Difference (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Difference(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns a new geometry representing the DIFFERENCE of both geometries
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeometryDifference (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else if (gaiaIsEmpty (result))
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_SymDifference (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ SymDifference(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns a new geometry representing the SYMMETRIC DIFFERENCE of both geometries
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaGeometrySymDifference (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else if (gaiaIsEmpty (result))
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Equals (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Equals(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries are "spatially equal"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollEquals (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Intersects (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Intersects(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries do "spatially intersects"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollIntersects (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Disjoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Disjoint(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries are "spatially disjoint"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollDisjoint (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Overlaps (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Overlaps(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries do "spatially overlaps"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollOverlaps (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Crosses (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Crosses(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries do "spatially crosses"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollCrosses (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Touches (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Touches(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if the two geometries do "spatially touches"
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollTouches (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Within (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Within(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if GEOM-1 is completely contained within GEOM-2
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollWithin (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Contains (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Contains(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if GEOM-1 completely contains GEOM-2
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollContains (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Relate (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Relate(BLOBencoded geom1, BLOBencoded geom2, string pattern)
/
/ returns:
/ 1 if GEOM-1 and GEOM-2 have a spatial relationship as specified by the patternMatrix 
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    const unsigned char *pattern;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[2]) != SQLITE_TEXT)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    pattern = sqlite3_value_text (argv[2]);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollRelate (geo1, geo2, (char *) pattern);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Distance (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Distance(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns the distance between GEOM-1 and GEOM-2
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    double dist;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  ret = gaiaGeomCollDistance (geo1, geo2, &dist);
	  if (!ret)
	      sqlite3_result_null (context);
	  sqlite3_result_double (context, dist);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_PtDistWithin (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ PtDistWithin(BLOBencoded geom1, BLOBencoded geom2, double dist 
/ [, boolen use_spheroid])
/
/ returns TRUE if the distance between GEOM-1 and GEOM-2
/ is less or equal to dist
/
/ - if both geom1 and geom2 are in the 4326 (WGS84) SRID,
/   (and does actually contains a single POINT each one)
/   dist is assumed to be measured in Meters
/ - in this case the optional arg use_spheroid is
/   checked to determine if geodesic distance has to be
/   computed on the sphere (quickest) or on the spheroid 
/   default: use_spheroid = FALSE
/ 
/ in any other case the "plain" distance is evaluated
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaPointPtr pt;
    gaiaLinestringPtr ln;
    gaiaPolygonPtr pg;
    double ref_dist;
    int use_spheroid = 0;
    double x0;
    double y0;
    double x1;
    double y1;
    int pt0 = 0;
    int ln0 = 0;
    int pg0 = 0;
    int pt1 = 0;
    int ln1 = 0;
    int pg1 = 0;
    double dist;
    double a;
    double b;
    double rf;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER
	|| sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	;
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (argc == 4)
      {
	  /* optional use_spheroid arg */
	  if (sqlite3_value_type (argv[3]) != SQLITE_INTEGER)
	    {
		sqlite3_result_null (context);
		return;
	    }
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int dst = sqlite3_value_int (argv[2]);
	  ref_dist = dst;
      }
    else
	ref_dist = sqlite3_value_double (argv[2]);
    if (argc == 4)
	use_spheroid = sqlite3_value_int (argv[3]);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  if (geo1->Srid == 4326 && geo2->Srid == 4326)
	    {
		/* checking for single points */
		pt = geo1->FirstPoint;
		while (pt)
		  {
		      x0 = pt->X;
		      y0 = pt->Y;
		      pt0++;
		      pt = pt->Next;
		  }
		ln = geo1->FirstLinestring;
		while (ln)
		  {
		      ln0++;
		      ln = ln->Next;
		  }
		pg = geo1->FirstPolygon;
		while (pg)
		  {
		      pg0++;
		      pg = pg->Next;
		  }
		pt = geo2->FirstPoint;
		while (pt)
		  {
		      x1 = pt->X;
		      y1 = pt->Y;
		      pt1++;
		      pt = pt->Next;
		  }
		ln = geo2->FirstLinestring;
		while (ln)
		  {
		      ln1++;
		      ln = ln->Next;
		  }
		pg = geo2->FirstPolygon;
		while (pg)
		  {
		      pg1++;
		      pg = pg->Next;
		  }
		if (pt0 == 1 && pt1 == 1 && ln0 == 0 && ln1 == 0 && pg0 == 0
		    && pg1 == 0)
		  {
		      /* using geodesic distance */
		      a = 6378137.0;
		      rf = 298.257223563;
		      b = (a * (1.0 - (1.0 / rf)));
		      if (use_spheroid)
			{
			    dist =
				gaiaGeodesicDistance (a, b, rf, y0, x0, y1, x1);
			    if (dist <= ref_dist)
				sqlite3_result_int (context, 1);
			    else
				sqlite3_result_int (context, 0);
			}
		      else
			{
			    dist =
				gaiaGreatCircleDistance (a, b, y0, x0, y1, x1);
			    if (dist <= ref_dist)
				sqlite3_result_int (context, 1);
			    else
				sqlite3_result_int (context, 0);
			}
		  }
		goto stop;
	    }
/* defaulting to flat distance */
	  ret = gaiaGeomCollDistance (geo1, geo2, &dist);
	  if (!ret)
	      sqlite3_result_null (context);
	  if (dist <= ref_dist)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
      }
  stop:
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
geos_error (const char *fmt, ...)
{
/* reporting some GEOS error */
    va_list ap;
    char msg[2048];
    va_start (ap, fmt);
    vsprintf (msg, fmt, ap);
    va_end (ap);
    fprintf (stderr, "GEOS error: %s\n", msg);
    gaiaSetGeosErrorMsg (msg);
}


static void
geos_warning (const char *fmt, ...)
{
/* reporting some GEOS warning */
    va_list ap;
    char msg[2048];
    va_start (ap, fmt);
    vsprintf (msg, fmt, ap);
    va_end (ap);
    fprintf (stderr, "GEOS warning: %s\n", msg);
    gaiaSetGeosWarningMsg (msg);
}

static void
fnct_aux_polygonize (sqlite3_context * context, gaiaGeomCollPtr geom_org,
		     int force_multipolygon, int allow_multipolygon)
{
/* a  common function performing any kind of polygonization op */
    gaiaGeomCollPtr geom_new = NULL;
    int len;
    unsigned char *p_result = NULL;
    gaiaPolygonPtr pg;
    int pgs = 0;
    if (!geom_org)
	goto invalid;
    geom_new = gaiaPolygonize (geom_org, force_multipolygon);
    if (!geom_new)
	goto invalid;
    gaiaFreeGeomColl (geom_org);
    pg = geom_new->FirstPolygon;
    while (pg)
      {
	  pgs++;
	  pg = pg->Next;
      }
    if (pgs > 1 && allow_multipolygon == 0)
      {
	  /* invalid: a POLYGON is expected !!! */
	  gaiaFreeGeomColl (geom_new);
	  sqlite3_result_null (context);
	  return;
      }
    gaiaToSpatiaLiteBlobWkb (geom_new, &p_result, &len);
    gaiaFreeGeomColl (geom_new);
    sqlite3_result_blob (context, p_result, len, free);
    return;
  invalid:
    if (geom_org)
	gaiaFreeGeomColl (geom_org);
    sqlite3_result_null (context);
}

/*
/ the following functions performs initial argument checking, 
/ and then readdressing the request to fnct_aux_polygonize()
/ for actual processing
*/

static void
fnct_BdPolyFromText1 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BdPolyFromText(WKT encoded MULTILINESTRING)
/
/ returns the current geometry [POLYGON] by parsing a WKT encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, -1);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = -1;
    fnct_aux_polygonize (context, geo, 0, 0);
    return;
}

static void
fnct_BdPolyFromText2 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BdPolyFromText(WKT encoded MULTILINESTRING, SRID)
/
/ returns the current geometry [POLYGON] by parsing a WKT encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, -1);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    fnct_aux_polygonize (context, geo, 0, 0);
    return;
}

static void
fnct_BdMPolyFromText1 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ BdMPolyFromText(WKT encoded MULTILINESTRING)
/
/ returns the current geometry [MULTIPOLYGON] by parsing a WKT encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, -1);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = -1;
    fnct_aux_polygonize (context, geo, 1, 1);
    return;
}

static void
fnct_BdMPolyFromText2 (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ BdMPolyFromText(WKT encoded MULTILINESTRING, SRID)
/
/ returns the current geometry [MULTIPOLYGON] by parsing a WKT encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    const unsigned char *text;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_TEXT)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    text = sqlite3_value_text (argv[0]);
    geo = gaiaParseWkt (text, -1);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    fnct_aux_polygonize (context, geo, 1, 1);
    return;
}

static void
fnct_BdPolyFromWKB1 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ BdPolyFromWKB(WKB encoded MULTILINESTRING)
/
/ returns the current geometry [POLYGON] by parsing a WKB encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    int n_bytes;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, -1))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = -1;
    fnct_aux_polygonize (context, geo, 0, 0);
    return;
}

static void
fnct_BdPolyFromWKB2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ BdPolyFromWKB(WKB encoded MULTILINESTRING)
/
/ returns the current geometry [POLYGON] by parsing a WKB encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    int n_bytes;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, -1))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    fnct_aux_polygonize (context, geo, 0, 0);
    return;
}

static void
fnct_BdMPolyFromWKB1 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BdMPolyFromWKB(WKB encoded MULTILINESTRING)
/
/ returns the current geometry [MULTIPOLYGON] by parsing a WKB encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    int n_bytes;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, -1))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = -1;
    fnct_aux_polygonize (context, geo, 1, 1);
    return;
}

static void
fnct_BdMPolyFromWKB2 (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ BdMPolyFromWKB(WKB encoded MULTILINESTRING)
/
/ returns the current geometry [MULTIPOLYGON] by parsing a WKB encoded MULTILINESTRING 
/ or NULL if any error is encountered
/
*/
    int n_bytes;
    const unsigned char *wkb;
    gaiaGeomCollPtr geo = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_INTEGER)
      {
	  sqlite3_result_null (context);
	  return;
      }
    wkb = sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (!check_wkb (wkb, n_bytes, -1))
	return;
    geo = gaiaFromWkb (wkb, n_bytes);
    if (geo == NULL)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (geo->DeclaredType != GAIA_MULTILINESTRING)
      {
	  sqlite3_result_null (context);
	  return;
      }
    geo->Srid = sqlite3_value_int (argv[1]);
    fnct_aux_polygonize (context, geo, 1, 1);
    return;
}

#ifdef GEOS_ADVANCED		/* GEOS advanced and experimental features */

static void
fnct_OffsetCurve (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ OffsetCurve(BLOBencoded geometry, radius, left-or-right-side)
/
/ returns a new geometry representing the OFFSET-CURVE for current geometry
/ [a LINESTRING is expected]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    double radius;
    int int_value;
    int left_right;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	radius = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  radius = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	left_right = sqlite3_value_int (argv[2]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaOffsetCurve (geo, radius, 16, left_right);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_SingleSidedBuffer (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ SingleSidedBuffer(BLOBencoded geometry, radius, left-or-right-side)
/
/ returns a new geometry representing the SingleSided BUFFER 
/ for current geometry [a LINESTRING is expected]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    double radius;
    int int_value;
    int left_right;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	radius = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  radius = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
	left_right = sqlite3_value_int (argv[2]);
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaSingleSidedBuffer (geo, radius, 16, left_right);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_HausdorffDistance (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ HausdorffDistance(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns the discrete Hausdorff distance between GEOM-1 and GEOM-2
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    double dist;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_null (context);
    else
      {
	  ret = gaiaHausdorffDistance (geo1, geo2, &dist);
	  if (!ret)
	      sqlite3_result_null (context);
	  sqlite3_result_double (context, dist);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_SharedPaths (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ SharedPaths(BLOBencoded geometry1, BLOBencoded geometry2)
/
/ returns a new geometry representing common (shared) Edges
/ [two LINESTRINGs/MULTILINESTRINGs are expected]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo1 == NULL || geo2 == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaSharedPaths (geo1, geo2);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo1->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Covers (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Covers(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if GEOM-1 "spatially covers" GEOM-2
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollCovers (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_CoveredBy (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ CoveredBy(BLOBencoded geom1, BLOBencoded geom2)
/
/ returns:
/ 1 if GEOM-1 is "spatially covered by" GEOM-2
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    int ret;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo1 || !geo2)
	sqlite3_result_int (context, -1);
    else
      {
	  ret = gaiaGeomCollCoveredBy (geo1, geo2);
	  sqlite3_result_int (context, ret);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_LineInterpolatePoint (sqlite3_context * context, int argc,
			   sqlite3_value ** argv)
{
/* SQL function:
/ LineInterpolatePoint(BLOBencoded geometry1, double fraction)
/
/ returns a new geometry representing a point interpolated along a line
/ [a LINESTRING is expected / fraction ranging from 0.0 to 1.0]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int int_value;
    double fraction;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	fraction = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  fraction = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaLineInterpolatePoint (geo, fraction);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_LineLocatePoint (sqlite3_context * context, int argc,
		      sqlite3_value ** argv)
{
/* SQL function:
/ LineLocatePoint(BLOBencoded geometry1, BLOBencoded geometry2)
/
/ return a number (between 0.0 and 1.0) representing the location 
/ of the closest point on LineString to the given Point, as a fraction 
/ of total 2d line length
/
/ - geom1 is expected to represent some LINESTRING
/ - geom2 is expected to represent some POINT
*/
    unsigned char *p_blob;
    int n_bytes;
    double fraction;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo1 == NULL || geo2 == NULL)
	sqlite3_result_null (context);
    else
      {
	  fraction = gaiaLineLocatePoint (geo1, geo2);
	  if (fraction >= 0.0 && fraction <= 1.0)
	      sqlite3_result_double (context, fraction);
	  else
	      sqlite3_result_null (context);
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_LineSubstring (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ LineSubstring(BLOBencoded geometry1, double start_fraction, double end_fraction)
/
/ Return a Linestring being a substring of the input one starting and ending at 
/ the given fractions of total 2d length [fractions ranging from 0.0 to 1.0]
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int int_value;
    double fraction1;
    double fraction2;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	fraction1 = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  fraction1 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	fraction2 = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  fraction2 = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaLineSubstring (geo, fraction1, fraction2);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

static void
fnct_ClosestPoint (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ClosestPoint(BLOBencoded geometry1, BLOBencoded geometry2)
/
/ Returns the Point on geom1 that is closest to geom2
/ NULL is returned for invalid arguments (or if distance is ZERO)
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo1 == NULL || geo2 == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaShortestLine (geo1, geo2);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else if (result->FirstLinestring == NULL)
	    {
		gaiaFreeGeomColl (result);
		sqlite3_result_null (context);
	    }
	  else
	    {
		/* builds the BLOB geometry to be returned */
		double x;
		double y;
		double z;
		double m;
		int len;
		unsigned char *p_result = NULL;
		gaiaGeomCollPtr pt = NULL;
		gaiaLinestringPtr ln = result->FirstLinestring;
		if (ln->DimensionModel == GAIA_XY_Z)
		    pt = gaiaAllocGeomCollXYZ ();
		else if (ln->DimensionModel == GAIA_XY_M)
		    pt = gaiaAllocGeomCollXYM ();
		else if (ln->DimensionModel == GAIA_XY_Z_M)
		    pt = gaiaAllocGeomCollXYZM ();
		else
		    pt = gaiaAllocGeomColl ();
		if (ln->DimensionModel == GAIA_XY_Z)
		  {
		      gaiaGetPointXYZ (ln->Coords, 0, &x, &y, &z);
		      gaiaAddPointToGeomCollXYZ (pt, x, y, z);
		  }
		else if (ln->DimensionModel == GAIA_XY_M)
		  {
		      gaiaGetPointXYM (ln->Coords, 0, &x, &y, &m);
		      gaiaAddPointToGeomCollXYM (pt, x, y, m);
		  }
		else if (ln->DimensionModel == GAIA_XY_Z_M)
		  {
		      gaiaGetPointXYZM (ln->Coords, 0, &x, &y, &z, &m);
		      gaiaAddPointToGeomCollXYZM (pt, x, y, z, m);
		  }
		else
		  {
		      gaiaGetPoint (ln->Coords, 0, &x, &y);
		      gaiaAddPointToGeomColl (pt, x, y);
		  }
		pt->Srid = geo1->Srid;
		gaiaToSpatiaLiteBlobWkb (pt, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
		gaiaFreeGeomColl (pt);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_ShortestLine (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ShortestLine(BLOBencoded geometry1, BLOBencoded geometry2)
/
/ Returns the shortest line between two geometries
/ NULL is returned for invalid arguments (or if distance is ZERO)
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo1 == NULL || geo2 == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaShortestLine (geo1, geo2);
	  sqlite3_result_null (context);
	  if (!result)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo1->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_Snap (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ Snap(BLOBencoded geometry1, BLOBencoded geometry2, double tolerance)
/
/ Returns a new Geometry corresponding to geom1 snapped to geom2
/ and using the given tolerance
/ NULL is returned for invalid arguments (or if distance is ZERO)
*/
    unsigned char *p_blob;
    int n_bytes;
    int int_value;
    double tolerance;
    gaiaGeomCollPtr geo1 = NULL;
    gaiaGeomCollPtr geo2 = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[2]) == SQLITE_FLOAT)
	tolerance = sqlite3_value_double (argv[2]);
    else if (sqlite3_value_type (argv[2]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[2]);
	  tolerance = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo1 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    p_blob = (unsigned char *) sqlite3_value_blob (argv[1]);
    n_bytes = sqlite3_value_bytes (argv[1]);
    geo2 = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo1 == NULL || geo2 == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaSnap (geo1, geo2, tolerance);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo1->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo1);
    gaiaFreeGeomColl (geo2);
}

static void
fnct_LineMerge (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ LineMerge(BLOBencoded geometry)
/
/ Assuming that Geometry represents a set of sparse Linestrings,
/ this function will attempt to reassemble a single line
/ (or a set of lines)
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaLineMerge (geo);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}


static void
fnct_UnaryUnion (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ UnaryUnion(BLOBencoded geometry)
/
/ exactly like Union, but using a single Collection
/ NULL is returned for invalid arguments
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geo = NULL;
    gaiaGeomCollPtr result;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (geo == NULL)
	sqlite3_result_null (context);
    else
      {
	  result = gaiaUnaryUnion (geo);
	  if (result == NULL)
	      sqlite3_result_null (context);
	  else
	    {
		/* builds the BLOB geometry to be returned */
		int len;
		unsigned char *p_result = NULL;
		result->Srid = geo->Srid;
		gaiaToSpatiaLiteBlobWkb (result, &p_result, &len);
		sqlite3_result_blob (context, p_result, len, free);
		gaiaFreeGeomColl (result);
	    }
      }
    gaiaFreeGeomColl (geo);
}

#endif /* end GEOS advanced and experimental features */

#endif /* end including GEOS */

#ifndef OMIT_MATHSQL		/* supporting SQL math functions */

static void
fnct_math_acos (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ acos(double X)
/
/ Returns the arc cosine of X, that is, the value whose cosine is X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    errno = 0;
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = acos (sqlite3_value_double (argv[0]));
	  if (errno == EDOM)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = acos (x);
	  if (errno == EDOM)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_asin (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ asin(double X)
/
/ Returns the arc sine of X, that is, the value whose sine is X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    errno = 0;
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = asin (sqlite3_value_double (argv[0]));
	  if (errno == EDOM)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = asin (x);
	  if (errno == EDOM)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_atan (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ atan(double X)
/
/ Returns the arc tangent of X, that is, the value whose tangent is X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = atan (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = atan (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_ceil (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ ceil(double X)
/
/ Returns the smallest integer value not less than X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = ceil (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = ceil (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_cos (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ cos(double X)
/
/ Returns the cosine of X, where X is given in radians
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = cos (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = cos (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_cot (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ cot(double X)
/
/ Returns the cotangent of X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    double tang;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    tang = tan (x);
    if (tang == 0.0)
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = 1.0 / tang;
    sqlite3_result_double (context, x);
}

static void
fnct_math_degrees (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ degrees(double X)
/
/ Returns the argument X, converted from radians to degrees
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = x * 57.29577951308232;
    sqlite3_result_double (context, x);
}

static void
fnct_math_exp (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ exp(double X)
/
/ Returns the value of e (the base of natural logarithms) raised to the power of X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = exp (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = exp (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_floor (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ floor(double X)
/
/ Returns the largest integer value not greater than X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = floor (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = floor (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_logn (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ log(double X)
/
/ Returns the natural logarithm of X; that is, the base-e logarithm of X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    errno = 0;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = log (sqlite3_value_double (argv[0]));
	  if (errno == EDOM || errno == ERANGE)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = log (x);
	  if (errno == EDOM || errno == ERANGE)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_logn2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ log(double B, double X)
/
/ Returns the logarithm of X to the base B
/ or NULL if any error is encountered
*/
    int int_value;
    double x = 0.0;
    double b = 1.0;
    double log1;
    double log2;
    errno = 0;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	b = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  b = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (x <= 0.0 || b <= 1.0)
      {
	  sqlite3_result_null (context);
	  return;
      }
    log1 = log (x);
    if (errno == EDOM || errno == ERANGE)
      {
	  sqlite3_result_null (context);
	  return;
      }
    log2 = log (b);
    if (errno == EDOM || errno == ERANGE)
      {
	  sqlite3_result_null (context);
	  return;
      }
    sqlite3_result_double (context, log1 / log2);
}

static void
fnct_math_log_2 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ log2(double X)
/
/ Returns the base-2 logarithm of X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    double log1;
    double log2;
    errno = 0;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    log1 = log (x);
    if (errno == EDOM || errno == ERANGE)
      {
	  sqlite3_result_null (context);
	  return;
      }
    log2 = log (2.0);
    sqlite3_result_double (context, log1 / log2);
}

static void
fnct_math_log_10 (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ log10(double X)
/
/ Returns the base-10 logarithm of X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    double log1;
    double log2;
    GAIA_UNUSED ();
    errno = 0;
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    log1 = log (x);
    if (errno == EDOM || errno == ERANGE)
      {
	  sqlite3_result_null (context);
	  return;
      }
    log2 = log (10.0);
    sqlite3_result_double (context, log1 / log2);
}

static void
fnct_math_pi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ pi(void)
/
/ Returns the value of (pi)
*/
    GAIA_UNUSED ();
    sqlite3_result_double (context, 3.14159265358979323846);
}

static void
fnct_math_pow (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ pow(double X, double Y)
/
/ Returns the value of X raised to the power of Y.
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    double y;
    double p;
    GAIA_UNUSED ();
    errno = 0;
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (sqlite3_value_type (argv[1]) == SQLITE_FLOAT)
	y = sqlite3_value_double (argv[1]);
    else if (sqlite3_value_type (argv[1]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[1]);
	  y = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    p = pow (x, y);
    if (errno == EDOM)
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, p);
}

static void
fnct_math_stddev_step (sqlite3_context * context, int argc,
		       sqlite3_value ** argv)
{
/* SQL function:
/ stddev_pop(double X)
/ stddev_samp(double X)
/ var_pop(double X)
/ var_samp(double X)
/
/ aggregate function - STEP
/
*/
    struct stddev_str *p;
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
	return;
    p = sqlite3_aggregate_context (context, sizeof (struct stddev_str));
    if (!(p->cleaned))
      {
	  p->cleaned = 1;
	  p->mean = x;
	  p->quot = 0.0;
	  p->count = 0.0;
      }
    p->count += 1.0;
    p->quot =
	p->quot +
	(((p->count - 1.0) * ((x - p->mean) * (x - p->mean))) / p->count);
    p->mean = p->mean + ((x - p->mean) / p->count);
}

static void
fnct_math_stddev_pop_final (sqlite3_context * context)
{
/* SQL function:
/ stddev_pop(double X)
/ aggregate function -  FINAL
/
*/
    double x;
    struct stddev_str *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = sqrt (p->quot / p->count);
    sqlite3_result_double (context, x);
}

static void
fnct_math_stddev_samp_final (sqlite3_context * context)
{
/* SQL function:
/ stddev_samp(double X)
/ aggregate function -  FINAL
/
*/
    double x;
    struct stddev_str *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = sqrt (p->quot / (p->count - 1.0));
    sqlite3_result_double (context, x);
}

static void
fnct_math_var_pop_final (sqlite3_context * context)
{
/* SQL function:
/ var_pop(double X)
/ aggregate function -  FINAL
/
*/
    double x;
    struct stddev_str *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = p->quot / p->count;
    sqlite3_result_double (context, x);
}

static void
fnct_math_var_samp_final (sqlite3_context * context)
{
/* SQL function:
/ var_samp(double X)
/ aggregate function -  FINAL
/
*/
    double x;
    struct stddev_str *p = sqlite3_aggregate_context (context, 0);
    if (!p)
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = p->quot / (p->count - 1.0);
    sqlite3_result_double (context, x);
}

static void
fnct_math_radians (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ radians(double X)
/
/ Returns the argument X, converted from degrees to radians
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    x = x * .0174532925199432958;
    sqlite3_result_double (context, x);
}


static void
fnct_math_round (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ round(double X)
/
/ Returns the the nearest integer, but round halfway cases away from zero
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = math_round (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = math_round (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_sign (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ sign(double X)
/
/ Returns the sign of the argument as -1, 0, or 1, depending on whether X is negative, zero, or positive
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	x = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (x > 0.0)
	sqlite3_result_double (context, 1.0);
    else if (x < 0.0)
	sqlite3_result_double (context, -1.0);
    else
	sqlite3_result_double (context, 0.0);
}

static void
fnct_math_sin (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ sin(double X)
/
/ Returns the sine of X, where X is given in radians
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = sin (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = sin (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_sqrt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ sqrt(double X)
/
/ Returns the square root of a non-negative number X
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    errno = 0;
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = sqrt (sqlite3_value_double (argv[0]));
	  if (errno)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = sqrt (x);
	  if (errno == EDOM)
	      sqlite3_result_null (context);
	  else
	      sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

static void
fnct_math_tan (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ tan(double X)
/
/ Returns the tangent of X, where X is given in radians
/ or NULL if any error is encountered
*/
    int int_value;
    double x;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
      {
	  x = tan (sqlite3_value_double (argv[0]));
	  sqlite3_result_double (context, x);
      }
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  x = int_value;
	  x = tan (x);
	  sqlite3_result_double (context, x);
      }
    else
	sqlite3_result_null (context);
}

#endif /* end supporting SQL math functions */

static void
fnct_GeomFromExifGpsBlob (sqlite3_context * context, int argc,
			  sqlite3_value ** argv)
{
/* SQL function:
/ GeomFromExifGpsBlob(BLOB encoded image)
/
/ returns:
/ a POINT geometry
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    gaiaGeomCollPtr geom;
    unsigned char *geoblob;
    int geosize;
    double longitude;
    double latitude;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    if (gaiaGetGpsCoords (p_blob, n_bytes, &longitude, &latitude))
      {
	  geom = gaiaAllocGeomColl ();
	  geom->Srid = 4326;
	  gaiaAddPointToGeomColl (geom, longitude, latitude);
	  gaiaToSpatiaLiteBlobWkb (geom, &geoblob, &geosize);
	  gaiaFreeGeomColl (geom);
	  sqlite3_result_blob (context, geoblob, geosize, free);
      }
    else
	sqlite3_result_null (context);
}

static void
blob_guess (sqlite3_context * context, int argc, sqlite3_value ** argv,
	    int request)
{
/* SQL function:
/ IsGifBlob(BLOB encoded image)
/ IsPngBlob, IsJpegBlob, IsExifBlob, IsExifGpsBlob, IsTiffBlob,
/ IsZipBlob, IsPdfBlob,IsGeometryBlob
/
/ returns:
/ 1 if the required BLOB_TYPE is TRUE
/ 0 otherwise
/ or -1 if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    int blob_type;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_int (context, -1);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    blob_type = gaiaGuessBlobType (p_blob, n_bytes);
    if (request == GAIA_GEOMETRY_BLOB)
      {
	  if (blob_type == GAIA_GEOMETRY_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_ZIP_BLOB)
      {
	  if (blob_type == GAIA_ZIP_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_PDF_BLOB)
      {
	  if (blob_type == GAIA_PDF_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_TIFF_BLOB)
      {
	  if (blob_type == GAIA_TIFF_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_GIF_BLOB)
      {
	  if (blob_type == GAIA_GIF_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_PNG_BLOB)
      {
	  if (blob_type == GAIA_PNG_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_JPEG_BLOB)
      {
	  if (blob_type == GAIA_JPEG_BLOB || blob_type == GAIA_EXIF_BLOB
	      || blob_type == GAIA_EXIF_GPS_BLOB)
	      sqlite3_result_int (context, 1);
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_EXIF_BLOB)
      {
	  if (blob_type == GAIA_EXIF_BLOB || blob_type == GAIA_EXIF_GPS_BLOB)
	    {
		sqlite3_result_int (context, 1);
	    }
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    if (request == GAIA_EXIF_GPS_BLOB)
      {
	  if (blob_type == GAIA_EXIF_GPS_BLOB)
	    {
		sqlite3_result_int (context, 1);
	    }
	  else
	      sqlite3_result_int (context, 0);
	  return;
      }
    sqlite3_result_int (context, -1);
}

/*
/ the following functions simply readdress the blob_guess()
/ setting the appropriate request mode
*/

static void
fnct_IsGeometryBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_GEOMETRY_BLOB);
}

static void
fnct_IsZipBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_ZIP_BLOB);
}

static void
fnct_IsPdfBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_PDF_BLOB);
}

static void
fnct_IsTiffBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_TIFF_BLOB);
}

static void
fnct_IsGifBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_GIF_BLOB);
}

static void
fnct_IsPngBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_PNG_BLOB);
}

static void
fnct_IsJpegBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_JPEG_BLOB);
}

static void
fnct_IsExifBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_EXIF_BLOB);
}

static void
fnct_IsExifGpsBlob (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    blob_guess (context, argc, argv, GAIA_EXIF_GPS_BLOB);
}

static void
fnct_GeodesicLength (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
/* SQL function:
/ GeodesicLength(BLOB encoded GEOMETRYCOLLECTION)
/
/ returns  the total Geodesic length for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double l;
    double length = 0.0;
    double a;
    double b;
    double rf;
    gaiaGeomCollPtr geo = NULL;
    gaiaLinestringPtr line;
    gaiaPolygonPtr polyg;
    gaiaRingPtr ring;
    int ib;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (get_ellipse_params (sqlite, geo->Srid, &a, &b, &rf))
	    {
		line = geo->FirstLinestring;
		while (line)
		  {
		      /* Linestrings */
		      l = gaiaGeodesicTotalLength (a, b, rf,
						   line->DimensionModel,
						   line->Coords, line->Points);
		      if (l < 0.0)
			{
			    length = -1.0;
			    break;
			}
		      length += l;
		      line = line->Next;
		  }
		if (length >= 0)
		  {
		      /* Polygons */
		      polyg = geo->FirstPolygon;
		      while (polyg)
			{
			    /* exterior Ring */
			    ring = polyg->Exterior;
			    l = gaiaGeodesicTotalLength (a, b, rf,
							 ring->DimensionModel,
							 ring->Coords,
							 ring->Points);
			    if (l < 0.0)
			      {
				  length = -1.0;
				  break;
			      }
			    length += l;
			    for (ib = 0; ib < polyg->NumInteriors; ib++)
			      {
				  /* interior Rings */
				  ring = polyg->Interiors + ib;
				  l = gaiaGeodesicTotalLength (a, b, rf,
							       ring->
							       DimensionModel,
							       ring->Coords,
							       ring->Points);
				  if (l < 0.0)
				    {
					length = -1.0;
					break;
				    }
				  length += l;
			      }
			    if (length < 0.0)
				break;
			    polyg = polyg->Next;
			}
		  }
		if (length < 0.0)
		    sqlite3_result_null (context);
		else
		    sqlite3_result_double (context, length);
	    }
	  else
	      sqlite3_result_null (context);
	  gaiaFreeGeomColl (geo);
      }
}

static void
fnct_GreatCircleLength (sqlite3_context * context, int argc,
			sqlite3_value ** argv)
{
/* SQL function:
/ GreatCircleLength(BLOB encoded GEOMETRYCOLLECTION)
/
/ returns  the total Great Circle length for current geometry 
/ or NULL if any error is encountered
*/
    unsigned char *p_blob;
    int n_bytes;
    double length = 0.0;
    double a;
    double b;
    double rf;
    gaiaGeomCollPtr geo = NULL;
    gaiaLinestringPtr line;
    gaiaPolygonPtr polyg;
    gaiaRingPtr ring;
    int ib;
    sqlite3 *sqlite = sqlite3_context_db_handle (context);
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) != SQLITE_BLOB)
      {
	  sqlite3_result_null (context);
	  return;
      }
    p_blob = (unsigned char *) sqlite3_value_blob (argv[0]);
    n_bytes = sqlite3_value_bytes (argv[0]);
    geo = gaiaFromSpatiaLiteBlobWkb (p_blob, n_bytes);
    if (!geo)
	sqlite3_result_null (context);
    else
      {
	  if (get_ellipse_params (sqlite, geo->Srid, &a, &b, &rf))
	    {
		line = geo->FirstLinestring;
		while (line)
		  {
		      /* Linestrings */
		      length +=
			  gaiaGreatCircleTotalLength (a, b,
						      line->DimensionModel,
						      line->Coords,
						      line->Points);
		      line = line->Next;
		  }
		if (length >= 0)
		  {
		      /* Polygons */
		      polyg = geo->FirstPolygon;
		      while (polyg)
			{
			    /* exterior Ring */
			    ring = polyg->Exterior;
			    length +=
				gaiaGreatCircleTotalLength (a, b,
							    ring->
							    DimensionModel,
							    ring->Coords,
							    ring->Points);
			    for (ib = 0; ib < polyg->NumInteriors; ib++)
			      {
				  /* interior Rings */
				  ring = polyg->Interiors + ib;
				  length +=
				      gaiaGreatCircleTotalLength (a, b,
								  ring->
								  DimensionModel,
								  ring->Coords,
								  ring->Points);
			      }
			    polyg = polyg->Next;
			}
		  }
		sqlite3_result_double (context, length);
	    }
	  else
	      sqlite3_result_null (context);
	  gaiaFreeGeomColl (geo);
      }
}

static void
convertUnit (sqlite3_context * context, int argc, sqlite3_value ** argv,
	     int unit_from, int unit_to)
{
/* SQL functions:
/ CvtToKm(), CvtToDm(), CvtToCm(), CvtToMm(), CvtToKmi(), CvtToIn(), CvtToFt(),
/ CvtToYd(), CvtToMi(), CvtToFath(), CvtToCh(), CvtToLink(), CvtToUsIn(), 
/ CvtToUsFt(), CvtToUsYd(), CvtToUsCh(), CvtToUsMi(), CvtToIndFt(), 
/ CvtToIndYd(), CvtToIndCh(), 
/ CvtFromKm(), CvtFromDm(), CvtFromCm(), CvtFromMm(), CvtFromKmi(), 
/ CvtFromIn(), CvtFromFt(), CvtFromYd(), CvtFromMi(), CvtFromFath(), 
/ CvtFromCh(), CvtFromLink(), CvtFromUsIn(), CvtFromUsFt(), CvtFromUsYd(), 
/ CvtFromUsCh(), CvtFromUsMi(), CvtFromIndFt(), CvtFromIndYd(), 
/ CvtFromIndCh()
/
/ converts a Length from one unit to a different one
/ or NULL if any error is encountered
*/
    double cvt;
    double value;
    int int_value;
    GAIA_UNUSED ();
    if (sqlite3_value_type (argv[0]) == SQLITE_FLOAT)
	value = sqlite3_value_double (argv[0]);
    else if (sqlite3_value_type (argv[0]) == SQLITE_INTEGER)
      {
	  int_value = sqlite3_value_int (argv[0]);
	  value = int_value;
      }
    else
      {
	  sqlite3_result_null (context);
	  return;
      }
    if (!gaiaConvertLength (value, unit_from, unit_to, &cvt))
	sqlite3_result_null (context);
    else
	sqlite3_result_double (context, cvt);
}

static void
fnct_cvtToKm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_KM);
}

static void
fnct_cvtToDm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_DM);
}

static void
fnct_cvtToCm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_CM);
}

static void
fnct_cvtToMm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_MM);
}

static void
fnct_cvtToKmi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_KMI);
}

static void
fnct_cvtToIn (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_IN);
}

static void
fnct_cvtToFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_FT);
}

static void
fnct_cvtToYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_YD);
}

static void
fnct_cvtToMi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_MI);
}

static void
fnct_cvtToFath (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_FATH);
}

static void
fnct_cvtToCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_CH);
}

static void
fnct_cvtToLink (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_LINK);
}

static void
fnct_cvtToUsIn (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_US_IN);
}

static void
fnct_cvtToUsFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_US_FT);
}

static void
fnct_cvtToUsYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_US_YD);
}

static void
fnct_cvtToUsCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_US_CH);
}

static void
fnct_cvtToUsMi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_US_MI);
}

static void
fnct_cvtToIndFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_IND_FT);
}

static void
fnct_cvtToIndYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_IND_YD);
}

static void
fnct_cvtToIndCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_M, GAIA_IND_CH);
}

static void
fnct_cvtFromKm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_KM, GAIA_M);
}

static void
fnct_cvtFromDm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_DM, GAIA_M);
}

static void
fnct_cvtFromCm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_CM, GAIA_M);
}

static void
fnct_cvtFromMm (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_MM, GAIA_M);
}

static void
fnct_cvtFromKmi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_KMI, GAIA_M);
}

static void
fnct_cvtFromIn (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_IN, GAIA_M);
}

static void
fnct_cvtFromFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_FT, GAIA_M);
}

static void
fnct_cvtFromYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_YD, GAIA_M);
}

static void
fnct_cvtFromMi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_MI, GAIA_M);
}

static void
fnct_cvtFromFath (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_FATH, GAIA_M);
}

static void
fnct_cvtFromCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_CH, GAIA_M);
}

static void
fnct_cvtFromLink (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_LINK, GAIA_M);
}

static void
fnct_cvtFromUsIn (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_US_IN, GAIA_M);
}

static void
fnct_cvtFromUsFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_US_FT, GAIA_M);
}

static void
fnct_cvtFromUsYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_US_YD, GAIA_M);
}

static void
fnct_cvtFromUsCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_US_CH, GAIA_M);
}

static void
fnct_cvtFromUsMi (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_US_MI, GAIA_M);
}

static void
fnct_cvtFromIndFt (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_IND_FT, GAIA_M);
}

static void
fnct_cvtFromIndYd (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_IND_YD, GAIA_M);
}

static void
fnct_cvtFromIndCh (sqlite3_context * context, int argc, sqlite3_value ** argv)
{
    convertUnit (context, argc, argv, GAIA_IND_CH, GAIA_M);
}

static void
register_spatialite_sql_functions (sqlite3 * db)
{
    sqlite3_create_function (db, "spatialite_version", 0, SQLITE_ANY, 0,
			     fnct_spatialite_version, 0, 0);
    sqlite3_create_function (db, "proj4_version", 0, SQLITE_ANY, 0,
			     fnct_proj4_version, 0, 0);
    sqlite3_create_function (db, "geos_version", 0, SQLITE_ANY, 0,
			     fnct_geos_version, 0, 0);
    sqlite3_create_function (db, "GeometryConstraints", 3, SQLITE_ANY, 0,
			     fnct_GeometryConstraints, 0, 0);
    sqlite3_create_function (db, "GeometryConstraints", 4, SQLITE_ANY, 0,
			     fnct_GeometryConstraints, 0, 0);
    sqlite3_create_function (db, "RTreeAlign", 3, SQLITE_ANY, 0,
			     fnct_RTreeAlign, 0, 0);
    sqlite3_create_function (db, "CheckSpatialMetaData", 0, SQLITE_ANY, 0,
			     fnct_CheckSpatialMetaData, 0, 0);
    sqlite3_create_function (db, "AutoFDOStart", 0, SQLITE_ANY, 0,
			     fnct_AutoFDOStart, 0, 0);
    sqlite3_create_function (db, "AutoFDOStop", 0, SQLITE_ANY, 0,
			     fnct_AutoFDOStop, 0, 0);
    sqlite3_create_function (db, "InitFDOSpatialMetaData", 0, SQLITE_ANY, 0,
			     fnct_InitFDOSpatialMetaData, 0, 0);
    sqlite3_create_function (db, "AddFDOGeometryColumn", 6, SQLITE_ANY, 0,
			     fnct_AddFDOGeometryColumn, 0, 0);
    sqlite3_create_function (db, "RecoverFDOGeometryColumn", 6, SQLITE_ANY, 0,
			     fnct_RecoverFDOGeometryColumn, 0, 0);
    sqlite3_create_function (db, "DiscardFDOGeometryColumn", 2, SQLITE_ANY, 0,
			     fnct_DiscardFDOGeometryColumn, 0, 0);
    sqlite3_create_function (db, "InitSpatialMetaData", 0, SQLITE_ANY, 0,
			     fnct_InitSpatialMetaData, 0, 0);
    sqlite3_create_function (db, "AddGeometryColumn", 5, SQLITE_ANY, 0,
			     fnct_AddGeometryColumn, 0, 0);
    sqlite3_create_function (db, "AddGeometryColumn", 6, SQLITE_ANY, 0,
			     fnct_AddGeometryColumn, 0, 0);
    sqlite3_create_function (db, "RecoverGeometryColumn", 5, SQLITE_ANY, 0,
			     fnct_RecoverGeometryColumn, 0, 0);
    sqlite3_create_function (db, "DiscardGeometryColumn", 2, SQLITE_ANY, 0,
			     fnct_DiscardGeometryColumn, 0, 0);
    sqlite3_create_function (db, "RecoverSpatialIndex", 0, SQLITE_ANY, 0,
			     fnct_RecoverSpatialIndex, 0, 0);
    sqlite3_create_function (db, "RecoverSpatialIndex", 1, SQLITE_ANY, 0,
			     fnct_RecoverSpatialIndex, 0, 0);
    sqlite3_create_function (db, "RecoverSpatialIndex", 2, SQLITE_ANY, 0,
			     fnct_RecoverSpatialIndex, 0, 0);
    sqlite3_create_function (db, "RecoverSpatialIndex", 3, SQLITE_ANY, 0,
			     fnct_RecoverSpatialIndex, 0, 0);
    sqlite3_create_function (db, "CheckSpatialIndex", 0, SQLITE_ANY, 0,
			     fnct_CheckSpatialIndex, 0, 0);
    sqlite3_create_function (db, "CheckSpatialIndex", 2, SQLITE_ANY, 0,
			     fnct_CheckSpatialIndex, 0, 0);
    sqlite3_create_function (db, "CreateSpatialIndex", 2, SQLITE_ANY, 0,
			     fnct_CreateSpatialIndex, 0, 0);
    sqlite3_create_function (db, "CreateMbrCache", 2, SQLITE_ANY, 0,
			     fnct_CreateMbrCache, 0, 0);
    sqlite3_create_function (db, "DisableSpatialIndex", 2, SQLITE_ANY, 0,
			     fnct_DisableSpatialIndex, 0, 0);
    sqlite3_create_function (db, "RebuildGeometryTriggers", 2, SQLITE_ANY, 0,
			     fnct_RebuildGeometryTriggers, 0, 0);
    sqlite3_create_function (db, "CreateTopologyTables", 2, SQLITE_ANY, 0,
			     fnct_CreateTopologyTables, 0, 0);
    sqlite3_create_function (db, "CreateTopologyTables", 3, SQLITE_ANY, 0,
			     fnct_CreateTopologyTables, 0, 0);
    sqlite3_create_function (db, "UpdateLayerStatistics", 0, SQLITE_ANY, 0,
			     fnct_UpdateLayerStatistics, 0, 0);
    sqlite3_create_function (db, "UpdateLayerStatistics", 1, SQLITE_ANY, 0,
			     fnct_UpdateLayerStatistics, 0, 0);
    sqlite3_create_function (db, "UpdateLayerStatistics", 2, SQLITE_ANY, 0,
			     fnct_UpdateLayerStatistics, 0, 0);
    sqlite3_create_function (db, "AsText", 1, SQLITE_ANY, 0, fnct_AsText, 0, 0);
    sqlite3_create_function (db, "ST_AsText", 1, SQLITE_ANY, 0, fnct_AsText,
			     0, 0);
    sqlite3_create_function (db, "AsWkt", 1, SQLITE_ANY, 0, fnct_AsWkt, 0, 0);
    sqlite3_create_function (db, "AsWkt", 2, SQLITE_ANY, 0, fnct_AsWkt, 0, 0);
    sqlite3_create_function (db, "AsSvg", 1, SQLITE_ANY, 0, fnct_AsSvg1, 0, 0);
    sqlite3_create_function (db, "AsSvg", 2, SQLITE_ANY, 0, fnct_AsSvg2, 0, 0);
    sqlite3_create_function (db, "AsSvg", 3, SQLITE_ANY, 0, fnct_AsSvg3, 0, 0);

#ifndef OMIT_PROJ		/* PROJ.4 is strictly required to support KML */
    sqlite3_create_function (db, "AsKml", 1, SQLITE_ANY, 0, fnct_AsKml, 0, 0);
    sqlite3_create_function (db, "AsKml", 2, SQLITE_ANY, 0, fnct_AsKml, 0, 0);
    sqlite3_create_function (db, "AsKml", 3, SQLITE_ANY, 0, fnct_AsKml, 0, 0);
    sqlite3_create_function (db, "AsKml", 4, SQLITE_ANY, 0, fnct_AsKml, 0, 0);
#endif /* end including PROJ.4 */

    sqlite3_create_function (db, "AsGml", 1, SQLITE_ANY, 0, fnct_AsGml, 0, 0);
    sqlite3_create_function (db, "AsGml", 2, SQLITE_ANY, 0, fnct_AsGml, 0, 0);
    sqlite3_create_function (db, "AsGml", 3, SQLITE_ANY, 0, fnct_AsGml, 0, 0);
    sqlite3_create_function (db, "GeomFromGml", 1, SQLITE_ANY, 0,
			     fnct_FromGml, 0, 0);
    sqlite3_create_function (db, "AsGeoJSON", 1, SQLITE_ANY, 0, fnct_AsGeoJSON,
			     0, 0);
    sqlite3_create_function (db, "AsGeoJSON", 2, SQLITE_ANY, 0, fnct_AsGeoJSON,
			     0, 0);
    sqlite3_create_function (db, "AsGeoJSON", 3, SQLITE_ANY, 0, fnct_AsGeoJSON,
			     0, 0);
    sqlite3_create_function (db, "GeomFromGeoJSON", 1, SQLITE_ANY, 0,
			     fnct_FromGeoJSON, 0, 0);
    sqlite3_create_function (db, "GeomFromKml", 1, SQLITE_ANY, 0,
			     fnct_FromKml, 0, 0);
    sqlite3_create_function (db, "AsFGF", 2, SQLITE_ANY, 0, fnct_AsFGF, 0, 0);
    sqlite3_create_function (db, "GeomFromEWKB", 1, SQLITE_ANY, 0,
			     fnct_FromEWKB, 0, 0);
    sqlite3_create_function (db, "AsEWKB", 1, SQLITE_ANY, 0, fnct_ToEWKB, 0, 0);
    sqlite3_create_function (db, "AsEWKT", 1, SQLITE_ANY, 0, fnct_ToEWKT, 0, 0);
    sqlite3_create_function (db, "GeomFromEWKT", 1, SQLITE_ANY, 0,
			     fnct_FromEWKT, 0, 0);
    sqlite3_create_function (db, "AsBinary", 1, SQLITE_ANY, 0, fnct_AsBinary, 0,
			     0);
    sqlite3_create_function (db, "ST_AsBinary", 1, SQLITE_ANY, 0, fnct_AsBinary,
			     0, 0);
    sqlite3_create_function (db, "GeomFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomFromText1, 0, 0);
    sqlite3_create_function (db, "GeomFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomFromText2, 0, 0);
    sqlite3_create_function (db, "GeometryFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomFromText1, 0, 0);
    sqlite3_create_function (db, "GeometryFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomFromText2, 0, 0);
    sqlite3_create_function (db, "GeomCollFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromText1, 0, 0);
    sqlite3_create_function (db, "GeomCollFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromText2, 0, 0);
    sqlite3_create_function (db, "GeometryCollectionFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromText1, 0, 0);
    sqlite3_create_function (db, "GeometryCollectionFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromText2, 0, 0);
    sqlite3_create_function (db, "PointFromText", 1, SQLITE_ANY, 0,
			     fnct_PointFromText1, 0, 0);
    sqlite3_create_function (db, "PointFromText", 2, SQLITE_ANY, 0,
			     fnct_PointFromText2, 0, 0);
    sqlite3_create_function (db, "LineFromText", 1, SQLITE_ANY, 0,
			     fnct_LineFromText1, 0, 0);
    sqlite3_create_function (db, "LineFromText", 2, SQLITE_ANY, 0,
			     fnct_LineFromText2, 0, 0);
    sqlite3_create_function (db, "LineStringFromText", 1, SQLITE_ANY, 0,
			     fnct_LineFromText1, 0, 0);
    sqlite3_create_function (db, "LineStringFromText", 2, SQLITE_ANY, 0,
			     fnct_LineFromText2, 0, 0);
    sqlite3_create_function (db, "PolyFromText", 1, SQLITE_ANY, 0,
			     fnct_PolyFromText1, 0, 0);
    sqlite3_create_function (db, "PolyFromText", 2, SQLITE_ANY, 0,
			     fnct_PolyFromText2, 0, 0);
    sqlite3_create_function (db, "PolygonFromText", 1, SQLITE_ANY, 0,
			     fnct_PolyFromText1, 0, 0);
    sqlite3_create_function (db, "PolygonFromText", 2, SQLITE_ANY, 0,
			     fnct_PolyFromText2, 0, 0);
    sqlite3_create_function (db, "MPointFromText", 1, SQLITE_ANY, 0,
			     fnct_MPointFromText1, 0, 0);
    sqlite3_create_function (db, "MPointFromText", 2, SQLITE_ANY, 0,
			     fnct_MPointFromText2, 0, 0);
    sqlite3_create_function (db, "MultiPointFromText", 1, SQLITE_ANY, 0,
			     fnct_MPointFromText1, 0, 0);
    sqlite3_create_function (db, "MultiPointFromText", 2, SQLITE_ANY, 0,
			     fnct_MPointFromText2, 0, 0);
    sqlite3_create_function (db, "MLineFromText", 1, SQLITE_ANY, 0,
			     fnct_MLineFromText1, 0, 0);
    sqlite3_create_function (db, "MLineFromText", 2, SQLITE_ANY, 0,
			     fnct_MLineFromText2, 0, 0);
    sqlite3_create_function (db, "MultiLineStringFromText", 1, SQLITE_ANY, 0,
			     fnct_MLineFromText1, 0, 0);
    sqlite3_create_function (db, "MultiLineStringFromText", 2, SQLITE_ANY, 0,
			     fnct_MLineFromText2, 0, 0);
    sqlite3_create_function (db, "MPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromText1, 0, 0);
    sqlite3_create_function (db, "MPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromText2, 0, 0);
    sqlite3_create_function (db, "MultiPolygonFromText", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromText1, 0, 0);
    sqlite3_create_function (db, "MultiPolygonFromText", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromText2, 0, 0);
    sqlite3_create_function (db, "GeomFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomFromWkb1, 0, 0);
    sqlite3_create_function (db, "GeomFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomFromWkb2, 0, 0);
    sqlite3_create_function (db, "GeometryFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomFromWkb1, 0, 0);
    sqlite3_create_function (db, "GeometryFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomFromWkb2, 0, 0);
    sqlite3_create_function (db, "GeomCollFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb1, 0, 0);
    sqlite3_create_function (db, "GeomCollFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb2, 0, 0);
    sqlite3_create_function (db, "GeometryCollectionFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb1, 0, 0);
    sqlite3_create_function (db, "GeometryCollectionFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb2, 0, 0);
    sqlite3_create_function (db, "PointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PointFromWkb1, 0, 0);
    sqlite3_create_function (db, "PointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PointFromWkb2, 0, 0);
    sqlite3_create_function (db, "LineFromWKB", 1, SQLITE_ANY, 0,
			     fnct_LineFromWkb1, 0, 0);
    sqlite3_create_function (db, "LineFromWKB", 2, SQLITE_ANY, 0,
			     fnct_LineFromWkb2, 0, 0);
    sqlite3_create_function (db, "LineStringFromWKB", 1, SQLITE_ANY, 0,
			     fnct_LineFromWkb1, 0, 0);
    sqlite3_create_function (db, "LineStringFromWKB", 2, SQLITE_ANY, 0,
			     fnct_LineFromWkb2, 0, 0);
    sqlite3_create_function (db, "PolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "PolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "PolygonFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "PolygonFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "MPointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPointFromWkb1, 0, 0);
    sqlite3_create_function (db, "MPointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPointFromWkb2, 0, 0);
    sqlite3_create_function (db, "MultiPointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPointFromWkb1, 0, 0);
    sqlite3_create_function (db, "MultiPointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPointFromWkb2, 0, 0);
    sqlite3_create_function (db, "MLineFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MLineFromWkb1, 0, 0);
    sqlite3_create_function (db, "MLineFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MLineFromWkb2, 0, 0);
    sqlite3_create_function (db, "MultiLineStringFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MLineFromWkb1, 0, 0);
    sqlite3_create_function (db, "MultiLineStringFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MLineFromWkb2, 0, 0);
    sqlite3_create_function (db, "MPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "MPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "MultiPolygonFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "MultiPolygonFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_GeomFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomFromText1, 0, 0);
    sqlite3_create_function (db, "ST_GeomFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomFromText2, 0, 0);
    sqlite3_create_function (db, "ST_GeometryFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomFromText1, 0, 0);
    sqlite3_create_function (db, "ST_GeometryFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomFromText2, 0, 0);
    sqlite3_create_function (db, "ST_GeomCollFromText", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromText1, 0, 0);
    sqlite3_create_function (db, "ST_GeomCollFromText", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromText2, 0, 0);
    sqlite3_create_function (db, "ST_GeometryCollectionFromText", 1, SQLITE_ANY,
			     0, fnct_GeomCollFromText1, 0, 0);
    sqlite3_create_function (db, "ST_GeometryCollectionFromText", 2, SQLITE_ANY,
			     0, fnct_GeomCollFromText2, 0, 0);
    sqlite3_create_function (db, "ST_PointFromText", 1, SQLITE_ANY, 0,
			     fnct_PointFromText1, 0, 0);
    sqlite3_create_function (db, "ST_PointFromText", 2, SQLITE_ANY, 0,
			     fnct_PointFromText2, 0, 0);
    sqlite3_create_function (db, "ST_LineFromText", 1, SQLITE_ANY, 0,
			     fnct_LineFromText1, 0, 0);
    sqlite3_create_function (db, "ST_LineFromText", 2, SQLITE_ANY, 0,
			     fnct_LineFromText2, 0, 0);
    sqlite3_create_function (db, "ST_LineStringFromText", 1, SQLITE_ANY, 0,
			     fnct_LineFromText1, 0, 0);
    sqlite3_create_function (db, "ST_LineStringFromText", 2, SQLITE_ANY, 0,
			     fnct_LineFromText2, 0, 0);
    sqlite3_create_function (db, "ST_PolyFromText", 1, SQLITE_ANY, 0,
			     fnct_PolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_PolyFromText", 2, SQLITE_ANY, 0,
			     fnct_PolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_PolygonFromText", 1, SQLITE_ANY, 0,
			     fnct_PolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_PolygonFromText", 2, SQLITE_ANY, 0,
			     fnct_PolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MPointFromText", 1, SQLITE_ANY, 0,
			     fnct_MPointFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MPointFromText", 2, SQLITE_ANY, 0,
			     fnct_MPointFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MultiPointFromText", 1, SQLITE_ANY, 0,
			     fnct_MPointFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MultiPointFromText", 2, SQLITE_ANY, 0,
			     fnct_MPointFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MLineFromText", 1, SQLITE_ANY, 0,
			     fnct_MLineFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MLineFromText", 2, SQLITE_ANY, 0,
			     fnct_MLineFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MultiLineStringFromText", 1, SQLITE_ANY, 0,
			     fnct_MLineFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MultiLineStringFromText", 2, SQLITE_ANY, 0,
			     fnct_MLineFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_MultiPolygonFromText", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_MultiPolygonFromText", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_GeomFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_GeomFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_GeometryFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_GeometryFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_GeomCollFromWKB", 1, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_GeomCollFromWKB", 2, SQLITE_ANY, 0,
			     fnct_GeomCollFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_GeometryCollectionFromWKB", 1, SQLITE_ANY,
			     0, fnct_GeomCollFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_GeometryCollectionFromWKB", 2, SQLITE_ANY,
			     0, fnct_GeomCollFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_PointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PointFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_PointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PointFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_LineFromWKB", 1, SQLITE_ANY, 0,
			     fnct_LineFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_LineFromWKB", 2, SQLITE_ANY, 0,
			     fnct_LineFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_LineStringFromWKB", 1, SQLITE_ANY, 0,
			     fnct_LineFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_LineStringFromWKB", 2, SQLITE_ANY, 0,
			     fnct_LineFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_PolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_PolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_PolygonFromWKB", 1, SQLITE_ANY, 0,
			     fnct_PolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_PolygonFromWKB", 2, SQLITE_ANY, 0,
			     fnct_PolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MPointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPointFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MPointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPointFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MultiPointFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPointFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MultiPointFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPointFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MLineFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MLineFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MLineFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MLineFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MultiLineStringFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MLineFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MultiLineStringFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MLineFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "ST_MultiPolygonFromWKB", 1, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb1, 0, 0);
    sqlite3_create_function (db, "ST_MultiPolygonFromWKB", 2, SQLITE_ANY, 0,
			     fnct_MPolyFromWkb2, 0, 0);
    sqlite3_create_function (db, "GeomFromFGF", 1, SQLITE_ANY, 0,
			     fnct_GeometryFromFGF1, 0, 0);
    sqlite3_create_function (db, "GeomFromFGF", 2, SQLITE_ANY, 0,
			     fnct_GeometryFromFGF2, 0, 0);
    sqlite3_create_function (db, "CompressGeometry", 1, SQLITE_ANY, 0,
			     fnct_CompressGeometry, 0, 0);
    sqlite3_create_function (db, "UncompressGeometry", 1, SQLITE_ANY, 0,
			     fnct_UncompressGeometry, 0, 0);
    sqlite3_create_function (db, "SanitizeGeometry", 1, SQLITE_ANY, 0,
			     fnct_SanitizeGeometry, 0, 0);
    sqlite3_create_function (db, "CastToPoint", 1, SQLITE_ANY, 0,
			     fnct_CastToPoint, 0, 0);
    sqlite3_create_function (db, "CastToLinestring", 1, SQLITE_ANY, 0,
			     fnct_CastToLinestring, 0, 0);
    sqlite3_create_function (db, "CastToPolygon", 1, SQLITE_ANY, 0,
			     fnct_CastToPolygon, 0, 0);
    sqlite3_create_function (db, "CastToMultiPoint", 1, SQLITE_ANY, 0,
			     fnct_CastToMultiPoint, 0, 0);
    sqlite3_create_function (db, "CastToMultiLinestring", 1, SQLITE_ANY, 0,
			     fnct_CastToMultiLinestring, 0, 0);
    sqlite3_create_function (db, "CastToMultiPolygon", 1, SQLITE_ANY, 0,
			     fnct_CastToMultiPolygon, 0, 0);
    sqlite3_create_function (db, "CastToGeometryCollection", 1, SQLITE_ANY, 0,
			     fnct_CastToGeometryCollection, 0, 0);
    sqlite3_create_function (db, "CastToMulti", 1, SQLITE_ANY, 0,
			     fnct_CastToMulti, 0, 0);
    sqlite3_create_function (db, "ST_Multi", 1, SQLITE_ANY, 0, fnct_CastToMulti,
			     0, 0);
    sqlite3_create_function (db, "CastToSingle", 1, SQLITE_ANY, 0,
			     fnct_CastToSingle, 0, 0);
    sqlite3_create_function (db, "CastToXY", 1, SQLITE_ANY, 0, fnct_CastToXY, 0,
			     0);
    sqlite3_create_function (db, "CastToXYZ", 1, SQLITE_ANY, 0, fnct_CastToXYZ,
			     0, 0);
    sqlite3_create_function (db, "CastToXYM", 1, SQLITE_ANY, 0, fnct_CastToXYM,
			     0, 0);
    sqlite3_create_function (db, "CastToXYZM", 1, SQLITE_ANY, 0,
			     fnct_CastToXYZM, 0, 0);
    sqlite3_create_function (db, "ExtractMultiPoint", 1, SQLITE_ANY, 0,
			     fnct_ExtractMultiPoint, 0, 0);
    sqlite3_create_function (db, "ExtractMultiLinestring", 1, SQLITE_ANY, 0,
			     fnct_ExtractMultiLinestring, 0, 0);
    sqlite3_create_function (db, "ExtractMultiPolygon", 1, SQLITE_ANY, 0,
			     fnct_ExtractMultiPolygon, 0, 0);
    sqlite3_create_function (db, "Dimension", 1, SQLITE_ANY, 0, fnct_Dimension,
			     0, 0);
    sqlite3_create_function (db, "ST_Dimension", 1, SQLITE_ANY, 0,
			     fnct_Dimension, 0, 0);
    sqlite3_create_function (db, "CoordDimension", 1, SQLITE_ANY, 0,
			     fnct_CoordDimension, 0, 0);
    sqlite3_create_function (db, "GeometryType", 1, SQLITE_ANY, 0,
			     fnct_GeometryType, 0, 0);
    sqlite3_create_function (db, "ST_GeometryType", 1, SQLITE_ANY, 0,
			     fnct_GeometryType, 0, 0);
    sqlite3_create_function (db, "GeometryAliasType", 1, SQLITE_ANY, 0,
			     fnct_GeometryAliasType, 0, 0);
    sqlite3_create_function (db, "SridFromAuthCRS", 2, SQLITE_ANY, 0,
			     fnct_SridFromAuthCRS, 0, 0);
    sqlite3_create_function (db, "SRID", 1, SQLITE_ANY, 0, fnct_SRID, 0, 0);
    sqlite3_create_function (db, "ST_SRID", 1, SQLITE_ANY, 0, fnct_SRID, 0, 0);
    sqlite3_create_function (db, "SetSRID", 2, SQLITE_ANY, 0, fnct_SetSRID, 0,
			     0);
    sqlite3_create_function (db, "IsEmpty", 1, SQLITE_ANY, 0, fnct_IsEmpty, 0,
			     0);
    sqlite3_create_function (db, "ST_IsEmpty", 1, SQLITE_ANY, 0, fnct_IsEmpty,
			     0, 0);
    sqlite3_create_function (db, "Envelope", 1, SQLITE_ANY, 0, fnct_Envelope,
			     0, 0);
    sqlite3_create_function (db, "ST_Envelope", 1, SQLITE_ANY, 0,
			     fnct_Envelope, 0, 0);
    sqlite3_create_function (db, "X", 1, SQLITE_ANY, 0, fnct_X, 0, 0);
    sqlite3_create_function (db, "Y", 1, SQLITE_ANY, 0, fnct_Y, 0, 0);
    sqlite3_create_function (db, "Z", 1, SQLITE_ANY, 0, fnct_Z, 0, 0);
    sqlite3_create_function (db, "M", 1, SQLITE_ANY, 0, fnct_M, 0, 0);
    sqlite3_create_function (db, "ST_X", 1, SQLITE_ANY, 0, fnct_X, 0, 0);
    sqlite3_create_function (db, "ST_Y", 1, SQLITE_ANY, 0, fnct_Y, 0, 0);
    sqlite3_create_function (db, "ST_Z", 1, SQLITE_ANY, 0, fnct_Z, 0, 0);
    sqlite3_create_function (db, "ST_M", 1, SQLITE_ANY, 0, fnct_M, 0, 0);
    sqlite3_create_function (db, "NumPoints", 1, SQLITE_ANY, 0,
			     fnct_NumPoints, 0, 0);
    sqlite3_create_function (db, "ST_NumPoints", 1, SQLITE_ANY, 0,
			     fnct_NumPoints, 0, 0);
    sqlite3_create_function (db, "StartPoint", 1, SQLITE_ANY, 0,
			     fnct_StartPoint, 0, 0);
    sqlite3_create_function (db, "EndPoint", 1, SQLITE_ANY, 0, fnct_EndPoint,
			     0, 0);
    sqlite3_create_function (db, "ST_StartPoint", 1, SQLITE_ANY, 0,
			     fnct_StartPoint, 0, 0);
    sqlite3_create_function (db, "ST_EndPoint", 1, SQLITE_ANY, 0,
			     fnct_EndPoint, 0, 0);
    sqlite3_create_function (db, "PointN", 2, SQLITE_ANY, 0, fnct_PointN, 0, 0);
    sqlite3_create_function (db, "ST_PointN", 2, SQLITE_ANY, 0, fnct_PointN,
			     0, 0);
    sqlite3_create_function (db, "ExteriorRing", 1, SQLITE_ANY, 0,
			     fnct_ExteriorRing, 0, 0);
    sqlite3_create_function (db, "ST_ExteriorRing", 1, SQLITE_ANY, 0,
			     fnct_ExteriorRing, 0, 0);
    sqlite3_create_function (db, "NumInteriorRing", 1, SQLITE_ANY, 0,
			     fnct_NumInteriorRings, 0, 0);
    sqlite3_create_function (db, "NumInteriorRings", 1, SQLITE_ANY, 0,
			     fnct_NumInteriorRings, 0, 0);
    sqlite3_create_function (db, "ST_NumInteriorRing", 1, SQLITE_ANY, 0,
			     fnct_NumInteriorRings, 0, 0);
    sqlite3_create_function (db, "InteriorRingN", 2, SQLITE_ANY, 0,
			     fnct_InteriorRingN, 0, 0);
    sqlite3_create_function (db, "ST_InteriorRingN", 2, SQLITE_ANY, 0,
			     fnct_InteriorRingN, 0, 0);
    sqlite3_create_function (db, "NumGeometries", 1, SQLITE_ANY, 0,
			     fnct_NumGeometries, 0, 0);
    sqlite3_create_function (db, "ST_NumGeometries", 1, SQLITE_ANY, 0,
			     fnct_NumGeometries, 0, 0);
    sqlite3_create_function (db, "GeometryN", 2, SQLITE_ANY, 0,
			     fnct_GeometryN, 0, 0);
    sqlite3_create_function (db, "ST_GeometryN", 2, SQLITE_ANY, 0,
			     fnct_GeometryN, 0, 0);
    sqlite3_create_function (db, "MBRContains", 2, SQLITE_ANY, 0,
			     fnct_MbrContains, 0, 0);
    sqlite3_create_function (db, "MbrDisjoint", 2, SQLITE_ANY, 0,
			     fnct_MbrDisjoint, 0, 0);
    sqlite3_create_function (db, "MBREqual", 2, SQLITE_ANY, 0, fnct_MbrEqual,
			     0, 0);
    sqlite3_create_function (db, "MbrIntersects", 2, SQLITE_ANY, 0,
			     fnct_MbrIntersects, 0, 0);
    sqlite3_create_function (db, "MBROverlaps", 2, SQLITE_ANY, 0,
			     fnct_MbrOverlaps, 0, 0);
    sqlite3_create_function (db, "MbrTouches", 2, SQLITE_ANY, 0,
			     fnct_MbrTouches, 0, 0);
    sqlite3_create_function (db, "MbrWithin", 2, SQLITE_ANY, 0,
			     fnct_MbrWithin, 0, 0);
    sqlite3_create_function (db, "ShiftCoords", 3, SQLITE_ANY, 0,
			     fnct_ShiftCoords, 0, 0);
    sqlite3_create_function (db, "ShiftCoordinates", 3, SQLITE_ANY, 0,
			     fnct_ShiftCoords, 0, 0);
    sqlite3_create_function (db, "ScaleCoords", 2, SQLITE_ANY, 0,
			     fnct_ScaleCoords, 0, 0);
    sqlite3_create_function (db, "ScaleCoordinates", 2, SQLITE_ANY, 0,
			     fnct_ScaleCoords, 0, 0);
    sqlite3_create_function (db, "ScaleCoords", 3, SQLITE_ANY, 0,
			     fnct_ScaleCoords, 0, 0);
    sqlite3_create_function (db, "ScaleCoordinates", 3, SQLITE_ANY, 0,
			     fnct_ScaleCoords, 0, 0);
    sqlite3_create_function (db, "RotateCoords", 2, SQLITE_ANY, 0,
			     fnct_RotateCoords, 0, 0);
    sqlite3_create_function (db, "RotateCoordinates", 2, SQLITE_ANY, 0,
			     fnct_RotateCoords, 0, 0);
    sqlite3_create_function (db, "ReflectCoords", 3, SQLITE_ANY, 0,
			     fnct_ReflectCoords, 0, 0);
    sqlite3_create_function (db, "ReflectCoordinates", 3, SQLITE_ANY, 0,
			     fnct_ReflectCoords, 0, 0);
    sqlite3_create_function (db, "SwapCoords", 1, SQLITE_ANY, 0,
			     fnct_SwapCoords, 0, 0);
    sqlite3_create_function (db, "SwapCoordinates", 1, SQLITE_ANY, 0,
			     fnct_SwapCoords, 0, 0);
    sqlite3_create_function (db, "BuildMbr", 4, SQLITE_ANY, 0, fnct_BuildMbr1,
			     0, 0);
    sqlite3_create_function (db, "BuildMbr", 5, SQLITE_ANY, 0, fnct_BuildMbr2,
			     0, 0);
    sqlite3_create_function (db, "BuildCircleMbr", 3, SQLITE_ANY, 0,
			     fnct_BuildCircleMbr1, 0, 0);
    sqlite3_create_function (db, "BuildCircleMbr", 4, SQLITE_ANY, 0,
			     fnct_BuildCircleMbr2, 0, 0);
    sqlite3_create_function (db, "Extent", 1, SQLITE_ANY, 0, 0,
			     fnct_Extent_step, fnct_Extent_final);
    sqlite3_create_function (db, "MbrMinX", 1, SQLITE_ANY, 0, fnct_MbrMinX, 0,
			     0);
    sqlite3_create_function (db, "MbrMaxX", 1, SQLITE_ANY, 0, fnct_MbrMaxX, 0,
			     0);
    sqlite3_create_function (db, "MbrMinY", 1, SQLITE_ANY, 0, fnct_MbrMinY, 0,
			     0);
    sqlite3_create_function (db, "MbrMaxY", 1, SQLITE_ANY, 0, fnct_MbrMaxY, 0,
			     0);
    sqlite3_create_function (db, "MakePoint", 2, SQLITE_ANY, 0, fnct_MakePoint1,
			     0, 0);
    sqlite3_create_function (db, "MakePoint", 3, SQLITE_ANY, 0, fnct_MakePoint2,
			     0, 0);
    sqlite3_create_function (db, "MakeLine", 1, SQLITE_ANY, 0, 0,
			     fnct_MakeLine_step, fnct_MakeLine_final);
    sqlite3_create_function (db, "MakeLine", 2, SQLITE_ANY, 0, fnct_MakeLine, 0,
			     0);
    sqlite3_create_function (db, "Collect", 1, SQLITE_ANY, 0, 0,
			     fnct_Collect_step, fnct_Collect_final);
    sqlite3_create_function (db, "Collect", 2, SQLITE_ANY, 0, fnct_Collect, 0,
			     0);
    sqlite3_create_function (db, "ST_Collect", 1, SQLITE_ANY, 0, 0,
			     fnct_Collect_step, fnct_Collect_final);
    sqlite3_create_function (db, "ST_Collect", 2, SQLITE_ANY, 0, fnct_Collect,
			     0, 0);
    sqlite3_create_function (db, "BuildMbrFilter", 4, SQLITE_ANY, 0,
			     fnct_BuildMbrFilter, 0, 0);
    sqlite3_create_function (db, "FilterMbrWithin", 4, SQLITE_ANY, 0,
			     fnct_FilterMbrWithin, 0, 0);
    sqlite3_create_function (db, "FilterMbrContains", 4, SQLITE_ANY, 0,
			     fnct_FilterMbrContains, 0, 0);
    sqlite3_create_function (db, "FilterMbrIntersects", 4, SQLITE_ANY, 0,
			     fnct_FilterMbrIntersects, 0, 0);
    sqlite3_create_function (db, "LinesFromRings", 1, SQLITE_ANY, 0,
			     fnct_LinesFromRings, 0, 0);
    sqlite3_create_function (db, "LinesFromRings", 2, SQLITE_ANY, 0,
			     fnct_LinesFromRings, 0, 0);

#ifndef OMIT_GEOS		/* including GEOS */
    sqlite3_create_function (db, "BuildArea", 1, SQLITE_ANY, 0, fnct_BuildArea,
			     0, 0);
    sqlite3_create_function (db, "ST_BuildArea", 1, SQLITE_ANY, 0,
			     fnct_BuildArea, 0, 0);
    sqlite3_create_function (db, "Polygonize", 1, SQLITE_ANY, 0, 0,
			     fnct_Polygonize_step, fnct_Polygonize_final);
    sqlite3_create_function (db, "ST_Polygonize", 1, SQLITE_ANY, 0, 0,
			     fnct_Polygonize_step, fnct_Polygonize_final);
#endif /* end including GEOS */

    sqlite3_create_function (db, "DissolveSegments", 1, SQLITE_ANY, 0,
			     fnct_DissolveSegments, 0, 0);
    sqlite3_create_function (db, "ST_DissolveSegments", 1, SQLITE_ANY, 0,
			     fnct_DissolveSegments, 0, 0);
    sqlite3_create_function (db, "DissolvePoints", 1, SQLITE_ANY, 0,
			     fnct_DissolvePoints, 0, 0);
    sqlite3_create_function (db, "ST_DissolvePoints", 1, SQLITE_ANY, 0,
			     fnct_DissolvePoints, 0, 0);
    sqlite3_create_function (db, "CollectionExtract", 2, SQLITE_ANY, 0,
			     fnct_CollectionExtract, 0, 0);
    sqlite3_create_function (db, "ST_CollectionExtract", 2, SQLITE_ANY, 0,
			     fnct_CollectionExtract, 0, 0);
#ifndef OMIT_GEOCALLBACKS	/* supporting RTree geometry callbacks */
    sqlite3_rtree_geometry_callback (db, "RTreeWithin", fnct_RTreeIntersects,
				     0);
    sqlite3_rtree_geometry_callback (db, "RTreeContains", fnct_RTreeIntersects,
				     0);
    sqlite3_rtree_geometry_callback (db, "RTreeIntersects",
				     fnct_RTreeIntersects, 0);
    sqlite3_rtree_geometry_callback (db, "RTreeDistWithin",
				     fnct_RTreeDistWithin, 0);
#endif /* end RTree geometry callbacks */

/* some BLOB/JPEG/EXIF functions */
    sqlite3_create_function (db, "IsGeometryBlob", 1, SQLITE_ANY, 0,
			     fnct_IsGeometryBlob, 0, 0);
    sqlite3_create_function (db, "IsZipBlob", 1, SQLITE_ANY, 0,
			     fnct_IsZipBlob, 0, 0);
    sqlite3_create_function (db, "IsPdfBlob", 1, SQLITE_ANY, 0,
			     fnct_IsPdfBlob, 0, 0);
    sqlite3_create_function (db, "IsTiffBlob", 1, SQLITE_ANY, 0,
			     fnct_IsTiffBlob, 0, 0);
    sqlite3_create_function (db, "IsGifBlob", 1, SQLITE_ANY, 0,
			     fnct_IsGifBlob, 0, 0);
    sqlite3_create_function (db, "IsPngBlob", 1, SQLITE_ANY, 0,
			     fnct_IsPngBlob, 0, 0);
    sqlite3_create_function (db, "IsJpegBlob", 1, SQLITE_ANY, 0,
			     fnct_IsJpegBlob, 0, 0);
    sqlite3_create_function (db, "IsExifBlob", 1, SQLITE_ANY, 0,
			     fnct_IsExifBlob, 0, 0);
    sqlite3_create_function (db, "IsExifGpsBlob", 1, SQLITE_ANY, 0,
			     fnct_IsExifGpsBlob, 0, 0);
    sqlite3_create_function (db, "GeomFromExifGpsBlob", 1, SQLITE_ANY, 0,
			     fnct_GeomFromExifGpsBlob, 0, 0);

/* some Geodesic functions */
    sqlite3_create_function (db, "GreatCircleLength", 1, SQLITE_ANY, 0,
			     fnct_GreatCircleLength, 0, 0);
    sqlite3_create_function (db, "GeodesicLength", 1, SQLITE_ANY, 0,
			     fnct_GeodesicLength, 0, 0);

/* some Length Unit conversion functions */
    sqlite3_create_function (db, "CvtToKm", 1, SQLITE_ANY, 0, fnct_cvtToKm, 0,
			     0);
    sqlite3_create_function (db, "CvtToDm", 1, SQLITE_ANY, 0, fnct_cvtToDm, 0,
			     0);
    sqlite3_create_function (db, "CvtToCm", 1, SQLITE_ANY, 0, fnct_cvtToCm, 0,
			     0);
    sqlite3_create_function (db, "CvtToMm", 1, SQLITE_ANY, 0, fnct_cvtToMm, 0,
			     0);
    sqlite3_create_function (db, "CvtToKmi", 1, SQLITE_ANY, 0, fnct_cvtToKmi,
			     0, 0);
    sqlite3_create_function (db, "CvtToIn", 1, SQLITE_ANY, 0, fnct_cvtToIn, 0,
			     0);
    sqlite3_create_function (db, "CvtToFt", 1, SQLITE_ANY, 0, fnct_cvtToFt, 0,
			     0);
    sqlite3_create_function (db, "CvtToYd", 1, SQLITE_ANY, 0, fnct_cvtToYd, 0,
			     0);
    sqlite3_create_function (db, "CvtToMi", 1, SQLITE_ANY, 0, fnct_cvtToMi, 0,
			     0);
    sqlite3_create_function (db, "CvtToFath", 1, SQLITE_ANY, 0,
			     fnct_cvtToFath, 0, 0);
    sqlite3_create_function (db, "CvtToCh", 1, SQLITE_ANY, 0, fnct_cvtToCh, 0,
			     0);
    sqlite3_create_function (db, "CvtToLink", 1, SQLITE_ANY, 0,
			     fnct_cvtToLink, 0, 0);
    sqlite3_create_function (db, "CvtToUsIn", 1, SQLITE_ANY, 0,
			     fnct_cvtToUsIn, 0, 0);
    sqlite3_create_function (db, "CvtToUsFt", 1, SQLITE_ANY, 0,
			     fnct_cvtToUsFt, 0, 0);
    sqlite3_create_function (db, "CvtToUsYd", 1, SQLITE_ANY, 0,
			     fnct_cvtToUsYd, 0, 0);
    sqlite3_create_function (db, "CvtToUsCh", 1, SQLITE_ANY, 0,
			     fnct_cvtToUsCh, 0, 0);
    sqlite3_create_function (db, "CvtToUsMi", 1, SQLITE_ANY, 0,
			     fnct_cvtToUsMi, 0, 0);
    sqlite3_create_function (db, "CvtToIndFt", 1, SQLITE_ANY, 0,
			     fnct_cvtToIndFt, 0, 0);
    sqlite3_create_function (db, "CvtToIndYd", 1, SQLITE_ANY, 0,
			     fnct_cvtToIndYd, 0, 0);
    sqlite3_create_function (db, "CvtToIndCh", 1, SQLITE_ANY, 0,
			     fnct_cvtToIndCh, 0, 0);
    sqlite3_create_function (db, "CvtFromKm", 1, SQLITE_ANY, 0,
			     fnct_cvtFromKm, 0, 0);
    sqlite3_create_function (db, "CvtFromDm", 1, SQLITE_ANY, 0,
			     fnct_cvtFromDm, 0, 0);
    sqlite3_create_function (db, "CvtFromCm", 1, SQLITE_ANY, 0,
			     fnct_cvtFromCm, 0, 0);
    sqlite3_create_function (db, "CvtFromMm", 1, SQLITE_ANY, 0,
			     fnct_cvtFromMm, 0, 0);
    sqlite3_create_function (db, "CvtFromKmi", 1, SQLITE_ANY, 0,
			     fnct_cvtFromKmi, 0, 0);
    sqlite3_create_function (db, "CvtFromIn", 1, SQLITE_ANY, 0,
			     fnct_cvtFromIn, 0, 0);
    sqlite3_create_function (db, "CvtFromFt", 1, SQLITE_ANY, 0,
			     fnct_cvtFromFt, 0, 0);
    sqlite3_create_function (db, "CvtFromYd", 1, SQLITE_ANY, 0,
			     fnct_cvtFromYd, 0, 0);
    sqlite3_create_function (db, "CvtFromMi", 1, SQLITE_ANY, 0,
			     fnct_cvtFromMi, 0, 0);
    sqlite3_create_function (db, "CvtFromFath", 1, SQLITE_ANY, 0,
			     fnct_cvtFromFath, 0, 0);
    sqlite3_create_function (db, "CvtFromCh", 1, SQLITE_ANY, 0,
			     fnct_cvtFromCh, 0, 0);
    sqlite3_create_function (db, "CvtFromLink", 1, SQLITE_ANY, 0,
			     fnct_cvtFromLink, 0, 0);
    sqlite3_create_function (db, "CvtFromUsIn", 1, SQLITE_ANY, 0,
			     fnct_cvtFromUsIn, 0, 0);
    sqlite3_create_function (db, "CvtFromUsFt", 1, SQLITE_ANY, 0,
			     fnct_cvtFromUsFt, 0, 0);
    sqlite3_create_function (db, "CvtFromUsYd", 1, SQLITE_ANY, 0,
			     fnct_cvtFromUsYd, 0, 0);
    sqlite3_create_function (db, "CvtFromUsCh", 1, SQLITE_ANY, 0,
			     fnct_cvtFromUsCh, 0, 0);
    sqlite3_create_function (db, "CvtFromUsMi", 1, SQLITE_ANY, 0,
			     fnct_cvtFromUsMi, 0, 0);
    sqlite3_create_function (db, "CvtFromIndFt", 1, SQLITE_ANY, 0,
			     fnct_cvtFromIndFt, 0, 0);
    sqlite3_create_function (db, "CvtFromIndYd", 1, SQLITE_ANY, 0,
			     fnct_cvtFromIndYd, 0, 0);
    sqlite3_create_function (db, "CvtFromIndCh", 1, SQLITE_ANY, 0,
			     fnct_cvtFromIndCh, 0, 0);

#ifndef OMIT_MATHSQL		/* supporting SQL math functions */

/* some extra math functions */
    sqlite3_create_function (db, "acos", 1, SQLITE_ANY, 0, fnct_math_acos, 0,
			     0);
    sqlite3_create_function (db, "asin", 1, SQLITE_ANY, 0, fnct_math_asin, 0,
			     0);
    sqlite3_create_function (db, "atan", 1, SQLITE_ANY, 0, fnct_math_atan, 0,
			     0);
    sqlite3_create_function (db, "ceil", 1, SQLITE_ANY, 0, fnct_math_ceil, 0,
			     0);
    sqlite3_create_function (db, "ceiling", 1, SQLITE_ANY, 0, fnct_math_ceil,
			     0, 0);
    sqlite3_create_function (db, "cos", 1, SQLITE_ANY, 0, fnct_math_cos, 0, 0);
    sqlite3_create_function (db, "cot", 1, SQLITE_ANY, 0, fnct_math_cot, 0, 0);
    sqlite3_create_function (db, "degrees", 1, SQLITE_ANY, 0,
			     fnct_math_degrees, 0, 0);
    sqlite3_create_function (db, "exp", 1, SQLITE_ANY, 0, fnct_math_exp, 0, 0);
    sqlite3_create_function (db, "floor", 1, SQLITE_ANY, 0, fnct_math_floor,
			     0, 0);
    sqlite3_create_function (db, "ln", 1, SQLITE_ANY, 0, fnct_math_logn, 0, 0);
    sqlite3_create_function (db, "log", 1, SQLITE_ANY, 0, fnct_math_logn, 0, 0);
    sqlite3_create_function (db, "log", 2, SQLITE_ANY, 0, fnct_math_logn2, 0,
			     0);
    sqlite3_create_function (db, "log2", 1, SQLITE_ANY, 0, fnct_math_log_2, 0,
			     0);
    sqlite3_create_function (db, "log10", 1, SQLITE_ANY, 0, fnct_math_log_10,
			     0, 0);
    sqlite3_create_function (db, "pi", 0, SQLITE_ANY, 0, fnct_math_pi, 0, 0);
    sqlite3_create_function (db, "pow", 2, SQLITE_ANY, 0, fnct_math_pow, 0, 0);
    sqlite3_create_function (db, "power", 2, SQLITE_ANY, 0, fnct_math_pow, 0,
			     0);
    sqlite3_create_function (db, "radians", 1, SQLITE_ANY, 0,
			     fnct_math_radians, 0, 0);
    sqlite3_create_function (db, "round", 1, SQLITE_ANY, 0, fnct_math_round,
			     0, 0);
    sqlite3_create_function (db, "sign", 1, SQLITE_ANY, 0, fnct_math_sign, 0,
			     0);
    sqlite3_create_function (db, "sin", 1, SQLITE_ANY, 0, fnct_math_sin, 0, 0);
    sqlite3_create_function (db, "stddev_pop", 1, SQLITE_ANY, 0, 0,
			     fnct_math_stddev_step, fnct_math_stddev_pop_final);
    sqlite3_create_function (db, "stddev_samp", 1, SQLITE_ANY, 0, 0,
			     fnct_math_stddev_step,
			     fnct_math_stddev_samp_final);
    sqlite3_create_function (db, "sqrt", 1, SQLITE_ANY, 0, fnct_math_sqrt, 0,
			     0);
    sqlite3_create_function (db, "tan", 1, SQLITE_ANY, 0, fnct_math_tan, 0, 0);
    sqlite3_create_function (db, "var_pop", 1, SQLITE_ANY, 0, 0,
			     fnct_math_stddev_step, fnct_math_var_pop_final);
    sqlite3_create_function (db, "var_samp", 1, SQLITE_ANY, 0, 0,
			     fnct_math_stddev_step, fnct_math_var_samp_final);

#endif /* end supporting SQL math functions */

#ifndef OMIT_PROJ		/* including PROJ.4 */

    sqlite3_create_function (db, "Transform", 2, SQLITE_ANY, 0,
			     fnct_Transform, 0, 0);
    sqlite3_create_function (db, "ST_Transform", 2, SQLITE_ANY, 0,
			     fnct_Transform, 0, 0);

#endif /* end including PROJ.4 */

#if !(defined(OMIT_GEOS) && defined(OMIT_BOOSTGEOMETRY))		/* including GEOS */

    initGEOS (geos_warning, geos_error);
    sqlite3_create_function (db, "Boundary", 1, SQLITE_ANY, 0, fnct_Boundary,
			     0, 0);
    sqlite3_create_function (db, "ST_Boundary", 1, SQLITE_ANY, 0,
			     fnct_Boundary, 0, 0);
    sqlite3_create_function (db, "IsClosed", 1, SQLITE_ANY, 0, fnct_IsClosed,
			     0, 0);
    sqlite3_create_function (db, "ST_IsClosed", 1, SQLITE_ANY, 0,
			     fnct_IsClosed, 0, 0);
    sqlite3_create_function (db, "IsSimple", 1, SQLITE_ANY, 0, fnct_IsSimple,
			     0, 0);
    sqlite3_create_function (db, "ST_IsSimple", 1, SQLITE_ANY, 0,
			     fnct_IsSimple, 0, 0);
    sqlite3_create_function (db, "IsRing", 1, SQLITE_ANY, 0, fnct_IsRing, 0, 0);
    sqlite3_create_function (db, "ST_IsRing", 1, SQLITE_ANY, 0, fnct_IsRing,
			     0, 0);
    sqlite3_create_function (db, "IsValid", 1, SQLITE_ANY, 0, fnct_IsValid, 0,
			     0);
    sqlite3_create_function (db, "ST_IsValid", 1, SQLITE_ANY, 0, fnct_IsValid,
			     0, 0);
    sqlite3_create_function (db, "GLength", 1, SQLITE_ANY, 0, fnct_Length, 0,
			     0);
    sqlite3_create_function (db, "ST_Length", 1, SQLITE_ANY, 0, fnct_Length,
			     0, 0);
    sqlite3_create_function (db, "Area", 1, SQLITE_ANY, 0, fnct_Area, 0, 0);
    sqlite3_create_function (db, "ST_Area", 1, SQLITE_ANY, 0, fnct_Area, 0, 0);
    sqlite3_create_function (db, "Centroid", 1, SQLITE_ANY, 0, fnct_Centroid,
			     0, 0);
    sqlite3_create_function (db, "ST_Centroid", 1, SQLITE_ANY, 0,
			     fnct_Centroid, 0, 0);
    sqlite3_create_function (db, "PointOnSurface", 1, SQLITE_ANY, 0,
			     fnct_PointOnSurface, 0, 0);
    sqlite3_create_function (db, "ST_PointOnSurface", 1, SQLITE_ANY, 0,
			     fnct_PointOnSurface, 0, 0);
    sqlite3_create_function (db, "Simplify", 2, SQLITE_ANY, 0, fnct_Simplify,
			     0, 0);
    sqlite3_create_function (db, "ST_Generalize", 2, SQLITE_ANY, 0,
			     fnct_Simplify, 0, 0);
    sqlite3_create_function (db, "SimplifyPreserveTopology", 2, SQLITE_ANY, 0,
			     fnct_SimplifyPreserveTopology, 0, 0);
    sqlite3_create_function (db, "ConvexHull", 1, SQLITE_ANY, 0,
			     fnct_ConvexHull, 0, 0);
    sqlite3_create_function (db, "ST_ConvexHull", 1, SQLITE_ANY, 0,
			     fnct_ConvexHull, 0, 0);
    sqlite3_create_function (db, "Buffer", 2, SQLITE_ANY, 0, fnct_Buffer, 0, 0);
    sqlite3_create_function (db, "ST_Buffer", 2, SQLITE_ANY, 0, fnct_Buffer,
			     0, 0);
    sqlite3_create_function (db, "Intersection", 2, SQLITE_ANY, 0,
			     fnct_Intersection, 0, 0);
    sqlite3_create_function (db, "ST_Intersection", 2, SQLITE_ANY, 0,
			     fnct_Intersection, 0, 0);
    sqlite3_create_function (db, "GUnion", 1, SQLITE_ANY, 0, 0,
			     fnct_Union_step, fnct_Union_final);
    sqlite3_create_function (db, "GUnion", 2, SQLITE_ANY, 0, fnct_Union, 0, 0);
    sqlite3_create_function (db, "ST_Union", 1, SQLITE_ANY, 0, 0,
			     fnct_Union_step, fnct_Union_final);
    sqlite3_create_function (db, "ST_Union", 2, SQLITE_ANY, 0, fnct_Union, 0,
			     0);
    sqlite3_create_function (db, "Difference", 2, SQLITE_ANY, 0,
			     fnct_Difference, 0, 0);
    sqlite3_create_function (db, "ST_Difference", 2, SQLITE_ANY, 0,
			     fnct_Difference, 0, 0);
    sqlite3_create_function (db, "SymDifference", 2, SQLITE_ANY, 0,
			     fnct_SymDifference, 0, 0);
    sqlite3_create_function (db, "ST_SymDifference", 2, SQLITE_ANY, 0,
			     fnct_SymDifference, 0, 0);
    sqlite3_create_function (db, "Equals", 2, SQLITE_ANY, 0, fnct_Equals, 0, 0);
    sqlite3_create_function (db, "ST_Equals", 2, SQLITE_ANY, 0, fnct_Equals,
			     0, 0);
    sqlite3_create_function (db, "Intersects", 2, SQLITE_ANY, 0,
			     fnct_Intersects, 0, 0);
    sqlite3_create_function (db, "ST_Intersects", 2, SQLITE_ANY, 0,
			     fnct_Intersects, 0, 0);
    sqlite3_create_function (db, "Disjoint", 2, SQLITE_ANY, 0, fnct_Disjoint,
			     0, 0);
    sqlite3_create_function (db, "ST_Disjoint", 2, SQLITE_ANY, 0,
			     fnct_Disjoint, 0, 0);
    sqlite3_create_function (db, "Overlaps", 2, SQLITE_ANY, 0, fnct_Overlaps,
			     0, 0);
    sqlite3_create_function (db, "ST_Overlaps", 2, SQLITE_ANY, 0,
			     fnct_Overlaps, 0, 0);
    sqlite3_create_function (db, "Crosses", 2, SQLITE_ANY, 0, fnct_Crosses, 0,
			     0);
    sqlite3_create_function (db, "ST_Crosses", 2, SQLITE_ANY, 0, fnct_Crosses,
			     0, 0);
    sqlite3_create_function (db, "Touches", 2, SQLITE_ANY, 0, fnct_Touches, 0,
			     0);
    sqlite3_create_function (db, "ST_Touches", 2, SQLITE_ANY, 0, fnct_Touches,
			     0, 0);
    sqlite3_create_function (db, "Within", 2, SQLITE_ANY, 0, fnct_Within, 0, 0);
    sqlite3_create_function (db, "ST_Within", 2, SQLITE_ANY, 0, fnct_Within,
			     0, 0);
    sqlite3_create_function (db, "Contains", 2, SQLITE_ANY, 0, fnct_Contains,
			     0, 0);
    sqlite3_create_function (db, "ST_Contains", 2, SQLITE_ANY, 0,
			     fnct_Contains, 0, 0);
    sqlite3_create_function (db, "Relate", 3, SQLITE_ANY, 0, fnct_Relate, 0, 0);
    sqlite3_create_function (db, "ST_Relate", 3, SQLITE_ANY, 0, fnct_Relate,
			     0, 0);
    sqlite3_create_function (db, "Distance", 2, SQLITE_ANY, 0, fnct_Distance,
			     0, 0);
    sqlite3_create_function (db, "ST_Distance", 2, SQLITE_ANY, 0,
			     fnct_Distance, 0, 0);
    sqlite3_create_function (db, "PtDistWithin", 3, SQLITE_ANY, 0,
			     fnct_PtDistWithin, 0, 0);
    sqlite3_create_function (db, "PtDistWithin", 4, SQLITE_ANY, 0,
			     fnct_PtDistWithin, 0, 0);
    sqlite3_create_function (db, "BdPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_BdPolyFromText1, 0, 0);
    sqlite3_create_function (db, "BdPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_BdPolyFromText2, 0, 0);
    sqlite3_create_function (db, "BdMPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_BdMPolyFromText1, 0, 0);
    sqlite3_create_function (db, "BdMPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_BdMPolyFromText2, 0, 0);
    sqlite3_create_function (db, "BdPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_BdPolyFromWKB1, 0, 0);
    sqlite3_create_function (db, "BdPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_BdPolyFromWKB2, 0, 0);
    sqlite3_create_function (db, "BdMPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_BdMPolyFromWKB1, 0, 0);
    sqlite3_create_function (db, "BdMPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_BdMPolyFromWKB2, 0, 0);
    sqlite3_create_function (db, "ST_BdPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_BdPolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_BdPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_BdPolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_BdMPolyFromText", 1, SQLITE_ANY, 0,
			     fnct_BdMPolyFromText1, 0, 0);
    sqlite3_create_function (db, "ST_BdMPolyFromText", 2, SQLITE_ANY, 0,
			     fnct_BdMPolyFromText2, 0, 0);
    sqlite3_create_function (db, "ST_BdPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_BdPolyFromWKB1, 0, 0);
    sqlite3_create_function (db, "ST_BdPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_BdPolyFromWKB2, 0, 0);
    sqlite3_create_function (db, "ST_BdMPolyFromWKB", 1, SQLITE_ANY, 0,
			     fnct_BdMPolyFromWKB1, 0, 0);
    sqlite3_create_function (db, "ST_BdMPolyFromWKB", 2, SQLITE_ANY, 0,
			     fnct_BdMPolyFromWKB2, 0, 0);

#ifdef GEOS_ADVANCED		/* GEOS advanced and experimental features */

    sqlite3_create_function (db, "OffsetCurve", 3, SQLITE_ANY, 0,
			     fnct_OffsetCurve, 0, 0);
    sqlite3_create_function (db, "ST_OffsetCurve", 3, SQLITE_ANY, 0,
			     fnct_OffsetCurve, 0, 0);
    sqlite3_create_function (db, "SingleSidedBuffer", 3, SQLITE_ANY, 0,
			     fnct_SingleSidedBuffer, 0, 0);
    sqlite3_create_function (db, "ST_SingleSidedBuffer", 3, SQLITE_ANY, 0,
			     fnct_SingleSidedBuffer, 0, 0);
    sqlite3_create_function (db, "HausdorffDistance", 2, SQLITE_ANY, 0,
			     fnct_HausdorffDistance, 0, 0);
    sqlite3_create_function (db, "ST_HausdorffDistance", 2, SQLITE_ANY, 0,
			     fnct_HausdorffDistance, 0, 0);
    sqlite3_create_function (db, "SharedPaths", 2, SQLITE_ANY, 0,
			     fnct_SharedPaths, 0, 0);
    sqlite3_create_function (db, "ST_SharedPaths", 2, SQLITE_ANY, 0,
			     fnct_SharedPaths, 0, 0);
    sqlite3_create_function (db, "Covers", 2, SQLITE_ANY, 0, fnct_Covers, 0, 0);
    sqlite3_create_function (db, "ST_Covers", 2, SQLITE_ANY, 0, fnct_Covers, 0,
			     0);
    sqlite3_create_function (db, "CoveredBy", 2, SQLITE_ANY, 0, fnct_CoveredBy,
			     0, 0);
    sqlite3_create_function (db, "ST_CoveredBy", 2, SQLITE_ANY, 0,
			     fnct_CoveredBy, 0, 0);
    sqlite3_create_function (db, "Line_Interpolate_Point", 2, SQLITE_ANY, 0,
			     fnct_LineInterpolatePoint, 0, 0);
    sqlite3_create_function (db, "ST_Line_Interpolate_Point", 2, SQLITE_ANY, 0,
			     fnct_LineInterpolatePoint, 0, 0);
    sqlite3_create_function (db, "Line_Locate_Point", 2, SQLITE_ANY, 0,
			     fnct_LineLocatePoint, 0, 0);
    sqlite3_create_function (db, "ST_Line_Locate_Point", 2, SQLITE_ANY, 0,
			     fnct_LineLocatePoint, 0, 0);
    sqlite3_create_function (db, "Line_Substring", 3, SQLITE_ANY, 0,
			     fnct_LineSubstring, 0, 0);
    sqlite3_create_function (db, "ST_Line_Substring", 3, SQLITE_ANY, 0,
			     fnct_LineSubstring, 0, 0);
    sqlite3_create_function (db, "ClosestPoint", 2, SQLITE_ANY, 0,
			     fnct_ClosestPoint, 0, 0);
    sqlite3_create_function (db, "ST_ClosestPoint", 2, SQLITE_ANY, 0,
			     fnct_ClosestPoint, 0, 0);
    sqlite3_create_function (db, "ShortestLine", 2, SQLITE_ANY, 0,
			     fnct_ShortestLine, 0, 0);
    sqlite3_create_function (db, "ST_ShortestLine", 2, SQLITE_ANY, 0,
			     fnct_ShortestLine, 0, 0);
    sqlite3_create_function (db, "Snap", 3, SQLITE_ANY, 0, fnct_Snap, 0, 0);
    sqlite3_create_function (db, "ST_Snap", 3, SQLITE_ANY, 0, fnct_Snap, 0, 0);
    sqlite3_create_function (db, "LineMerge", 1, SQLITE_ANY, 0,
			     fnct_LineMerge, 0, 0);
    sqlite3_create_function (db, "ST_LineMerge", 1, SQLITE_ANY, 0,
			     fnct_LineMerge, 0, 0);
    sqlite3_create_function (db, "UnaryUnion", 1, SQLITE_ANY, 0,
			     fnct_UnaryUnion, 0, 0);
    sqlite3_create_function (db, "ST_UnaryUnion", 1, SQLITE_ANY, 0,
			     fnct_UnaryUnion, 0, 0);

#endif /* end GEOS advanced and experimental features */

#endif /* end including GEOS */
}

static void
init_spatialite_virtualtables (sqlite3 * db)
{
#ifndef OMIT_ICONV		/* when ICONV is disabled SHP/DBF/TXT cannot be supported */
/* initializing the VirtualShape  extension */
    virtualshape_extension_init (db);
/* initializing the VirtualDbf  extension */
    virtualdbf_extension_init (db);
/* initializing the VirtualText extension */
    virtualtext_extension_init (db);
#ifndef OMIT_FREEXL
/* initializing the VirtualXL  extension */
    virtualXL_extension_init (db);
#endif /* FreeXL enabled/disable */
#endif /* ICONV enabled/disabled */

/* initializing the VirtualNetwork  extension */
    virtualnetwork_extension_init (db);
/* initializing the MbrCache  extension */
    mbrcache_extension_init (db);
/* initializing the VirtualFDO  extension */
    virtualfdo_extension_init (db);
/* initializing the VirtualSpatialIndex  extension */
    virtual_spatialindex_extension_init (db);
}

static void
init_static_spatialite (sqlite3 * db, char **pzErrMsg,
			const sqlite3_api_routines * pApi)
{
    SQLITE_EXTENSION_INIT2 (pApi);
/* setting the POSIX locale for numeric */
    setlocale (LC_NUMERIC, "POSIX");
    *pzErrMsg = NULL;

    register_spatialite_sql_functions (db);

    init_spatialite_virtualtables (db);

/* setting a timeout handler */
    sqlite3_busy_timeout (db, 5000);
}

void
spatialite_init (int verbose)
{
/* used when SQLite initializes SpatiaLite via statically linked lib */
    sqlite3_auto_extension ((void (*)(void)) init_static_spatialite);
    if (isatty (1))
      {
	  /* printing "hello" message only when stdout is on console */
	  if (verbose)
	    {
		printf ("SpatiaLite version ..: %s", spatialite_version ());
		printf ("\tSupported Extensions:\n");
#ifndef OMIT_ICONV		/* ICONV is required by SHP/DBF/TXT */
		printf ("\t- 'VirtualShape'\t[direct Shapefile access]\n");
		printf ("\t- 'VirtualDbf'\t\t[direct DBF access]\n");
#ifndef OMIT_FREEXL
		printf ("\t- 'VirtualXL'\t\t[direct XLS access]\n");
#endif /* end FreeXL conditional */
		printf ("\t- 'VirtualText'\t\t[direct CSV/TXT access]\n");
#endif /* end ICONV conditional */
		printf ("\t- 'VirtualNetwork'\t[Dijkstra shortest path]\n");
		printf ("\t- 'RTree'\t\t[Spatial Index - R*Tree]\n");
		printf ("\t- 'MbrCache'\t\t[Spatial Index - MBR cache]\n");
		printf ("\t- 'VirtualSpatialIndex'\t[R*Tree metahandler]\n");
		printf ("\t- 'VirtualFDO'\t\t[FDO-OGR interoperability]\n");
		printf ("\t- 'SpatiaLite'\t\t[Spatial SQL - OGC]\n");
	    }
#ifndef OMIT_PROJ		/* PROJ.4 version */
	  if (verbose)
	      printf ("PROJ.4 version ......: %s\n", pj_get_release ());
#endif /* end including PROJ.4 */
#ifndef OMIT_GEOS		/* GEOS version */
	  if (verbose)
	      printf ("GEOS version ........: %s\n", GEOSversion ());
#endif /* end GEOS version */
      }
}

void
spatialite_cleanup ()
{
#ifndef OMIT_GEOS
    finishGEOS ();
#endif
}

SPATIALITE_DECLARE int
sqlite3_extension_init (sqlite3 * db, char **pzErrMsg,
			const sqlite3_api_routines * pApi)
{
/* SQLite invokes this routine once when it dynamically loads the extension. */
    SQLITE_EXTENSION_INIT2 (pApi);
    setlocale (LC_NUMERIC, "POSIX");
    *pzErrMsg = NULL;

    register_spatialite_sql_functions (db);

    init_spatialite_virtualtables (db);

    /* setting a timeout handler */
    sqlite3_busy_timeout (db, 5000);
    return 0;
}

SPATIALITE_DECLARE sqlite3_int64
math_llabs (sqlite3_int64 value)
{
/* replacing the C99 llabs() function */
    return value < 0 ? -value : value;
}

SPATIALITE_DECLARE double
math_round (double value)
{
/* replacing the C99 round() function */
    double min = floor (value);
    if (fabs (value - min) < 0.5)
	return min;
    return min + 1.0;
}

static int
check_layer_statistics (sqlite3 * sqlite)
{
/*
/ checks the LAYER_STATISTICS table for validity;
/ if the table doesn't exist, attempts to create
*/
    char sql[8192];
    char **results;
    int rows;
    int columns;
    int ret;
    int raster_layer = 0;
    int table_name = 0;
    int geometry_column = 0;
    int row_count = 0;
    int extent_min_x = 0;
    int extent_min_y = 0;
    int extent_max_x = 0;
    int extent_max_y = 0;
    int i;
    const char *name;

/* checking the LAYER_STATISTICS table */
    strcpy (sql, "PRAGMA table_info(layer_statistics)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return 0;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "raster_layer") == 0)
		    raster_layer = 1;
		if (strcasecmp (name, "table_name") == 0)
		    table_name = 1;
		if (strcasecmp (name, "geometry_column") == 0)
		    geometry_column = 1;
		if (strcasecmp (name, "row_count") == 0)
		    row_count = 1;
		if (strcasecmp (name, "extent_min_x") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_min_y") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_max_x") == 0)
		    extent_max_x = 1;
		if (strcasecmp (name, "extent_max_y") == 0)
		    extent_max_y = 1;
	    }
      }
    sqlite3_free_table (results);

/* LAYER_STATISTICS already exists and has a valid layout */
    if (raster_layer && table_name && geometry_column && row_count
	&& extent_min_x && extent_max_x && extent_max_y)
	return 1;
/* LAYER_STATISTICS already exists, but has an invalid layout */
    if (raster_layer || table_name || geometry_column || row_count
	|| extent_min_x || extent_max_x || extent_max_y)
	return 0;

/* attempting to create LAYER_STATISTICS */
    strcpy (sql, "CREATE TABLE layer_statistics (\n");
    strcat (sql, "raster_layer INTEGER NOT NULL,\n");
    strcat (sql, "table_name TEXT NOT NULL,\n");
    strcat (sql, "geometry_column TEXT NOT NULL,\n");
    strcat (sql, "row_count INTEGER,\n");
    strcat (sql, "extent_min_x DOUBLE,\n");
    strcat (sql, "extent_min_y DOUBLE,\n");
    strcat (sql, "extent_max_x DOUBLE,\n");
    strcat (sql, "extent_max_y DOUBLE,\n");
    strcat (sql, "CONSTRAINT pk_layer_statistics PRIMARY KEY ");
    strcat (sql, "(raster_layer, table_name, geometry_column),\n");
    strcat (sql, "CONSTRAINT fk_layer_statistics FOREIGN KEY ");
    strcat (sql, "(table_name, geometry_column) REFERENCES ");
    strcat (sql, "geometry_columns (f_table_name, f_geometry_column) ");
    strcat (sql, "ON DELETE CASCADE)");
    ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
    if (ret != SQLITE_OK)
	return 0;
    return 1;
}

static int
do_update_layer_statistics (sqlite3 * sqlite, const char *table,
			    const char *column, int count, int has_coords,
			    double min_x, double min_y, double max_x,
			    double max_y)
{
/* update LAYER_STATISTICS [single table/geometry] */
    char sql[8192];
    int ret;
    int error = 0;
    sqlite3_stmt *stmt;

    if (!check_layer_statistics (sqlite))
	return 0;
    strcpy (sql, "INSERT OR REPLACE INTO layer_statistics ");
    strcat (sql, "(raster_layer, table_name, geometry_column, ");
    strcat (sql, "row_count, extent_min_x, extent_min_y, ");
    strcat (sql, "extent_max_x, extent_max_y) ");
    strcat (sql, "VALUES (0, ?, ?, ?, ?, ?, ?, ?)");

/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
	return 0;

/* binding INSERT params */
    sqlite3_reset (stmt);
    sqlite3_clear_bindings (stmt);
    sqlite3_bind_text (stmt, 1, table, strlen (table), SQLITE_STATIC);
    sqlite3_bind_text (stmt, 2, column, strlen (column), SQLITE_STATIC);
    sqlite3_bind_int (stmt, 3, count);
    if (has_coords)
      {
	  sqlite3_bind_double (stmt, 4, min_x);
	  sqlite3_bind_double (stmt, 5, min_y);
	  sqlite3_bind_double (stmt, 6, max_x);
	  sqlite3_bind_double (stmt, 7, max_y);
      }
    else
      {
	  sqlite3_bind_null (stmt, 4);
	  sqlite3_bind_null (stmt, 5);
	  sqlite3_bind_null (stmt, 6);
	  sqlite3_bind_null (stmt, 7);
      }
    ret = sqlite3_step (stmt);
    if (ret == SQLITE_DONE || ret == SQLITE_ROW)
	;
    else
	error = 1;
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
	return 0;
    if (error)
	return 0;
    return 1;
}

static int
check_views_layer_statistics (sqlite3 * sqlite)
{
/*
/ checks the VIEWS_LAYER_STATISTICS table for validity;
/ if the table doesn't exist, attempts to create
*/
    char sql[8192];
    char **results;
    int rows;
    int columns;
    int ret;
    int view_name = 0;
    int view_geometry = 0;
    int row_count = 0;
    int extent_min_x = 0;
    int extent_min_y = 0;
    int extent_max_x = 0;
    int extent_max_y = 0;
    int i;
    const char *name;

/* checking the VIEWS_LAYER_STATISTICS table */
    strcpy (sql, "PRAGMA table_info(views_layer_statistics)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return 0;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "view_name") == 0)
		    view_name = 1;
		if (strcasecmp (name, "view_geometry") == 0)
		    view_geometry = 1;
		if (strcasecmp (name, "row_count") == 0)
		    row_count = 1;
		if (strcasecmp (name, "extent_min_x") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_min_y") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_max_x") == 0)
		    extent_max_x = 1;
		if (strcasecmp (name, "extent_max_y") == 0)
		    extent_max_y = 1;
	    }
      }
    sqlite3_free_table (results);

/* VIEWS_LAYER_STATISTICS already exists and has a valid layout */
    if (view_name && view_geometry && row_count && extent_min_x && extent_max_x
	&& extent_max_y)
	return 1;
/* VIEWS_LAYER_STATISTICS already exists, but has an invalid layout */
    if (view_name || view_geometry || row_count || extent_min_x || extent_max_x
	|| extent_max_y)
	return 0;

/* attempting to create VIEWS_LAYER_STATISTICS */
    strcpy (sql, "CREATE TABLE views_layer_statistics (\n");
    strcat (sql, "view_name TEXT NOT NULL,\n");
    strcat (sql, "view_geometry TEXT NOT NULL,\n");
    strcat (sql, "row_count INTEGER,\n");
    strcat (sql, "extent_min_x DOUBLE,\n");
    strcat (sql, "extent_min_y DOUBLE,\n");
    strcat (sql, "extent_max_x DOUBLE,\n");
    strcat (sql, "extent_max_y DOUBLE,\n");
    strcat (sql, "CONSTRAINT pk_views_layer_statistics PRIMARY KEY ");
    strcat (sql, "(view_name, view_geometry),\n");
    strcat (sql, "CONSTRAINT fk_views_layer_statistics FOREIGN KEY ");
    strcat (sql, "(view_name, view_geometry) REFERENCES ");
    strcat (sql, "views_geometry_columns (view_name, view_geometry) ");
    strcat (sql, "ON DELETE CASCADE)");
    ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
    if (ret != SQLITE_OK)
	return 0;
    return 1;
}

static int
do_update_views_layer_statistics (sqlite3 * sqlite, const char *table,
				  const char *column, int count, int has_coords,
				  double min_x, double min_y, double max_x,
				  double max_y)
{
/* update VIEWS_LAYER_STATISTICS [single table/geometry] */
    char sql[8192];
    int ret;
    int error = 0;
    sqlite3_stmt *stmt;

    if (!check_views_layer_statistics (sqlite))
	return 0;
    strcpy (sql, "INSERT OR REPLACE INTO views_layer_statistics ");
    strcat (sql, "(view_name, view_geometry, ");
    strcat (sql, "row_count, extent_min_x, extent_min_y, ");
    strcat (sql, "extent_max_x, extent_max_y) ");
    strcat (sql, "VALUES (?, ?, ?, ?, ?, ?, ?)");

/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
	return 0;

/* binding INSERT params */
    sqlite3_reset (stmt);
    sqlite3_clear_bindings (stmt);
    sqlite3_bind_text (stmt, 1, table, strlen (table), SQLITE_STATIC);
    sqlite3_bind_text (stmt, 2, column, strlen (column), SQLITE_STATIC);
    sqlite3_bind_int (stmt, 3, count);
    if (has_coords)
      {
	  sqlite3_bind_double (stmt, 4, min_x);
	  sqlite3_bind_double (stmt, 5, min_y);
	  sqlite3_bind_double (stmt, 6, max_x);
	  sqlite3_bind_double (stmt, 7, max_y);
      }
    else
      {
	  sqlite3_bind_null (stmt, 4);
	  sqlite3_bind_null (stmt, 5);
	  sqlite3_bind_null (stmt, 6);
	  sqlite3_bind_null (stmt, 7);
      }
    ret = sqlite3_step (stmt);
    if (ret == SQLITE_DONE || ret == SQLITE_ROW)
	;
    else
	error = 1;
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
	return 0;
    if (error)
	return 0;
    return 1;
}

static int
check_virts_layer_statistics (sqlite3 * sqlite)
{
/*
/ checks the VIRTS_LAYER_STATISTICS table for validity;
/ if the table doesn't exist, attempts to create
*/
    char sql[8192];
    char **results;
    int rows;
    int columns;
    int ret;
    int virt_name = 0;
    int virt_geometry = 0;
    int row_count = 0;
    int extent_min_x = 0;
    int extent_min_y = 0;
    int extent_max_x = 0;
    int extent_max_y = 0;
    int i;
    const char *name;

/* checking the VIRTS_LAYER_STATISTICS table */
    strcpy (sql, "PRAGMA table_info(virts_layer_statistics)");
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return 0;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		name = results[(i * columns) + 1];
		if (strcasecmp (name, "virt_name") == 0)
		    virt_name = 1;
		if (strcasecmp (name, "virt_geometry") == 0)
		    virt_geometry = 1;
		if (strcasecmp (name, "row_count") == 0)
		    row_count = 1;
		if (strcasecmp (name, "extent_min_x") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_min_y") == 0)
		    extent_min_x = 1;
		if (strcasecmp (name, "extent_max_x") == 0)
		    extent_max_x = 1;
		if (strcasecmp (name, "extent_max_y") == 0)
		    extent_max_y = 1;
	    }
      }
    sqlite3_free_table (results);

/* VIRTS_LAYER_STATISTICS already exists and has a valid layout */
    if (virt_name && virt_geometry && row_count && extent_min_x && extent_max_x
	&& extent_max_y)
	return 1;
/* VIRTS_LAYER_STATISTICS already exists, but has an invalid layout */
    if (virt_name || virt_geometry || row_count || extent_min_x || extent_max_x
	|| extent_max_y)
	return 0;

/* attempting to create VIRTS_LAYER_STATISTICS */
    strcpy (sql, "CREATE TABLE virts_layer_statistics (\n");
    strcat (sql, "virt_name TEXT NOT NULL,\n");
    strcat (sql, "virt_geometry TEXT NOT NULL,\n");
    strcat (sql, "row_count INTEGER,\n");
    strcat (sql, "extent_min_x DOUBLE,\n");
    strcat (sql, "extent_min_y DOUBLE,\n");
    strcat (sql, "extent_max_x DOUBLE,\n");
    strcat (sql, "extent_max_y DOUBLE,\n");
    strcat (sql, "CONSTRAINT pk_virts_layer_statistics PRIMARY KEY ");
    strcat (sql, "(virt_name, virt_geometry),\n");
    strcat (sql, "CONSTRAINT fk_virts_layer_statistics FOREIGN KEY ");
    strcat (sql, "(virt_name, virt_geometry) REFERENCES ");
    strcat (sql, "virts_geometry_columns (virt_name, virt_geometry) ");
    strcat (sql, "ON DELETE CASCADE)");
    ret = sqlite3_exec (sqlite, sql, NULL, 0, NULL);
    if (ret != SQLITE_OK)
	return 0;
    return 1;
}

static int
do_update_virts_layer_statistics (sqlite3 * sqlite, const char *table,
				  const char *column, int count, int has_coords,
				  double min_x, double min_y, double max_x,
				  double max_y)
{
/* update VIRTS_LAYER_STATISTICS [single table/geometry] */
    char sql[8192];
    int ret;
    int error = 0;
    sqlite3_stmt *stmt;

    if (!check_virts_layer_statistics (sqlite))
	return 0;
    strcpy (sql, "INSERT OR REPLACE INTO virts_layer_statistics ");
    strcat (sql, "(virt_name, virt_geometry, ");
    strcat (sql, "row_count, extent_min_x, extent_min_y, ");
    strcat (sql, "extent_max_x, extent_max_y) ");
    strcat (sql, "VALUES (?, ?, ?, ?, ?, ?, ?)");

/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
	return 0;

/* binding INSERT params */
    sqlite3_reset (stmt);
    sqlite3_clear_bindings (stmt);
    sqlite3_bind_text (stmt, 1, table, strlen (table), SQLITE_STATIC);
    sqlite3_bind_text (stmt, 2, column, strlen (column), SQLITE_STATIC);
    sqlite3_bind_int (stmt, 3, count);
    if (has_coords)
      {
	  sqlite3_bind_double (stmt, 4, min_x);
	  sqlite3_bind_double (stmt, 5, min_y);
	  sqlite3_bind_double (stmt, 6, max_x);
	  sqlite3_bind_double (stmt, 7, max_y);
      }
    else
      {
	  sqlite3_bind_null (stmt, 4);
	  sqlite3_bind_null (stmt, 5);
	  sqlite3_bind_null (stmt, 6);
	  sqlite3_bind_null (stmt, 7);
      }
    ret = sqlite3_step (stmt);
    if (ret == SQLITE_DONE || ret == SQLITE_ROW)
	;
    else
	error = 1;
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
	return 0;
    if (error)
	return 0;
    return 1;
}

#define SPATIALITE_STATISTICS_GENUINE	1
#define SPATIALITE_STATISTICS_VIEWS	2
#define SPATIALITE_STATISTICS_VIRTS	3

static int
do_compute_layer_statistics (sqlite3 * sqlite, const char *table,
			     const char *column, int stat_type)
{
/* computes LAYER_STATISTICS [single table/geometry] */
    char xtable[1024];
    char xgeom[1024];
    char sql[8192];
    char sql2[2048];
    int ret;
    int error = 0;
    int count;
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    int has_coords = 1;
    sqlite3_stmt *stmt;

    strcpy (xtable, table);
    double_quoted_sql (xtable);
    strcpy (xgeom, column);
    double_quoted_sql (xgeom);
    strcpy (sql, "SELECT Count(*), ");
    sprintf (sql2, "Min(MbrMinX(%s)), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "Min(MbrMinY(%s)), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "Max(MbrMaxX(%s)), ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "Max(MbrMaxY(%s)) ", xgeom);
    strcat (sql, sql2);
    sprintf (sql2, "FROM %s", xtable);
    strcat (sql, sql2);

/* compiling SQL prepared statement */
    ret = sqlite3_prepare_v2 (sqlite, sql, strlen (sql), &stmt, NULL);
    if (ret != SQLITE_OK)
	return 0;
    while (1)
      {
	  /* scrolling the result set rows */
	  ret = sqlite3_step (stmt);
	  if (ret == SQLITE_DONE)
	      break;		/* end of result set */
	  if (ret == SQLITE_ROW)
	    {
		count = sqlite3_column_int (stmt, 0);
		if (sqlite3_column_type (stmt, 1) == SQLITE_NULL)
		    has_coords = 0;
		else
		    min_x = sqlite3_column_double (stmt, 1);
		if (sqlite3_column_type (stmt, 2) == SQLITE_NULL)
		    has_coords = 0;
		else
		    min_y = sqlite3_column_double (stmt, 2);
		if (sqlite3_column_type (stmt, 3) == SQLITE_NULL)
		    has_coords = 0;
		else
		    max_x = sqlite3_column_double (stmt, 3);
		if (sqlite3_column_type (stmt, 4) == SQLITE_NULL)
		    has_coords = 0;
		else
		    max_y = sqlite3_column_double (stmt, 4);
		switch (stat_type)
		  {
		  case SPATIALITE_STATISTICS_GENUINE:
		      if (!do_update_layer_statistics
			  (sqlite, table, column, count, has_coords, min_x,
			   min_y, max_x, max_y))
			  error = 1;
		      break;
		  case SPATIALITE_STATISTICS_VIEWS:
		      if (!do_update_views_layer_statistics
			  (sqlite, table, column, count, has_coords, min_x,
			   min_y, max_x, max_y))
			  error = 1;
		      break;
		  case SPATIALITE_STATISTICS_VIRTS:
		      if (!do_update_virts_layer_statistics
			  (sqlite, table, column, count, has_coords, min_x,
			   min_y, max_x, max_y))
			  error = 1;
		      break;
		  };
	    }
	  else
	      error = 1;
      }
    ret = sqlite3_finalize (stmt);
    if (ret != SQLITE_OK)
	return 0;
    if (error)
	return 0;
    return 1;
}

static int
genuine_layer_statistics (sqlite3 * sqlite, const char *table,
			  const char *column)
{
/* updating genuine LAYER_STATISTICS metadata */
    int count = 0;
    char xtable[1024];
    char xgeom[1024];
    char sql[8192];
    char sql2[2048];
    int ret;
    const char *f_table_name;
    const char *f_geometry_column;
    int i;
    char **results;
    int rows;
    int columns;
    int error = 0;

    if (table == NULL && column == NULL)
      {
	  /* processing any table/geometry found in GEOMETRY_COLUMNS */
	  strcpy (sql, "SELECT f_table_name, f_geometry_column ");
	  strcat (sql, "FROM geometry_columns");
      }
    else if (column == NULL)
      {
	  /* processing any geometry belonging to this table */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (sql, "SELECT f_table_name, f_geometry_column ");
	  strcat (sql, "FROM geometry_columns ");
	  sprintf (sql2, "WHERE f_table_name LIKE '%s'", xtable);
	  strcat (sql, sql2);
      }
    else
      {
	  /* processing a single table/geometry entry */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (xgeom, column);
	  clean_sql_string (xgeom);
	  strcpy (sql, "SELECT f_table_name, f_geometry_column ");
	  strcat (sql, "FROM geometry_columns ");
	  sprintf (sql2, "WHERE f_table_name LIKE '%s' ", xtable);
	  strcat (sql, sql2);
	  sprintf (sql2, "AND f_geometry_column LIKE '%s'", xgeom);
	  strcat (sql2, sql);
      }
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return -1;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		f_table_name = results[(i * columns) + 0];
		f_geometry_column = results[(i * columns) + 1];
		if (!do_compute_layer_statistics
		    (sqlite, f_table_name, f_geometry_column,
		     SPATIALITE_STATISTICS_GENUINE))
		  {
		      error = 1;
		      break;
		  }
		else
		    count++;
	    }
      }
    sqlite3_free_table (results);
    if (error)
	return -1;
    return count;
}

static int
views_layer_statistics (sqlite3 * sqlite, const char *table, const char *column)
{
/* updating VIEWS_LAYER_STATISTICS metadata */
    int count = 0;
    char xtable[1024];
    char xgeom[1024];
    char sql[8192];
    char sql2[2048];
    int ret;
    const char *view_name;
    const char *view_geometry;
    int i;
    char **results;
    int rows;
    int columns;
    int error = 0;

    if (table == NULL && column == NULL)
      {
	  /* processing any table/geometry found in VIEWS_GEOMETRY_COLUMNS */
	  strcpy (sql, "SELECT view_name, view_geometry ");
	  strcat (sql, "FROM views_geometry_columns");
      }
    else if (column == NULL)
      {
	  /* processing any geometry belonging to this table */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (sql, "SELECT view_name, view_geometry ");
	  strcat (sql, "FROM views_geometry_columns ");
	  sprintf (sql2, "WHERE view_name LIKE '%s'", xtable);
	  strcat (sql, sql2);
      }
    else
      {
	  /* processing a single table/geometry entry */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (xgeom, column);
	  clean_sql_string (xgeom);
	  strcpy (sql, "SELECT view_name, view_geometry ");
	  strcat (sql, "FROM views_geometry_columns ");
	  sprintf (sql2, "WHERE view_name LIKE '%s' ", xtable);
	  strcat (sql, sql2);
	  sprintf (sql2, "AND view_geometry LIKE '%s'", xgeom);
	  strcat (sql2, sql);
      }
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return -1;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		view_name = results[(i * columns) + 0];
		view_geometry = results[(i * columns) + 1];
		if (!do_compute_layer_statistics
		    (sqlite, view_name, view_geometry,
		     SPATIALITE_STATISTICS_VIEWS))
		  {
		      error = 1;
		      break;
		  }
		else
		    count++;
	    }
      }
    sqlite3_free_table (results);
    if (error)
	return -1;
    return count;
}

static int
virts_layer_statistics (sqlite3 * sqlite, const char *table, const char *column)
{
/* updating VIRTS_LAYER_STATISTICS metadata */
    int count = 0;
    char xtable[1024];
    char xgeom[1024];
    char sql[8192];
    char sql2[2048];
    int ret;
    const char *f_table_name;
    const char *f_geometry_column;
    int i;
    char **results;
    int rows;
    int columns;
    int error = 0;

    if (table == NULL && column == NULL)
      {
	  /* processing any table/geometry found in GEOMETRY_COLUMNS */
	  strcpy (sql, "SELECT virt_name, virt_geometry ");
	  strcat (sql, "FROM virts_geometry_columns");
      }
    else if (column == NULL)
      {
	  /* processing any geometry belonging to this table */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (sql, "SELECT virt_name, virt_geometry ");
	  strcat (sql, "FROM virts_geometry_columns ");
	  sprintf (sql2, "WHERE virt_name LIKE '%s'", xtable);
	  strcat (sql, sql2);
      }
    else
      {
	  /* processing a single table/geometry entry */
	  strcpy (xtable, table);
	  clean_sql_string (xtable);
	  strcpy (xgeom, column);
	  clean_sql_string (xgeom);
	  strcpy (sql, "SELECT virt_name, virt_geometry ");
	  strcat (sql, "FROM virts_geometry_columns ");
	  sprintf (sql2, "WHERE virt_name LIKE '%s' ", xtable);
	  strcat (sql, sql2);
	  sprintf (sql2, "AND virt_geometry LIKE '%s'", xgeom);
	  strcat (sql2, sql);
      }
    ret = sqlite3_get_table (sqlite, sql, &results, &rows, &columns, NULL);
    if (ret != SQLITE_OK)
	return -1;
    if (rows < 1)
	;
    else
      {
	  for (i = 1; i <= rows; i++)
	    {
		f_table_name = results[(i * columns) + 0];
		f_geometry_column = results[(i * columns) + 1];
		if (!do_compute_layer_statistics
		    (sqlite, f_table_name, f_geometry_column,
		     SPATIALITE_STATISTICS_VIRTS))
		  {
		      error = 1;
		      break;
		  }
		else
		    count++;
	    }
      }
    sqlite3_free_table (results);
    if (error)
	return -1;
    return count;
}

SPATIALITE_DECLARE int
update_layer_statistics (sqlite3 * sqlite, const char *table,
			 const char *column)
{
/* updating LAYER_STATISTICS metadata [main] */
    int ret;
    int count = 0;
    ret = genuine_layer_statistics (sqlite, table, column);
    if (ret < 0)
	return 0;
    count += ret;
    ret = views_layer_statistics (sqlite, table, column);
    if (ret < 0)
	return 0;
    count += ret;
    ret = virts_layer_statistics (sqlite, table, column);
    if (ret < 0)
	return 0;
    count += ret;
    return count;
}
