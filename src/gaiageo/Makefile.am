
AM_CFLAGS = -D_LARGE_FILE=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE=1 \
        -DNDEBUG=1 \
	@OMIT_GEOS_FLAGS@ @GEOSADVANCED_FLAGS@ @OMIT_BOOSTGEOMETRY_FLAGS@ @OMIT_PROJ_FLAGS@ \
	@OMIT_ICONV_FLAGS@ @OMIT_MATHSQL_FLAGS@ \
	@OMIT_EPSG_FLAGS@ @OMIT_GEOCALLBACKS_FLAGS@ @OMIT_FREEXL_FLAGS@
	
AM_CXXFLAGS = $(AM_CFLAGS)

INCLUDES = @CFLAGS@
INCLUDES += -I$(top_srcdir)/src/headers

noinst_LTLIBRARIES = libgaiageo.la

libgaiageo_la_SOURCES = gg_advanced.c \
	gg_endian.c \
	gg_geodesic.c \
	gg_geometries.c \
	gg_geoscvt_boost.cpp \
	gg_geoscvt.c \
	gg_relations_boost.cpp \
	gg_relations.c \
	gg_shape.c \
	gg_transform.c \
	gg_wkb.c \
	gg_wkt.c \
	gg_vanuatu.c \
	gg_ewkt.c \
	gg_geoJSON.c \
	gg_kml.c \
	gg_gml.c

MOSTLYCLEANFILES = *.gcna *.gcno *.gcda
