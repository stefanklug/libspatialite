#include "gg_geometry_boost.h"

#ifdef SPL_AMALGAMATION		/* spatialite-amalgamation */
#include <spatialite/sqlite3ext.h>
#else
#include <sqlite3ext.h>
#endif

#include <spatialite/gaiageo.h>

gaiaGeomCollPtr
fromBoostGeometry (const geometry_collection& geom);

geometry_collection
toBoostGeometry (const gaiaGeomCollPtr gaia);
