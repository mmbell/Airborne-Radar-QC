/*
 *  DEM.h
 *
 *  Created by Michael Bell on 2/4/11.
 *  Based on example code from GeoTIFF library
 *  Copyright 2011. All rights reserved.
 *
 */

#ifndef DEM_H
#define DEM_H

#include "geotiff.h"
#include "geo_normalize.h"
#include "tiffio.h"

class DEM  
{

public:
	DEM();
	~DEM();
	
	bool readDem(char* fname);
	int getElevation(const double& lat, const double& lon);
	bool dumpAscii(int skip);
	int getMaxElevation();
	
private:
	int GTIFReportACorner( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
						const char * corner_name,
						double x, double y, int inv_flag, int dec_flag );
	void GTIFPrintCorners( GTIF *, GTIFDefn *, FILE *, int, int, int, int );

	char * demFilename;
	size_t npixels;
	uint32 xsize, ysize;
	int16* elevations;
	int dx, dy;
	double refLat, refLon;
	

};


#endif
