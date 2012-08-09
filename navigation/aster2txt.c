/*
 * aster2txt.c -- Conevert ASTER data to ASCII format for Airborne Radar QC navigation code
 * Michael M. Bell
 * Based on listgeo.c by Niles D. Ritter
 *
 */

#include "geotiff.h"
#include "xtiffio.h"
#include "geo_normalize.h"
#include "geo_simpletags.h"
#include "geovalues.h"
#include "tiffio.h"
#include "cpl_serv.h"
#include <stdio.h>
#include <string.h>

static void GTIFPrintCorners( GTIF *, GTIFDefn *, FILE *, int, int, int, int );
void Usage()

{
    printf( 
        "%s", 
        "Usage: aster2txt <filename>\n"
        "\n"
        "  filename: Name of the GeoTIFF file to dump.\n" );
        
    exit( 1 );
}

int main(int argc, char *argv[])
{
    char	*fname = NULL;
    char        *outfile = NULL;
    TIFF 	*tif=(TIFF*)0;  /* TIFF-level descriptor */
    GTIF	*gtif=(GTIF*)0; /* GeoKey-level descriptor */
    int		i, j, norm_print_flag = 0, proj4_print_flag = 0;
    int		tfw_flag = 0, inv_flag = 0, dec_flag = 1;
    int         st_test_flag = 0;
    size_t npixels;
    int16* elevations;
    int x, y, skip;
    GTIFDefn	defn;
    FILE* out;

    /*
     * Handle command line options.
     */
    for( i = 1; i < argc; i++ )
    {
        if( fname == NULL && argv[i][0] != '-' )
            fname = argv[i];
        else
        {
            Usage();
        }
    }
    if( fname == NULL)
        Usage();

    skip = 5;

    /*
     * Open the file, read the GeoTIFF information, and print to stdout. 
     */
    int flength = strlen(fname);
    outfile = (char *) malloc(flength);
    strncpy(outfile,fname,flength-4);
    strcat(outfile, ".asc\0");
    out = fopen(outfile, "w");
    printf("Writing to %s\n\n",outfile);
    free(outfile);

    tif=XTIFFOpen(fname,"r");
    if (!tif) goto failure;
    
    gtif = GTIFNew(tif);
    if (!gtif)
      {
	fprintf(stderr,"failed in GTIFNew\n");
	goto failure;
      }
	
    /* dump the GeoTIFF metadata to std out */

    GTIFPrint(gtif,0,0);

    if( GTIFGetDefn( gtif, &defn ) )
      {
	uint32		xsize, ysize;
        
	printf( "\n" );
	GTIFPrintDefn( &defn, stdout );
	
	if( proj4_print_flag )
	  {
	    printf( "\n" );
	    printf( "PROJ.4 Definition: %s\n", GTIFGetProj4Defn(&defn));
	  }
	int count, orient;
	double* data;
	TIFFGetField( tif, TIFFTAG_IMAGEWIDTH, &xsize );
	TIFFGetField( tif, TIFFTAG_IMAGELENGTH, &ysize );
	TIFFGetField( tif, TIFFTAG_GEOPIXELSCALE, &count, &data);
	TIFFGetField( tif, TIFFTAG_ORIENTATION, &orient );
	printf("Orientation:%d\n",orient);

	GTIFPrintCorners( gtif, &defn, stdout, xsize, ysize, inv_flag, dec_flag );
	const char* project = "ASTERGDEM";
	const char* yymmdd = "090629";
	double originx, originy;
	GTIFImageToPCS( gtif, &originx, &originy);
	int lat = (int)originy*1000;
	int lon = (int)originx*1000;
	int xmin = 0;
	int ymin = 0;
	int nx = xsize/skip + 1;
	int ny = ysize/skip + 1;
	int dx = 30 * skip;
	int dy = 30 * skip;
	int16 elev = 0;
	fprintf(out, "%12s%12s%7d%7d%7d%7d%7d%7d%7d%7d\n",
	       project, yymmdd, lat, lon, xmin, ymin, nx, ny, dx,dy);
	npixels = xsize * ysize;
	int16* buf;
	tsample_t sample;
	int nbytes;
	nbytes = TIFFScanlineSize(tif);
	buf = (int16*) _TIFFmalloc(nbytes);
	elevations = (int16*) _TIFFmalloc(npixels * sizeof (int16));
	if (elevations != NULL) {
	    for( y = 0; y < ysize; y++ )
	      {
		uint32 row = ysize - 1 - y;
		TIFFReadScanline(tif, buf, row, sample);
		for( x = 0; x < xsize; x++ ) 
		  {
		    elevations[y*xsize + x] = buf[x];
		  }
	      }	    
	  _TIFFfree(buf);
	}
	for( y = 0; y < ysize; y+=skip ) {
	  for( x = 0; x < xsize; x+=skip ) {
	    fprintf(out, "%6d", elevations[y*xsize + x]);
	  }
	  fprintf(out, "\n");
	}
	_TIFFfree(elevations);

      }
    
  Success:
    GTIFFree(gtif);
    if( st_test_flag )
        ST_Destroy( (ST_TIFF *) tif );
    else
        XTIFFClose(tif);
    GTIFDeaccessCSV();
    return 0;
		
  failure:
    fprintf(stderr,"failure in listgeo\n");
    if (tif) XTIFFClose(tif);
    if (gtif) GTIFFree(gtif);
    GTIFDeaccessCSV();
    return 1;
}

static int GTIFReportACorner( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
                              const char * corner_name,
                              double x, double y, int inv_flag, int dec_flag )

{
    double	x_saved, y_saved;

    /* Try to transform the coordinate into PCS space */
    if( !GTIFImageToPCS( gtif, &x, &y ) )
        return FALSE;
    
    x_saved = x;
    y_saved = y;

    fprintf( fp_out, "%-13s ", corner_name );

    if( defn->Model == ModelTypeGeographic )
    {
	if (dec_flag) 
	{
	    fprintf( fp_out, "(%.7f,", x );
	    fprintf( fp_out, "%.7f)\n", y );
	} 
	else 
	{
	    fprintf( fp_out, "(%s,", GTIFDecToDMS( x, "Long", 2 ) );
	    fprintf( fp_out, "%s)\n", GTIFDecToDMS( y, "Lat", 2 ) );
	}
    }
    else
    {
        fprintf( fp_out, "(%12.3f,%12.3f)", x, y );

        if( GTIFProj4ToLatLong( defn, 1, &x, &y ) )
        {
	    if (dec_flag) 
	    {
                fprintf( fp_out, "  (%.7f,", x );
                fprintf( fp_out, "%.7f)", y );
	    } 
	    else 
	    {
		fprintf( fp_out, "  (%s,", GTIFDecToDMS( x, "Long", 2 ) );
		fprintf( fp_out, "%s)", GTIFDecToDMS( y, "Lat", 2 ) );
	    }
        }

        fprintf( fp_out, "\n" );
    }

    if( inv_flag && GTIFPCSToImage( gtif, &x_saved, &y_saved ) )
    {
        fprintf( fp_out, "      inverse (%11.3f,%11.3f)\n", x_saved, y_saved );
    }
    
    return TRUE;
}

static void GTIFPrintCorners( GTIF *gtif, GTIFDefn *defn, FILE * fp_out,
                              int xsize, int ysize, int inv_flag, int dec_flag )

{
    printf( "\nCorner Coordinates:\n" );
    if( !GTIFReportACorner( gtif, defn, fp_out,
                            "Upper Left", 0.0, 0.0, inv_flag, dec_flag ) )
    {
        printf( " ... unable to transform points between pixel/line and PCS space\n" );
        return;
    }

    GTIFReportACorner( gtif, defn, fp_out, "Lower Left", 0.0, ysize, 
                       inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Upper Right", xsize, 0.0,
                       inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Lower Right", xsize, ysize,
                       inv_flag, dec_flag );
    GTIFReportACorner( gtif, defn, fp_out, "Center", xsize/2.0, ysize/2.0,
                       inv_flag, dec_flag );
}
