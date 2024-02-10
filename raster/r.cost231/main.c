
/****************************************************************************
 *
 * MODULE:       r.cost231
 * AUTHOR(S):    Andrej Hrovat, Jozef Stefan Institute                
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Calculates radio coverage from a single base station 
 *               according to Cost231 model
 *             
 *
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:

23-sep-2013 (I.O.)
 - minor modification

30-sep-2013 (I.O.)
 - corrected half-pixel shift in the output map

1-aug-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
     option name changed: inputDEM -> input_dem
 - minor corrections (to remove compiler warnings, redundant code, etc.)
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>


/*
 * Calculation of Cost231 for a pair of points:
 * -----------------------------------------
 * tr_height_eff ... effective height of Tx: difference between total Tx and total Rx height
 * d ... distance between Rx and Tx
 * freq ... radio frequency
 * rec_height ... AGL height of Rx
 * limit ... distance up to which cost231 should be calculated (input parameter of the model)
 * area_type ... type of area (urban, suburban, open)
 */
FCELL calc_cost231( double tr_height_eff, double d, double freq, double rec_height, double limit, char *area_type)
{ 
  double Lm;   /* Path loss in metropolitan area (in dB) */
  double Lusu;  /* Path loss in medium cities and suburban area (in dB) */
  double ahr;  /* Correction factor */
  FCELL x;     /* return value */

  /* get absolute value of effective height */
  tr_height_eff = fabs( tr_height_eff);

  d = d / 1000;  /* in cost231, distance has to be given in km */    

  /* If Rx and Tx are closer than 10m, then do not calculate */
  if (d < 0.01 || d > limit)
    {    
      Rast_set_f_null_value( &x, 1);
      return x;
    }

  /* Correction factor ahr */
  ahr = ( 1.1 * log10( freq) - 0.7) * rec_height - ( 1.56 * log10( freq) - 0.8);

  /* Lm */
  Lm = 46.33 + 33.9 * log10( freq) - 13.82 * log10( tr_height_eff) - ahr + ( 44.9 - 6.55 * log10( tr_height_eff)) * log10( d) + 3;

  /* Lusu */
  Lusu = 46.33 + 33.9 * log10( freq) - 13.82 * log10( tr_height_eff) - ahr + ( 44.9-6.55 * log10( tr_height_eff)) * log10( d);
  
  /* return x */
  if ( strcmp( area_type, "metropolitan") == 0)       x = (FCELL)Lm;
  else if ( strcmp( area_type, "medium_cities") == 0) x = (FCELL)Lusu;
  else G_fatal_error( _("Unknown area type: [%s]."), area_type);
  
  return x;
}
 

/*
 * main function
 */
int main( int argc, char *argv[])
{
  double east;
  double north;
  double ant_height, frequency, radius; //, dem_height;
  double rec_height = 1.5; /* height of receiver from the ground */
    
  struct Cell_head window;  /*   database window         */
  struct Cell_head cellhd;  /* it stores region information,
                               and header information of rasters */
  char *name;               /* input raster name */
  char *result;             /* output raster name */
  const char *mapset;       /* mapset name */
  void *inrast;             /* input buffer */
  unsigned char *outrast;   /* output buffer */
  int nrows, ncols;
  int row, col;
  int tr_row, tr_col;
  int infd, outfd;          /* file descriptor */
  int verbose;
  struct History history;   /* holds meta-data (title, comments,..) */

  struct GModule *module;   /* GRASS module for parsing arguments */

  struct Option *input,*opt1, *opt2, *opt4, *opt5, *opt3,  *output;  /* options */
  struct Flag *flag1;  /* flags */


  /* initialize GIS environment */
  G_gisinit( argv[ 0]);  /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword( _("raster"));
  G_add_keyword( _("cost231"));
  module->description = _("RaPlaT - Cost231 module (v01aug2017)");

  /* Define the different options as defined in gis.h */
  input = G_define_standard_option(G_OPT_R_INPUT);  
  input->key = "input_dem";      
  output = G_define_standard_option(G_OPT_R_OUTPUT);

  /* Define the different flags */
  flag1 = G_define_flag();
  flag1->key = 'q';
  flag1->description = _("Quiet");

  opt1 = G_define_option();
  opt1->key = "coordinate";
  opt1->type = TYPE_STRING;
  opt1->required = YES;
  opt1->key_desc = "x,y";
  opt1->description = _("Base station coordinates");

  opt2 = G_define_option();
  opt2->key = "ant_height";
  opt2->type = TYPE_DOUBLE;
  opt2->required = NO;
  opt2->answer = "10";
  opt2->description = _("Transmitter antenna height [m]");

  opt4 = G_define_option();
  opt4->key = "radius";
  opt4->type = TYPE_DOUBLE;
  opt4->required = NO;
  opt4->answer = "10";
  opt4->description = _("Computation radius [km]");

  opt5 = G_define_option();
  opt5->key = "area_type";
  opt5->type = TYPE_STRING;
  opt5->required = NO;
  opt5->description = _("Area type");
  opt5->options = "metropolitan,medium_cities";
  opt5->answer = "medium_cities";

  opt3 = G_define_option();
  opt3->key = "frequency";
  opt3->type = TYPE_DOUBLE;
  opt3->required = YES;
  opt3->description = _("Frequency [MHz]");

  /* options and flags parser */
  if (G_parser( argc, argv)) exit( EXIT_FAILURE);

  /* stores options and flags to variables */
  name = input->answer;
  result = output->answer;
  verbose = !flag1->answer;
  G_scan_easting( opt1->answers[ 0], &east, G_projection());
  G_scan_northing( opt1->answers[ 1], &north, G_projection());
  sscanf(opt2->answer, "%lf", &ant_height);
  sscanf(opt4->answer, "%lf", &radius);
  sscanf(opt3->answer, "%lf", &frequency);  


  /* returns NULL if the map was not found in any mapset, mapset name otherwise */
  mapset = G_find_raster(name, "");
  if ( mapset == NULL)
    G_fatal_error( _("Raster map <%s> not found"), name);

  if ( G_legal_filename( result) < 0)
    G_fatal_error( _("<%s> is an illegal file name"), result);


  /* Rast_open_old - returns file destriptor (>0) */
  if ( ( infd = Rast_open_old( name, mapset)) < 0)
    G_fatal_error( _("Unable to open raster map <%s>"), name);


  /* open input raster */
  Rast_get_cellhd( name, mapset, &cellhd);

  G_debug( 3, "number of rows %d", cellhd.rows);

  G_get_window ( &window);

  /* Allocate input buffer */
  inrast = Rast_allocate_buf( FCELL_TYPE);

  /* Allocate output buffer, use input map data_type */
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf( FCELL_TYPE);

  /* controlling, if we can write the raster */
  if ( ( outfd = Rast_open_new( result, FCELL_TYPE)) < 0)
    G_fatal_error( _("Unable to create raster map <%s>"), result);

  /* check if specified transmitter location inside window   */
  if ( east < window.west || east > window.east ||
       north > window.north || north < window.south)
    G_fatal_error(_("Specified base station  coordinates are outside current region bounds."));  

  /* map array coordinates for transmitter */
  tr_row = ( window.north - north) / window.ns_res;
  tr_col = ( east - window.west) / window.ew_res;

  /* total height of transmitter */
  FCELL trans_elev;
  double trans_total_height;
  Rast_get_row( infd, inrast, tr_row, FCELL_TYPE);
  trans_elev = ( (FCELL *) inrast)[tr_col];
  trans_total_height = (double)trans_elev + ant_height;
    
  // check if transmitter is on DEM
  if ( isnan( (double)trans_elev))
    G_fatal_error( _("Transmitter outside raster DEM map."));
    
  /* do COST231 */
  double height_diff_Tx_Rx; /* difference of total heights (Tx - Rx) */
  double dist_Tx_Rx; /* distance between receiver and transmitter */
  double rec_east, rec_north;  /* receiver coordinates */

  /* for each row */
  for ( row = 0; row < nrows; row++) 
  {  
    if ( verbose) G_percent( row, nrows, 2);

    FCELL f_in, f_out;
 
    /* read input map */
    Rast_get_row( infd, inrast, row, FCELL_TYPE);
  
    /* process the data */
    for ( col = 0; col < ncols; col++) 
    { 
      f_in = ( (FCELL *)inrast)[col];

      /* calculate receiver coordinates */
      rec_east = window.west + window.ew_res / 2.0 + ( col * window.ew_res);
      rec_north = window.north - window.ns_res / 2.0 - ( row * window.ns_res);


      /* calculate distance */
      dist_Tx_Rx = sqrt( pow( ( east - rec_east), 2) + pow( ( north - rec_north), 2));

      /* calculate height difference */
      if ( (double)trans_elev > (double)f_in)
        height_diff_Tx_Rx = trans_total_height - (double)f_in - rec_height;
      else
        height_diff_Tx_Rx = ant_height;

      /* calculate cost231 */
      f_out = calc_cost231( height_diff_Tx_Rx, dist_Tx_Rx, frequency, rec_height, radius, opt5->answer);

      ( (FCELL *)outrast)[ col] = f_out;
    }
      
    /* write raster row to output raster map */
    Rast_put_row( outfd, outrast, FCELL_TYPE);
  }
         

  /* memory cleanup */
  G_free( inrast);
  G_free( outrast);

  /* closing raster maps */
  Rast_close( infd);
  Rast_close( outfd);

  /* add command line incantation to history file */
  Rast_short_history( result, "raster", &history);
  Rast_command_history( &history);
  Rast_write_history( result, &history);


  exit( EXIT_SUCCESS);
}
