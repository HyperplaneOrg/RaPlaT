
/****************************************************************************
 *
 * MODULE:       r.hata
 * AUTHOR(S):    Andrej Vilhar (original version), Jozef Stefan Institute
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Calculates radio coverage from a single base station 
 *               according to Hata model
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

02-jun-2014 (I.O.)
 - additional parameter rx_ant_height (default = 1.5)
 - "inverse" mode added (RX and TX roles exchanged, needed for transmitter localization)
 - C source style cleaned-up

20-jul-2017 (I.O.)
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
 * Calculation of Hata for a pair of points:
 * -----------------------------------------
 * tr_height_eff ... effective height of Tx: difference between total Tx and total Rx height
 * d ... distance between Rx and Tx
 * freq ... radio frequency
 * rec_height ... AGL height of Rx
 * limit ... distance up to which hata should be calculated (input parameter of the model)
 * area_type ... type of area (urban, suburban, open)
 */
FCELL calc_hata(double tr_height_eff, double d, double freq, double rec_height, double limit, char *area_type)
{ 
  double Lu;   /* Path loss in urban area (in dB) */
  double Lo;   /* Path loss in open area (in dB) */
  double Lsu;  /* Path loss in suburban area (in dB) */
  double ahr;  /* Correction factor */
  FCELL x;     /* return value */

  /* get absolute value of effective height */
  tr_height_eff = fabs( tr_height_eff);

  d = d/1000;  /* in Hata, distance has to be given in km */    

  /* If Rx and Tx are closer than 10m, then do not calculate */
  if ( d < 0.01 || d > limit)
  {    
    Rast_set_f_null_value( &x, 1);
    return x;
  }

  /* Correction factor ahr */
  ahr = ( 1.1 * log10( freq) - 0.7) * rec_height - (1.56 * log10( freq) - 0.8);

  /* Lu */
  Lu = 69.55 + 26.16 * log10( freq) - 13.82 * log10( tr_height_eff) - ahr + ( 44.9 - 6.55 * log10( tr_height_eff)) * log10( d);

  /* Lo */
  Lo = Lu - 4.78 * pow( log10( freq), 2) + 18.33 * log10(freq) - 40.94;

  /* Lsu */
  Lsu = Lu - 2 * pow( log10( freq / 28), 2) - 5.4;

  /* return x */
  if ( strcmp( area_type, "urban") == 0)
    x = (FCELL) Lu;
  else if ( strcmp( area_type, "suburban") == 0)
    x = (FCELL) Lsu;
  else if ( strcmp( area_type, "open") == 0)
    x = (FCELL) Lo;
  else
    G_fatal_error( _("Unknown area type: [%s]."), area_type);

  return x;
}
 

/*
 * main function
 */
int main( int argc, char *argv[])
{
  double east;
  double north;
  double ant_height, frequency, radius;
  double rec_ant_height;

  struct Cell_head window;  /*  database window         */
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
  int inverse_mode_f;

  struct History history;   /* holds meta-data (title, comments,..) */

  struct GModule *module;   /* GRASS module for parsing arguments */

  struct Option *input, *opt1, *opt2, *opt3, *opt4, *opt5, *opt6, *output;  /* options */
  struct Flag *flag1, *flag2;  /* flags */

  /* initialize GIS environment */
  G_gisinit( argv[ 0]);  /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("hata"));
  module->description = _("RaPlaT - Hata module (v20jul2017)");

  /* Define the different options as defined in gis.h */
  input = G_define_standard_option( G_OPT_R_INPUT);
  input->key = "input_dem";
  output = G_define_standard_option( G_OPT_R_OUTPUT);

  /* Define the different flags */
  flag1 = G_define_flag();
  flag1->key = 'q';
  flag1->description = _("Quiet");

  flag2 = G_define_flag();
  flag2->key = 'i';
  flag2->description = _("Inverse mode (RX and TX roles exchanged)");

  opt5 = G_define_option();
  opt5->key = "area_type";
  opt5->type = TYPE_STRING;
  opt5->required = NO;
  opt5->description = _("Area type");
  opt5->options = "urban,suburban,open";
  opt5->answer = "urban";

  opt1 = G_define_option();
  opt1->key = "coordinate";
  opt1->type = TYPE_STRING;
  opt1->required = YES;
  opt1->key_desc = "x,y";
  opt1->description = _("Base station coordinates, or receiver location in inverse mode");

  opt4 = G_define_option();
  opt4->key = "radius";
  opt4->type = TYPE_DOUBLE;
  opt4->required = NO;
  opt4->answer = "10";
  opt4->description = _("Computation radius [km]");

  opt2 = G_define_option();
  opt2->key = "ant_height";
  opt2->type = TYPE_DOUBLE;
  opt2->required = NO;
  opt2->answer = "10";
  opt2->description = _("Transmitter antenna height [m]");

  opt6 = G_define_option();
  opt6->key = "rx_ant_height";
  opt6->type = TYPE_DOUBLE;
  opt6->required = NO;
  opt6->answer = "1.5";
  opt6->description = _("Receiver antenna height [m]");

  opt3 = G_define_option();
  opt3->key = "frequency";
  opt3->type = TYPE_DOUBLE;
  opt3->required = YES;
  opt3->description = _("Frequency (MHz)");

   /* options and flags parser */
  if ( G_parser(argc, argv))
    exit( EXIT_FAILURE);

   /* stores options and flags to variables */
  name = input->answer;
  result = output->answer;
  verbose = !flag1->answer;
  inverse_mode_f = flag2->answer;
  G_scan_easting( opt1->answers[0], &east, G_projection());
  G_scan_northing( opt1->answers[1], &north, G_projection());
  sscanf( opt2->answer, "%lf", &ant_height);
  sscanf( opt4->answer, "%lf", &radius);
  sscanf( opt3->answer, "%lf", &frequency);

  sscanf( opt6->answer, "%lf", &rec_ant_height);


  /* returns NULL if the map was not found in any mapset, 
   * mapset name otherwise */
  mapset = G_find_raster(name, "");
  if ( mapset == NULL)
    G_fatal_error( _("Raster map <%s> not found"), name);

  if ( G_legal_filename(result) < 0)
    G_fatal_error( _("<%s> is an illegal file name"), result);


  /* Rast_open_old - returns file destriptor (>0) */
  if ( ( infd = Rast_open_old( name, mapset)) < 0)
    G_fatal_error( _("Unable to open raster map <%s>"), name);


  /* open input raster */
  Rast_get_cellhd( name, mapset, &cellhd);

  G_debug( 3, "number of rows %d", cellhd.rows);

  G_get_window( &window);

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
  {
    if ( !inverse_mode_f)
      G_fatal_error( _("Specified base station coordinates are outside current region bounds."));
    else
      G_fatal_error( _("Specified receive location coordinates are outside current region bounds."));
  }
    

  /* map array coordinates for transmitter */
  tr_row = ( window.north - north) / window.ns_res;
  tr_col = ( east - window.west) / window.ew_res;


  /* total height of transmitter */
  FCELL trans_elev;
  double trans_total_height;
  Rast_get_row( infd, inrast, tr_row, FCELL_TYPE);
  trans_elev = ( (FCELL *) inrast)[ tr_col];

  //G_message( _("Transmiter elevatino [%f]"), trans_elev);

  trans_total_height = (double) trans_elev + ant_height;

  // check if transmitter is on DEM
  if ( isnan( (double) trans_elev))
  {
    if ( !inverse_mode_f)
      G_fatal_error( _("Transmitter outside raster DEM map."));
    else
      G_fatal_error( _("Receiver outside raster DEM map."));
  }
    
  /* do HATA */

  double height_diff_Tx_Rx; /* difference of total heights (Tx - Rx) */
  double dist_Tx_Rx; /* distance between receiver and transmitter */
  double rec_east, rec_north;  /* receiver coordinates */
   

  /* for each row */
  for ( row = 0; row < nrows; row++) 
  {

    if ( verbose)
      G_percent( row, nrows, 10);

    FCELL f_in, f_out;

    /* read input map */
    Rast_get_row( infd, inrast, row, FCELL_TYPE);
  
    /* process the data */
    for ( col = 0; col < ncols; col++) 
    { 
      f_in = ( (FCELL *) inrast)[ col];
  
      /* calculate receiver coordinates */
      rec_east = window.west + window.ew_res / 2.0 + ( col * window.ew_res);
      rec_north = window.north - window.ns_res / 2.0 - ( row * window.ns_res);


      /* calculate distance */
      // Inverse mode: distance components change sign, but this can be ignored
      dist_Tx_Rx = sqrt( pow( ( east - rec_east), 2) + pow( ( north - rec_north), 2));


      if ( !inverse_mode_f)
      {
       /* calculate height difference */
        if ( (double) trans_elev > (double) f_in)
          height_diff_Tx_Rx = trans_total_height - (double) f_in;
        else
          height_diff_Tx_Rx = ant_height;
      }
      else
      {
        // inverse mode - setup height_diff_Tx_Rx properly
        // (trans_elev is actually receiver elevation)
        double trans_elev_imode = (double) f_in;
        double rec_elev_imode = (double) trans_elev;

        trans_total_height = trans_elev_imode + ant_height;
        if ( trans_elev_imode > rec_elev_imode)
          height_diff_Tx_Rx = trans_total_height - rec_elev_imode;
        else
          height_diff_Tx_Rx = ant_height;
      }

      /* calculate hata */
      f_out = calc_hata( height_diff_Tx_Rx, dist_Tx_Rx, frequency, rec_ant_height, radius, opt5->answer);

      ( (FCELL *) outrast)[ col] = f_out;
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

