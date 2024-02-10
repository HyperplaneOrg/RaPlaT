/****************************************************************************
 *
 * MODULE:       r.sector
 * AUTHOR(S):    Andrej Vilhar, Jozef Stefan Institute
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute             
 *
 * PURPOSE:      Takes propagation pathloss and calculates additional
 *               gain/pathloss according to antenna's directional diagram.
 *                            
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:
4-jan-2012:
 - corrections of angle calculations (by A. Vilhar)

5-feb-2013 (Igor Ozimek):
 - added support for dBd/dBi in the .MPI file GAIN parameter (default is dBd)
 - diagnose error if GAIN parameter missing in the .MPI file
 - skip calculations of raster points with undefined height (DEM)
   (set to NULL instead, avoiding module crash with segmentation fault)
 - some "housekeeping" (cleaning and reformatting of source code)

23-sep-2013 (I.O.)
 - minor modification

2-oct-2013 (I.O.)
 - corrected half-pixel shift in the output map
 - corrected 2x divide by zero (if d_north == 0, if dist_Tx_Rx == 0))

15-july-2014 (I.O.)
 - additional parameter rx_ant_height (default = 1.5)

28-jul-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
     option name changed: inputDEM -> input_dem
 - minor corrections (to remove compiler warnings, redundant code, etc.)
 - correction of a computation error giving undefined results (e.g. inf or null in computed raster map); caused by vert_diag_angle possibly >360.0 (and consequently temp_angle>360) after mechanical tilt applied (positiv or negative)
 - TODO:
     beam direction (input parameter) MAY NOT BE negative (it is not checked!)

6-dec-2018
 - three minor changes (workarounds) to avoid C compiler warnings (seems to be a compiler bug)
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>

#define PI 3.14159265


/*
 * main function
 */
int main(int argc, char *argv[])
{
  /* antenna data */
  int beamDirection;
  int mechanicalAntennaTilt;
  double heightAGL;
  double east;
  double north;
  
  //char antennaFolder[80];
  char antdata_file[100];

  double totalHeight;
  double gain;
  double horizontal[360];
  double vertical[360];
  /* end of antenna data */

  double rec_height;
  double radius;

  FILE *in;
  char buffer[256], text[32];

  struct Cell_head window; /* database window         */
  struct Cell_head cellhd; /* it stores region information and header information of rasters */

  char *name, *name2;     /* input raster name */
  char *result;           /* output raster name */
  const char *mapset, *mapset2; /* mapset name */
  void *inrast, *inrast2; /* input buffer */
  unsigned char *outrast; /* output buffer */
  int nrows, ncols;
  int row, col, j;
  int tr_row, tr_col;
  int infd, infd2, outfd; /* file descriptor */
  int verbose;
  struct History history; /* holds meta-data (title, comments,..) */

  struct GModule *module; /* GRASS module for parsing arguments */

  struct Option *input, *input2, *opt1, *opt2, *opt3, *opt4, *opt5, *opt6, *opt7, *opt9, *output; /* options */
  struct Flag *flag1; /* flags */

  //int dem_defined = 0;

  /* default path to antenna diagram */
  char def_path[1000];
  strcpy(def_path,getenv("GISBASE"));
  strcat(def_path,"/etc/radio_coverage/antenna_diagrams/");


  /* initialize GIS environment */
  G_gisinit(argv[0]);  /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("directional diagram"));
  G_add_keyword(_("antenna sector"));
  module->description = _("RaPlaT - Sector module (v06dec2018)");
  
  /* Define the different options as defined in gis.h */
  input = G_define_standard_option(G_OPT_R_INPUT);
  input->key = "pathloss_raster";
  input->description = _("Name of isotropic antenna pathloss raster map");

  input2 = G_define_standard_option(G_OPT_R_INPUT);
  input2->key = "input_dem";
  input2->required = YES;
  input2->description = _("Name of elevation model raster map - required for transmitter height determination");
      
  output = G_define_standard_option(G_OPT_R_OUTPUT);
  
  /* Define the different flags */
  flag1 = G_define_flag();
  flag1->key = 'q';
  flag1->description = _("Quiet");

  opt6 = G_define_option();
  opt6->key = "east";
  opt6->type = TYPE_DOUBLE;
  opt6->required = YES;
  opt6->label = _("Easting coordinate");    

  opt7 = G_define_option();
  opt7->key = "north";
  opt7->type = TYPE_DOUBLE;
  opt7->required = YES;
  opt7->label = _("Northing coordinate"); 
  
  opt9 = G_define_option();
  opt9->key = "radius";
  opt9->type = TYPE_DOUBLE;
  opt9->required = NO;
  opt9->answer = "10";
  opt9->description = _("Computation radius [km]");

  opt2 = G_define_option();
  opt2->key = "ant_data_file";
  opt2->type = TYPE_STRING;
  opt2->required = YES;
  opt2->label = _("Antenna data file");	

  opt4 = G_define_option();
  opt4->key = "height_agl";
  opt4->type = TYPE_DOUBLE;
  opt4->required = YES;
  opt4->label = _("Above ground level height [m]");

  opt1 = G_define_option();
  opt1->key = "beam_direction";
  opt1->type = TYPE_INTEGER;
  opt1->required = YES;
  opt1->label = _("Beam direction [deg]");

  opt3 = G_define_option();
  opt3->key = "mech_tilt";
  opt3->type = TYPE_INTEGER;
  opt3->required = YES;
  opt3->label = _("Mechanical antenna tilt [deg]");

  opt5 = G_define_option();
  opt5->key = "rx_ant_height";
  opt5->type = TYPE_DOUBLE;
  opt5->required = NO;
  opt5->answer = "1.5";
  opt5->label = _("Receiver antenna height [m]");

    
  /* options and flags parser */
  if (G_parser(argc, argv))
    exit(EXIT_FAILURE);

  /* stores options and flags to variables */
  name = input->answer;
  name2 = input2->answer;
  result = output->answer;
  verbose = (!flag1->answer);
  sscanf( opt1->answer, "%d", &beamDirection);
  sscanf( opt2->answer, "%s", antdata_file);
  sscanf( opt3->answer, "%d", &mechanicalAntennaTilt);
  sscanf( opt4->answer, "%lf", &heightAGL);
  sscanf( opt5->answer, "%lf", &rec_height);
  sscanf( opt6->answer, "%lf", &east);
  sscanf( opt7->answer, "%lf", &north);  
  sscanf( opt9->answer, "%lf", &radius);


  /* returns NULL if the map was not found in any mapset, mapset name otherwise */
  mapset = G_find_raster(name, "");
  if (mapset == NULL)
    G_fatal_error(_("Raster pathloss map <%s> not found"), name);

  mapset2 = G_find_raster(name2, "");
  if (mapset2 == NULL)
    G_fatal_error(_("Raster map <%s> not found"), name2);

  if (G_legal_filename(result) < 0)
    G_fatal_error(_("<%s> is an illegal file name"), result);


  /* Rast_open_old - returns file destriptor (>0) */
  if ((infd = Rast_open_old(name, mapset)) < 0)
    G_fatal_error(_("Unable to open raster map <%s>"), name);

  if ((infd2 = Rast_open_old(name2, mapset2)) < 0)
    G_fatal_error(_("Unable to open raster map <%s>"), name2);


  /* open input raster */
  Rast_get_cellhd(name, mapset, &cellhd);

  G_debug(3, "number of rows %d", cellhd.rows);

  G_get_window (&window);


  /* Allocate input buffer */
  inrast = Rast_allocate_buf(FCELL_TYPE);
  inrast2 = Rast_allocate_buf(FCELL_TYPE);

  /* Allocate output buffer, use input map data_type */
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf(FCELL_TYPE);

  /* controlling, if we can write the raster */
  if ((outfd = Rast_open_new(result, FCELL_TYPE)) < 0)
    G_fatal_error(_("Unable to create raster map <%s>"), result);


  /* (1) calculate total height of the antenna */

  /* check if specified transmitter location inside window   */
  if (east < window.west || east > window.east || north > window.north || north < window.south)
    G_fatal_error(_("r.sector - specified base station  coordinates are outside current region bounds."));
  
  /* map array coordinates for transmitter */
  tr_row = (window.north - north) / window.ns_res;
  tr_col = (east - window.west) / window.ew_res;

  /* total height of transmitter */
  FCELL trans_elev;
  Rast_get_row(infd2, inrast2, tr_row, FCELL_TYPE);
  trans_elev = ((FCELL *) inrast2)[tr_col];
  totalHeight = (double)trans_elev + heightAGL;

  // check if transmitter is on DEM
  if ( isnan((double)trans_elev))							
    G_fatal_error(_("r.sector - transmitter outside raster DEM map."));
	

  /* (2) get antenna's gain and directional diagrams */
  char fileName[1000];

  /* set correct filename */

  strcpy (fileName, antdata_file);

  if ( strncmp(fileName, "/", 1) != 0)
  {
    strcat(def_path, fileName);
    strcpy(fileName, def_path);
  }	

  /* open file */
  if( (in = fopen(fileName,"r")) == NULL )
    G_fatal_error(_("Unable to open antenna diagram in file <%s>"), fileName);  

  /* get gain and find the beginning of horizontal diagram */
  /* gain can be specified as dBd (default) or dBi         */
  double temp_gain;
  char dBstr[32];
  gain = 9999;
  while (1)
  {
    if (!fgets (buffer, 250, in)) 
      G_fatal_error(_("Empty or corrupted antenna diagram file <%s>"), fileName); 

    int ret_scan = sscanf (buffer, "%s %lf %s", text, &temp_gain, dBstr);
    if (strcmp (text, "GAIN") == 0)
    {
      gain = temp_gain;
      if ( ret_scan == 2) gain += 2.15; /* GAIN is dBd by default */
      else
      {
        if (ret_scan != 3) G_fatal_error(_("Bad GAIN parameter in .MSI file"));
        if (strcmp (dBstr, "dBd") == 0) gain += 2.15; /* GAIN is dbd */
        else if (strcmp (dBstr, "dBi")) G_fatal_error(_("Bad GAIN parameter in .MSI file - should be dBd or dBi"));
      }
    }
    if (strcmp (text, "HORIZONTAL") == 0) break; /* we have reached the beggining of HOR data*/
  }
  if ( gain == 9999) G_fatal_error(_("Missing GAIN parameter in .MSI file"));

  double angle, loss;
  /* read horizontal data - one angle degree per step */
  for (j = 0; j < 360; j++)
  {
//    (void) fgets (buffer, 250, in);
    if( fgets (buffer, 250, in)); // (void) not working, using if instead (workaround)
    sscanf (buffer, "%lf %lf", &angle, &loss);
    if (j != (int)angle) G_fatal_error(_("Bad antenna diagram format")); 
    horizontal[j] = loss;
  }

  /* skip one line ("VERTICAL 360")*/
//  (void) fgets (buffer, 250, in);
    if( fgets (buffer, 250, in)); // (void) not working, using if instead (workaround)

  /* read vertical data - one angle degree per step */
  for (j = 0; j < 360; j++)
  {
//    (void) fgets (buffer, 250, in);
    if( fgets (buffer, 250, in)); // (void) not working, using if instead (workaround)
    sscanf (buffer, "%lf %lf", &angle, &loss);
    if (j != (int)angle) G_fatal_error(_("Bad antenna diagram format")); 
    vertical[j] = loss;
  }

  fclose(in);
  

  /* (3) process the input pathloss data */

  double height_diff_Tx_Rx; /* difference of total heights (Tx - Rx) */
  double dist_Tx_Rx; /* distance between receiver and transmitter */
  double rec_east, rec_north;  /* receiver coordinates */
  double d_east, d_north; /* differences between receiver and transmitter coordinates (e.g. rec_east - tr_east)*/
  double hor_coor_angle, vert_coor_angle; /* angle between x-axis of coordinate system and the line between RX and TX - (hor)izontal and (vert)ical*/
  double hor_diag_angle, vert_diag_angle; /* angle in the antenna's diagram (antenna's direction is taken into account) - (hor)izontal and (vert)ical*/
  double horizontal_loss, vertical_loss; /* pathloss due to antenna's diagram */
  double temp_angle; /* temporary angle */


  /* for each row */
  for (row = 0; row < nrows; row++) 
  {
	      
    if (verbose)
      G_percent(row, nrows, 2);

    FCELL f_in, f_in_dem;  
    /* read input map */
    Rast_get_row(infd, inrast, row, FCELL_TYPE);
	
    Rast_get_row(infd2, inrast2, row, FCELL_TYPE);
	  
    /* process the data */
    for (col = 0; col < ncols; col++) 
    { 
      f_in = ((FCELL *) inrast)[col];
      f_in_dem = ((FCELL *) inrast2)[col];
	  

      /* Skip calculation if height undefined */
      if (Rast_is_f_null_value(&f_in_dem))
      {
        Rast_set_f_null_value(&((FCELL *) outrast)[col], 1);   
        continue;
      }


      /* calculate receiver coordinates */
//      rec_east = window.west + (col * window.ew_res);
//      rec_north = window.north - (row * window.ns_res);
      rec_east = window.west + window.ew_res/2.0 + (col * window.ew_res);
      rec_north = window.north - window.ns_res/2.0 - (row * window.ns_res);	

      /* calculate differences between receiver and transmitter coordinates */
      d_north = rec_north - north;
      d_east = rec_east - east;
	    
      /* calculate distance between Tx and Rx */
      dist_Tx_Rx = sqrt( pow(d_east,2) + pow(d_north,2) );
      dist_Tx_Rx=dist_Tx_Rx/1000;

	   
      /* If distance between Rx and Tx exceeds given radius, continue with other cells */
      if (dist_Tx_Rx > radius)
      {
        Rast_set_f_null_value(&((FCELL *) outrast)[col], 1);   
        continue;
      }


      /* determine horizontal angle and loss */			
      if ( d_north == 0.0) temp_angle = PI / 2.0;
      else temp_angle = atan (d_east / d_north);

      if (temp_angle < 0) temp_angle = - temp_angle;
	    	      
      if      (d_north >= 0 && d_east >= 0)    hor_coor_angle = temp_angle;
      else if (d_north >= 0 && d_east < 0 )    hor_coor_angle = 2*PI - temp_angle;
      else if (d_north < 0  && d_east < 0 )    hor_coor_angle = PI + temp_angle;
      else /* (d_north < 0  && d_east >= 0) */ hor_coor_angle = PI - temp_angle;

      hor_coor_angle = hor_coor_angle * 180 / PI;  /* convert from radians to degrees */

      hor_diag_angle = hor_coor_angle - (double)beamDirection;

	
      if (hor_diag_angle < 0) hor_diag_angle = 360 + hor_diag_angle;

      /* to prevent reading unallocated data (diagram comprises values 0 - 359)*/
      temp_angle = ceil(hor_diag_angle);
      if (temp_angle == 360) temp_angle = 0;

      horizontal_loss = horizontal[(int)floor(hor_diag_angle)] +
                        ((horizontal[(int)temp_angle] - horizontal[(int)floor(hor_diag_angle)])*(hor_diag_angle - floor(hor_diag_angle))); /* interpolation */
	    
      /* determine vertical angle and loss */
      height_diff_Tx_Rx = totalHeight - (double)f_in_dem - rec_height;

      if ( dist_Tx_Rx == 0.0) vert_coor_angle = PI / 2.0;
      else vert_coor_angle = atan (height_diff_Tx_Rx / (dist_Tx_Rx * 1000));
      vert_coor_angle = vert_coor_angle * 180 / PI;	

      if (vert_coor_angle < 0) vert_coor_angle = 360 + vert_coor_angle;

// |-->
// 3.1.2012 - Vilhar
// Calculate the impact of mechanical tilt with respect to horizontal angle. At 0 degrees, the value is the same as the input. At 180 degrees, the input value is negative. Inbetween, we interpolate. This correction does not contribute essentially, but nevertheless should improve the final result slightly. 
 
      double mechanicalAntennaTilt_Corrected;

      if (hor_diag_angle >= 0 && hor_diag_angle <= 180)
          mechanicalAntennaTilt_Corrected = (double)mechanicalAntennaTilt * (1 - (hor_diag_angle / 90));
      else if (hor_diag_angle > 180 && hor_diag_angle <= 360)
          mechanicalAntennaTilt_Corrected = (double)mechanicalAntennaTilt * ((hor_diag_angle / 90) - 3);
      else G_fatal_error(_("Horizontal angle is not between 0 and 360 degrees.")); 

// -->|		

      vert_diag_angle = vert_coor_angle - (double)mechanicalAntennaTilt_Corrected;

      if (vert_diag_angle < 0) vert_diag_angle = 360 + vert_diag_angle;
      if (vert_diag_angle > 360.0) vert_diag_angle = vert_diag_angle - 360.0; // I.O. jul2017


      /* to prevent reading unallocated data (diagram comprises values 0 - 359)*/
      temp_angle = ceil( vert_diag_angle);
      if (temp_angle == 360) temp_angle = 0;

      vertical_loss = vertical[ (int)floor( vert_diag_angle)] +
                      ( ( vertical[ (int)temp_angle] - vertical[ (int)floor( vert_diag_angle)]) * (vert_diag_angle - floor( vert_diag_angle))); /* interpolation */ 
	      
      /* finally take into account pathloss for determined diagram angles and antenna gain*/
      ((FCELL *) outrast)[col] = (FCELL)((double)f_in + horizontal_loss + vertical_loss - gain); 		
    } // for col ...

 
    /* write raster row to output raster map */
    Rast_put_row(outfd, outrast, FCELL_TYPE);
  } // for row ...


  /* memory cleanup */
  G_free(inrast);
  G_free(outrast);

  /* closing raster maps */
  Rast_close(infd);
  Rast_close(outfd);

  /* add command line incantation to history file */
  Rast_short_history(result, "raster", &history);
  Rast_command_history(&history);
  Rast_write_history(result, &history);


  exit(EXIT_SUCCESS);
}

