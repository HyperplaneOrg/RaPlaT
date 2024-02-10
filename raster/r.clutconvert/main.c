/****************************************************************************
 *
 * MODULE:       r.clutconvert
 * AUTHOR(S):    Andrej Hrovat, Jozef Stefan Institute
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Convert land usage code from clutter file to 
 *               path loss factors
 *             
 *
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:

04-nov-2012
 - new format of the Path_loss_values file (lines contain <terrain_type:terrain_loss> value pairs)

04-feb-2013 (I. Ozimek):
 - modified processing of Path_loss_values file lines (to avoid segmentation fault errors)
 - any line starting with '#' is ignored (used for comments);
   the first line is now processed normally (previously ignored and used as a comment line)

17-aug-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
 - minor corrections (to remove compiler warnings, redundant code, etc.)

18-aug-2017 (I.O.)
 - code clenaup and minor corrections
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>

/*
 * main function
 */
int main( int argc, char *argv[])
{

  struct Cell_head cellhd;  /* it stores region information, and header information of rasters */
  char *name, *nameMapfile; /* input raster name, landuse mapping file name */
  char *result;             /* output raster name */
  const char *mapset;       /* mapset name */
  void *inrast;             /* input buffer */
  unsigned char *outrast;   /* output buffer */
  int nrows, ncols;
  int row, col;
  int infd, outfd;          /* file descriptor */
  struct History history;   /* holds meta-data (title, comments,..) */

  struct GModule *module;   /* GRASS module for parsing arguments */

  struct Option *input, *output, *input2; /* options */
  FILE *inf;                 /*file for Path loss factors*/


  /* initialize GIS environment */
  G_gisinit( argv[ 0]); /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword( _("raster"));
  G_add_keyword( _("clutter"));
  module->description = _("Clutter convert module (v18aug2017)");

  /* Define the different options as defined in gis.h */
  input = G_define_standard_option( G_OPT_R_INPUT);
  input->description = _("Input raster map - land usage categories");

  input2 = G_define_standard_option( G_OPT_R_INPUT);
  input2->key = "landuse_to_pathloss";
  input2->type = TYPE_STRING;
/*
  char mapfile_path[1000];
  strcpy( mapfile_path, getenv("GISBASE"));
  strcat( mapfile_path, "/etc/radio_coverage/lossfactors.txt");
  input2->answer = mapfile_path;
*/
  input2->required = YES;
  input2->gisprompt = "old_file,file,input";
  input2->description = _("Input text file - mapping 'land usage' -> 'RaPlaT path loss' (e.g. clutter map for hataDEM model)");
 
  output = G_define_standard_option( G_OPT_R_OUTPUT);
  output->description = _("Output raster map - RaPlaT path loss (e.g. clutter map for hataDEM model)");


  /* options and flags parser */
  if ( G_parser( argc, argv))
    exit(EXIT_FAILURE);

  /* stores options and flags to variables */
  name = input->answer;
  nameMapfile = input2->answer;
  result = output->answer;

  /* open landusage mapping file */
  if ( ( inf = fopen( nameMapfile, "r")) == NULL )
    G_fatal_error( _("Unable to open file <%s>"), nameMapfile);  


  /* returns NULL if the map was not found in any mapset, mapset name otherwise */
  mapset = G_find_raster( name, "");
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

  G_set_window( &cellhd);
  G_get_set_window( &cellhd);

  /* Allocate input buffer */
  inrast = Rast_allocate_buf( FCELL_TYPE);

  /* Allocate output buffer, use input map data_type */
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf( FCELL_TYPE);

  /* controlling, if we can write the raster */
  if ( ( outfd = Rast_open_new( result, FCELL_TYPE)) < 0)
    G_fatal_error( _("Unable to create raster map <%s>"), result);


  char buffer[ 1024];
  char buffer2[ 1024];
  char *token; // substring with numeric value land usage number or path loss value
  char *eptr;  // pointer to the (sub)string after the numerical value (read by strod), should be end of token ('\0')
  
  double terr_category[ 100];
  double terr_pathloss[ 100];

  int counter = 0;  

  while ( fgets( buffer, 1024, inf) != NULL)
  {
    if ( *buffer == '\n' || *buffer == '#') continue;  // skip empty and comment lines

    if ( counter > 99)
      G_fatal_error( _("Maximum number of categories exceeded (100) "));

    /* land usage number */
    strncpy( buffer2, buffer, sizeof( buffer));
    token = strtok( buffer2, ":");
    if ( token == NULL)
      G_fatal_error( _("Land usage number or ':' missing in the following landuse mapping file line:\n %s"), buffer);
    terr_category[ counter] = strtod( token, &eptr);
    if ( *eptr != '\0')
      G_fatal_error( _("Wrong land usage numeric value in the following landuse mapping file line:\n %s"), buffer);

    /* path loss number */
    token = strtok( NULL, "\n\0");
    if ( token == NULL)
      G_fatal_error( _("Path loss value missing in the following landuse mapping file line:\n %s"), buffer);
    terr_pathloss[ counter] = strtod( token, &eptr);
    if ( *eptr != '\0')
      G_fatal_error( _("Wrong pathloss numeric value in the following landuse mapping file line:\n %s"), buffer);

    counter++;
  }

    
  /* for each row */
  for ( row = 0; row < nrows; row++) 
  {  
    FCELL f_in;

    /* read input map */
    Rast_get_row( infd, inrast, row, FCELL_TYPE);

    /* process the data */
    for ( col = 0; col < ncols; col++) 
    { 
      f_in = ( (FCELL *)inrast)[col];
      FCELL f_out = 999;

      // NULL value in the input clutter
      if ( Rast_is_null_value( &f_in, FCELL_TYPE))
      {
        ( (FCELL *)outrast)[ col] = f_in;
        continue;
      }

      int count_array = 0;
      for ( count_array = 0; count_array < counter; count_array++)
      {
        if ( terr_category[ count_array] == (double)f_in)
        {
          f_out = terr_pathloss[ count_array];
          break;
        }
      }

      if ( f_out == 999)
        G_fatal_error( _("Land usage value %f not found in landuse mapping file"), (double)f_in);
    
      ( (FCELL *)outrast)[ col] = f_out;
    }
      
    /* write raster row to output raster map */
    Rast_put_row( outfd, outrast, FCELL_TYPE);
  }
         

  /* memory cleanup */
  G_free( outrast);
  G_free( inrast);

  /* closing raster maps */
  Rast_close( infd);
  Rast_close( outfd);

  /* add command line incantation to history file */
  Rast_short_history( result, "raster", &history);
  Rast_command_history( &history);
  Rast_write_history( result, &history);

  exit( EXIT_SUCCESS);
}

