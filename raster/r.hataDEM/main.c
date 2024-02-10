/****************************************************************************
 *
 * MODULE:       r.hataDEM
 * AUTHOR(S):    Andrej Vilhar, Tomaz Javornik, Andrej Hrovat Jozef Stefan Institute                
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Calculates radio coverage from a single base station 
 *               according to model 9999
 *
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:
4-feb-2013 (I.O.)
 - minor modification

23-sep-2013 (I.O.)
 - minor modification

2-oct-2013 (I.O.)
 - corrected one pixel shift (right- and downwards) in the output map
 - redundant code disabled (commented out)

15-nov-2013 (I.O.)
 - additional parameter rx_ant_height (default = 1.5)
 - additional parameter clut_mode (rx/tx/none, default = tx)

9-june-2014 (I.O.)
 - redundant code removed (previously commented out)
 - C source style cleaned-up
 - corrected an error in elevation angle calculation (for NLOS conditions) - but not used any more - 
    - corrections related to the elevation angle disabled (not needed / wrong - http://www.mike-willis.com/Tutorial/PF7.htm)
 - "inverse" mode added (RX and TX roles exchanged, needed for transmitter localization)

12-june-2014 (I.O.)
  - changed processing of Zeff (effective antenna height) for small/negative values:
    Zeff is know downwards limited to AntHeightBS (transmitter antenna height)

4-aug-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
     option name changed: inputDEM -> input_dem; A0..A3 -> a0..a3 
 - minor corrections (to remove compiler warnings, redundant code, etc.)

7-dec-2018
 - minor changes to avoid C compiler warnings
*/


/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>
*/

#include "local_proto.h"


struct StructhataDEM
{
  double BSxIndex;     /* normalized position of BS -> UTMx/resolution */
  double BSyIndex;     /* normalized position of BS -> UTMx/resolution */
  double BSAntHeight;  /* Antenna height of BS [m] */
  double MSAntHeight;  /* Antenna height of MS [m] */
  int xN;              /* dimension of teh input(Raster) and output (PathLoss) */
  int yN;              /* dimension of teh input(Raster) and output (PathLoss) */
  double scale;        /* Resolution of DEM file */
  double freq;         /* Carrier frequency in MHz */
  double A0;           /* Model 9999 parameters */
  double A1;           /* Model 9999 parameters */
  double A2;           /* Model 9999 parameters */
  double A3;           /* Model 9999 parameters */
  double ResDist;      /* Resolution model 9999 profile calc */
  double radi;         /* Radius of calculation [km] */
};

int hataDEMPathLossSub( double**, double**, double**, struct StructhataDEM, int, int);


/*
 * main function
 */
int main(int argc, char *argv[])
{
  double east;
  double north;
  double ant_height, frequency, radius;
//  double rec_height = 1.5; /* height of receiver from the ground */
  double rec_height;

  struct Cell_head window;  /*  database window         */
  struct Cell_head cellhd;  /* it stores region information, and header information of rasters */
  char *name,  *name2;      /* input raster name */
  char *result;             /* output raster name */
  const char *mapset, *mapset2;   /* mapset name */
  void *inrast, *inrast2 = NULL;   /* input buffer */
  unsigned char *outrast;   /* output buffer */
  int nrows, ncols;
  int row, col;
  int tr_row, tr_col;
  int infd, outfd, infd2 = 0;   /* file descriptor */
  int verbose;
  int inverse_mode_f;

  struct History history;   /* holds meta-data (title, comments,..) */

  struct GModule *module;   /* GRASS module for parsing arguments */

  struct Option *input, *opt1, *opt2, *opt4, *opt3, *opt5, *opt6, *opt7, *opt8, *opt9, *opt10, *output, *input2; /* options */
  struct Flag *flag1, *flag2;       /* flags */

  double A0_main, A1_main, A2_main, A3_main;  /* model 9999 parameters*/


  /* initialize GIS environment */
  G_gisinit(argv[0]);  /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("hataDEM"));
  module->description = _("RaPlaT - HataDEM module (v07dec2018)");

  /* define different options as defined in gis.h */
  input = G_define_standard_option( G_OPT_R_INPUT); 
  input->key = "input_dem";       
    
  opt9 = G_define_option();
  opt9->key = "clut_mode";
  opt9->type = TYPE_STRING;
  opt9->required = NO;
  opt9->options = "rx,tx,none";
  opt9->answer = "rx";
  opt9->description = _("Clutter usage");

  input2 = G_define_standard_option( G_OPT_R_INPUT);
  input2->key = "clutter";
  input2->required = NO;
  input2->answer = "";
  input2->description = _("Name of clutter raster map with path loss coefficients");

  output = G_define_standard_option( G_OPT_R_OUTPUT);

  /* define  different flags */
  flag1 = G_define_flag();
  flag1->key = 'q';
  flag1->description = _("Quiet");

  flag2 = G_define_flag();
  flag2->key = 'i';
  flag2->description = _("Inverse mode (RX and TX roles exchanged)");


  /* hataDEM parameters */
  opt5 = G_define_option();
  opt5->key = "a0";
  opt5->type = TYPE_DOUBLE;
  opt5->required = YES;
  opt5->description = _("Parameter a0");

  opt6 = G_define_option();
  opt6->key = "a1";
  opt6->type = TYPE_DOUBLE;
  opt6->required = YES;
  opt6->description = _("Parameter a1");

  opt7 = G_define_option();
  opt7->key = "a2";
  opt7->type = TYPE_DOUBLE;
  opt7->required = YES;
  opt7->description = _("Parameter a2");

  opt8 = G_define_option();
  opt8->key = "a3";
  opt8->type = TYPE_DOUBLE;
  opt8->required = YES;
  opt8->description = _("Parameter a3");
  /* end of hataDEM parameters */


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

  opt10 = G_define_option();
  opt10->key = "rx_ant_height";
  opt10->type = TYPE_DOUBLE;
  opt10->required = NO;
  opt10->answer = "1.5";
  opt10->description = _("Receiver antenna height [m]");

  opt3 = G_define_option();
  opt3->key = "frequency";
  opt3->type = TYPE_DOUBLE;
  opt3->required = YES;
  opt3->description = _("Frequency [MHz]");

  /* options and flags parser */
  if (G_parser(argc, argv))
    exit(EXIT_FAILURE);

  /* stores options and flags to variables */
  name = input->answer;
  name2 = input2->answer;
  result = output->answer;
  verbose = !flag1->answer;
  inverse_mode_f = flag2->answer;

  G_scan_easting( opt1->answers[0], &east, G_projection());
  G_scan_northing( opt1->answers[1], &north, G_projection());
  sscanf( opt2->answer, "%lf", &ant_height);
  sscanf( opt4->answer, "%lf", &radius);
  sscanf( opt3->answer, "%lf", &frequency); 

  sscanf( opt5->answer, "%lf", &A0_main);
  sscanf( opt6->answer, "%lf", &A1_main);
  sscanf( opt7->answer, "%lf", &A2_main);
  sscanf( opt8->answer, "%lf", &A3_main);

  sscanf( opt10->answer, "%lf", &rec_height);


  // clutter mode (default tx)
  int clutmode = 0;
  if ( strcmp( opt9->answer, "rx") == 0) clutmode = 1;
  else if ( strcmp( opt9->answer, "tx") == 0) clutmode = 2;

  if ( clutmode)
    if ( strlen( name2) == 0)
      G_fatal_error( _("No clutter map specified"));


  /* returns NULL if the map was not found in any mapset, mapset name otherwise */
  mapset = G_find_raster( name, "");
  if ( mapset == NULL)
    G_fatal_error( _("Raster map <%s> not found"), name);
   
  if ( clutmode)
  {
    mapset2 = G_find_raster( name2, "");
    if ( mapset2 == NULL)
      G_fatal_error( _("Raster map <%s> not found"), name2);
  }

  if ( G_legal_filename(result) < 0)
    G_fatal_error( _("<%s> is an illegal file name"), result);


  /* Rast_open_old - returns file destriptor (>0) */
  if ( ( infd = Rast_open_old( name, mapset)) < 0)
    G_fatal_error( _("Unable to open DEM raster map <%s>"), name);

  if ( clutmode)
    if ( ( infd2 = Rast_open_old( name2, mapset2)) < 0)
      G_fatal_error( _("Unable to open clutter raster map <%s>"), name2);

  /* open input raster */   
  Rast_get_cellhd( name, mapset, &cellhd);


  G_debug( 3, "number of rows %d", cellhd.rows);

  G_get_window( &window);

  /* allocate input GRASS buffers */
  inrast = Rast_allocate_buf( FCELL_TYPE);
  if ( clutmode)
    inrast2 = Rast_allocate_buf( FCELL_TYPE);

  /* allocate output GRASS buffer, use input map data_type */
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf( FCELL_TYPE);

  /* check if we can write the raster */
  if ( ( outfd = Rast_open_new( result, FCELL_TYPE)) < 0)
    G_fatal_error( _("Unable to create raster map <%s>"), result);

  /* check if specified transmitter location lies inside GRASS region window */
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
  Rast_get_row( infd, inrast, tr_row, FCELL_TYPE);
  trans_elev = ( (FCELL *) inrast)[ tr_col];
    
  // check if transmitter is on DEM
  if ( isnan( (double) trans_elev))
  {
    if ( !inverse_mode_f)
      G_fatal_error( _("Transmitter outside raster DEM map."));
    else
      G_fatal_error( _("Receiver outside raster DEM map."));
  }


  /* define structure variables */
  double BSAntHeight = ant_height;     
  double MSAntHeight = rec_height;     
  int xN = window.rows;
  int yN = window.cols;
  double scale = window.ew_res;
  double freq = frequency;
  double A0 = A0_main;
  double A1 = A1_main;
  double A2 = A2_main;
  double A3 = A3_main;
  double ResDist = 1;

  double BSyIndex = ( east - window.west) / scale - 0.5;
  double BSxIndex = ( window.north - north) / scale - 0.5;

  double radi = radius;

  struct StructhataDEM InihataDEM = { BSxIndex, BSyIndex, BSAntHeight, MSAntHeight, xN, yN, scale, freq, A0, A1, A2, A3, ResDist, radi};


  /* do hataDEM */


  // allocate and initialize buffers for dem, clutter and output (path loss) raster maps
  // dem (input)
  double *tmp_rast; // pointer to raster the first pixel
  double **m_rast;  // array of pointers to raster lines (the first pixel of eachline)
  int i;
  tmp_rast = (double *) G_calloc( nrows * ncols, sizeof( double));
  m_rast = (double **) G_calloc( nrows, sizeof( double *));
  memset( tmp_rast, 0, nrows * ncols * sizeof( double));  // set all pixels values to zero
  for ( i = 0; i < nrows ; i++)
    m_rast [ i] = tmp_rast + i * ncols; // initialize pointers to raster lines

  // clutter (input)
  double *tmp_clut;
  double **m_clut;
  tmp_clut = (double *) G_calloc( nrows * ncols, sizeof( double));
  m_clut = (double **) G_calloc( nrows, sizeof( double *));
  memset( tmp_clut, 0, nrows * ncols * sizeof( double));
  for ( i = 0; i < nrows ; i++)
    m_clut[ i] = tmp_clut + i * ncols;

  // path loss (output)
  double *tmp_loss;
  double **m_loss;
  tmp_loss = (double *) G_calloc( nrows * ncols, sizeof( double));
  m_loss = (double **) G_calloc( nrows, sizeof( double *));
  memset( tmp_loss, 0, nrows * ncols * sizeof( double));
  for ( i = 0; i < nrows ; i++)
    m_loss[ i] = tmp_loss + i * ncols;
    

  // read GRASS maps (dem and clutter) to the application's raster buffers
  for ( row = 0; row < nrows; row++) 
  {
    //G_message( _("row: %d, col: %d"), row, col);  
    if ( verbose)
      G_percent( row, nrows, 2);

    FCELL f_in, f_in2;
    /* read one raster line of dem map */
    Rast_get_row( infd, inrast, row, FCELL_TYPE);
    /* read one raster line of clutter map */
    if ( clutmode)
      Rast_get_row( infd2, inrast2, row, FCELL_TYPE);

    /* copy the raster lines to the application's raster buffers */
    for ( col = 0; col < ncols; col++) 
    { 
      f_in = ( (FCELL *) inrast)[ col];
      m_rast[ row][ col] = (double) f_in;

      if ( clutmode)
      {
        f_in2 = ( (FCELL *) inrast2)[ col];
        m_clut[ row][ col] = (double) f_in2;
      }
    }
  }


  hataDEMPathLossSub (m_rast, m_clut, m_loss, InihataDEM, clutmode, inverse_mode_f); // hataDEM model computation


  // write the computed output raster (path loss) to the GRASS map (convert zero values to GRASS null values)
  double path_loss_num;
  FCELL  null_f_out;

  Rast_set_f_null_value( &null_f_out, 1);   

  for ( row = 0; row < nrows; row++)
  {
    G_percent( row, nrows, 2);
    for ( col = 0; col < ncols; col++) 
    {
      path_loss_num = m_loss[ row][ col];
      if ( path_loss_num == 0)
        ( (FCELL *) outrast)[ col] = null_f_out;
      else
        ( (FCELL *) outrast)[ col]  = (FCELL) path_loss_num;
    }
    /* write raster row to output raster map */
    Rast_put_row( outfd, outrast, FCELL_TYPE);
  }

  /* memory cleanup */
  G_free( inrast);
  G_free( outrast);
  if ( clutmode) G_free( inrast2);


  /* closing raster maps */
  Rast_close( infd);
  Rast_close( outfd);
  if ( clutmode) Rast_close( infd2);

  /* add command line incantation to history file */
  Rast_short_history( result, "raster", &history);
  Rast_command_history( &history);
  Rast_write_history( result, &history);

  exit( EXIT_SUCCESS);
}



/***********************************************************************************************
 *
 *  Function hataDEMPathLossSub calculates PathLoss in dB using hataDEM 9999 path loss formula
 *
 *    **PathLoss:     array of path loss in dB
 *    **Raster:       input DEM file
 *    **Clutter:      input clutter file
 *
 *    T.Javornik, Jan. 2010
 *
 ***********************************************************************************************/

int hataDEMPathLossSub( double** Raster, double** Clutter, double** PathLoss, struct StructhataDEM InihataDEM, int clutmode, int inverse_mode_f)
{
  // hataDEM model constants and variables
  double BSxIndex = 0;                          // normalized position of BS -> UTMx/resolution 
  double BSyIndex = 0;                          // normalized position of BS -> UTMy/resolution
  double AntHeightBS = InihataDEM.BSAntHeight;  // antenna height of BS [m]
  double AntHeightMS = InihataDEM.MSAntHeight;  // antenna height of MS [m]
  int xN = InihataDEM.xN;                       // dimension of the input(Raster) and output (PathLoss)
  int yN = InihataDEM.yN;                       // dimension of the input(Raster) and output (PathLoss)
  double scale = InihataDEM.scale;              // Model 9999 resolution
  double A0 = InihataDEM.A0;                    // Model 9999 parameters
  double A1 = InihataDEM.A1;                    // Model 9999 parameters
  double A2  = InihataDEM.A2;                   // Model 9999 parameters
  double A3  = InihataDEM.A3;                   // Model 9999 parameters
  double freq  = InihataDEM.freq;               // carrier frequency
  double ResDist = InihataDEM.ResDist;          // distance BS - MS sampling rate [normalized with scale
  double Lambda = 300.0 / freq;                 // wavelenght
  double radi = InihataDEM.radi;                // radius of calculation

  double ZoBS;                  // BS height about the sea level
  double ZObs2LOS = 0;
  double DistObs2BS = 0;
  double ZoTransBS,ZoTransMS;   // BS and MS antenna heights about the sea level
  double log10Zeff;
  double log10DistBS2MSKm;
  double tiltBS2MS;             // (ZoBS-ZoMS)/distBS2MSNorm
  double PathLossFreq = 0;      // path loss due to carrier frequency
  double PathLossTmp = 0;       // tmp path loss
  int ix; int iy;
  double DiffX, DiffY, Zeff;    // Difference in X and Y direction
  double PathLossAntHeightBS;
  double DistBS2MSNorm;         // normalized distance between MS and BS in xy plan sqrt(x2+y2)
  double DistBS2MSKm;           // distance between MS and BS in Km sqrt(x2+y2+z2) * scale / 1000
//  double ElevAngCos, Hdot, Ddot, Ddotdot, PathLossDiff;
  double Hdot, Ddot, Ddotdot, PathLossDiff;

  double MSxIndex = 0, MSyIndex = 0;  // normalized position of MS -> UTMx/resolution 
  int intMSxIndex, intMSyIndex;       // integer (rounded) values (for referencing the raster buffer)
  int intBSxIndex, intBSyIndex;       // integer (rounded) values (for referencing the raster buffer)


  if ( !inverse_mode_f)
  {
    BSxIndex = InihataDEM.BSxIndex;
    BSyIndex = InihataDEM.BSyIndex;
  }
  else
  {
    MSxIndex = InihataDEM.BSxIndex;
    MSyIndex = InihataDEM.BSyIndex;
  }


  PathLossFreq = 44.49 * log10( freq) - 4.78 * pow( log10( freq), 2);  // loss due to carrier frequency
  PathLossAntHeightBS = 3.2 * pow( log10( 11.75 * AntHeightMS), 2);    // negative loss due to MS antenna height

  for ( ix = 0; ix < xN; ix++)
  {
    G_percent( ix, xN, 2);
    for ( iy = 0; iy < yN; iy++)
    {

      if ( !inverse_mode_f)
      {
        MSxIndex = (double) ix;
        MSyIndex = (double) iy;
      }
      else
      {
        BSxIndex = (double) ix;
        BSyIndex = (double) iy;
      }

      // path Loss due to Hata model
      DiffX = BSxIndex - MSxIndex;
      DiffY = BSyIndex - MSyIndex;
      DistBS2MSNorm = sqrt( DiffX * DiffX + DiffY * DiffY);
      DistBS2MSKm = DistBS2MSNorm * scale / 1000;
      if ( ( DistBS2MSKm) > radi)
        continue;
      if ( DistBS2MSKm < 0.01)
        DistBS2MSKm = 0.01;

      intBSxIndex = (int)( BSxIndex + 0.5);
      intBSyIndex = (int)( BSyIndex + 0.5);
      intMSxIndex = (int)( MSxIndex + 0.5);
      intMSyIndex = (int)( MSyIndex + 0.5);

      ZoBS = Raster[ intBSxIndex][ intBSyIndex]; // BS height above the sea level calculated from raster DEM file
      ZoTransBS = ZoBS + AntHeightBS;            // BS transmitter antenna height above the sea level
      ZoTransMS = Raster[ intMSxIndex][ intMSyIndex] + AntHeightMS; // MS receiver antenna height above the sea level
      Zeff = ZoTransBS - ZoTransMS;
//      if ( ZoBS <= Raster[ intMSxIndex][ intMSyIndex]) 
      if ( Zeff < AntHeightBS)
        Zeff = AntHeightBS;
      log10Zeff = log10( Zeff);

      log10DistBS2MSKm = log10( DistBS2MSKm);

      PathLossTmp = A0 + A1 * log10DistBS2MSKm; 
      PathLossTmp += A2 * log10Zeff + A3 * log10DistBS2MSKm * log10Zeff;
      PathLossTmp += PathLossFreq - PathLossAntHeightBS;

      // calc position of the height and position of the highest obstacle
      tiltBS2MS = ZoTransBS - ZoTransMS;
      if ( DistBS2MSNorm > 0)
        tiltBS2MS = -tiltBS2MS / DistBS2MSNorm;
      else
        tiltBS2MS = 0;

      DoProfile( &ZObs2LOS, &DistObs2BS, ResDist, Raster, BSxIndex, BSyIndex, ZoTransBS, MSxIndex, MSyIndex, tiltBS2MS);

      // calc path loss due to NLOS conditions

      // ElevAngCos - related corrections not needed (wrong) (http://www.mike-willis.com/Tutorial/PF7.htm)
      Hdot = ZObs2LOS;
      Ddot = DistObs2BS; 
      Ddotdot = DistBS2MSNorm - Ddot;
      /*
      ElevAngCos = cos( atan( tiltBS2MS / scale));
      Hdot = ZObs2LOS * ElevAngCos;
      Ddot = DistObs2BS; 
      Ddotdot = DistBS2MSNorm - Ddot;
      if ( ElevAngCos != 0)
      {
        Ddot /= ElevAngCos;
        Ddotdot /= ElevAngCos;
      }
      */

      PathLossDiff = 0;
      if ( Ddot > 0 && Ddotdot > 0)
      {
        PathLossDiff = Hdot * sqrt( 2 * ( Ddot + Ddotdot) / ( Lambda * Ddot * Ddotdot * scale));
        if ( PathLossDiff < -0.75 )
          PathLossDiff = 0; 
        else
        {
          PathLossDiff = PathLossDiff - 0.1;  // (http://www.mike-willis.com/Tutorial/PF7.htm)
          PathLossDiff = 6.9 + 20 * log10( sqrt( PathLossDiff * PathLossDiff + 1) + PathLossDiff);
        }
      }
      PathLossTmp += PathLossDiff;

      // add clutter and write data to pathloss
      if      ( clutmode == 1)  PathLossTmp += Clutter[ intMSxIndex][ intMSyIndex];
      else if ( clutmode == 2)  PathLossTmp += Clutter[ intBSxIndex][ intBSyIndex];

      PathLoss[ ix][ iy] = PathLossTmp;

    } // end irow
  } // end icol

  return 0;
}

