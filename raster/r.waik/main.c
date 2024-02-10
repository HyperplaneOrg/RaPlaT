
/****************************************************************************
 *
 * MODULE:       r.waik
 * AUTHOR(S):    Tomaz Javornik, Andrej Hrovat Jozef Stefan Institute                
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Calculates radio coverage from a single base station 
 *               according to Walfish-Ikegami model
 *             
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:
8-nov-2011
 - input parameters added (comparable to TEMS)

23-sep-2013 (I.O.)
 - minor modification

2-oct-2013 (I.O.)
 - corrected one pixel shift (right- and downwards) in the output map
 - redundant code disabled (commented out)

16-aug-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
     option name changed: inputDEM -> input_dem, PHI_Street -> phi_street
 - minor corrections (to remove compiler warnings, redundant code, etc.)

7-dec-2018
 - identation correected to avoid C compiler warnings
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


#define PI 3.14159265

struct StructWaIk
{
  double BSxIndex;    /* normalized position of BS -> UTMx/resolution */  
  double BSyIndex;    /* normalized position of BS -> UTMx/resolution */
  double BSAntHeight; /* Antenna height of BS [m] */
  double MSAntHeight; /* Antenna height of MS [m] */
  int xN;             /* dimension of teh input(Raster) and output (PathLoss) */
  int yN;             /* dimension of teh input(Raster) and output (PathLoss) */
  double scale;       /* Resolution of DEM file */
  double freq;        /* Carrier frequency in MHz */  
    
  double W0;          /* Free space loss correction*/
  double W1;          /* Reduced base antenna height correction*/
  double W2;          /* Range correction*/
  double W3;          /* Street width correction*/
  double W4;          /* Frequency correction*/
  double W5;          /* Building height correction*/
  double W6;          /* Street width [m] */
  double W7;          /* Distance between buildings [m] */
  double W8;          /* Building height [m] */
  double PHIStreet;   /* Street orientation*/

  double ResDist;     /* Resolution Walfish-Ikegami model profile calc */
  double radi;        /* Radius of calculation [km] */    
};


int WaIkPathLossSub(double**, double**, struct StructWaIk, char*);


/*
 * main function
 */
int main(int argc, char *argv[])
{
  double east;
  double north;
  double ant_height, frequency, radius; //, dem_height, clut_value;
  double rec_height = 1.5; /* height of receiver from the ground */
    
  struct Cell_head window; /* database window */
  struct Cell_head cellhd; /* it stores region information, and header information of rasters */
  char *name;              /* input raster name */
  char *result;            /* output raster name */
  const char *mapset;      /* mapset name */
  void *inrast;            /* input buffer */
  unsigned char *outrast;  /* output buffer */
  int nrows, ncols;
  int row, col;
  int tr_row, tr_col;
  int infd, outfd;         /* file descriptor */
  int verbose;
  struct History history;  /* holds meta-data (title, comments,..) */

  struct GModule *module;  /* GRASS module for parsing arguments */

  struct Option *input,*opt1, *opt2, *opt4, *opt3, *opt5, *opt6, *opt7, *opt8, *output, *opt9, *opt10, *opt11, *opt12, *opt13, *opt14, *opt15;    /* options */
  struct Flag *flag1;      /* flags */
  //FILE *in;              /*file for WaIk param*/
  double W0_main, W1_main, W2_main, W3_main, W4_main, W5_main, W6_main, W7_main, W8_main, PHI_Street_main; /* WaIk model parameters */

  /* initialize GIS environment */
  G_gisinit(argv[0]); /* reads grass env., stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("waik"));
  module->description = _("RaPlaT - Walfish-Ikegami module (v07dec2018)");

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
 
  opt3 = G_define_option();
  opt3->key = "frequency";
  opt3->type = TYPE_DOUBLE;
  opt3->required = YES;
  opt3->description = _("Frequency [MHz]");

  opt4 = G_define_option();
  opt4->key = "radius";
  opt4->type = TYPE_DOUBLE;
  opt4->required = NO;
  opt4->answer = "10";
  opt4->description = _("Computation radius [km]");

  /*WA-IK parameters ------------------*/
  opt10 = G_define_option();
  opt10->key = "free_space_loss_correction";
  opt10->type = TYPE_DOUBLE;
  opt10->required = NO;
  opt10->answer = "32.5";
  opt10->description = _("Free space loss correction W0 (20-60)");

  opt11 = G_define_option();
  opt11->key = "bs_correction";
  opt11->type = TYPE_DOUBLE;
  opt11->required = NO;
  opt11->answer = "54";
  opt11->description = _("Reduced base antenna height correction W1 (30-70)");

  opt12 = G_define_option();
  opt12->key = "range_correction";
  opt12->type = TYPE_DOUBLE;
  opt12->required = NO;
  opt12->answer = "10";
  opt12->description = _("Range correction W2 (5-35)");

  opt13 = G_define_option();
  opt13->key = "street_width_correction";
  opt13->type = TYPE_DOUBLE;
  opt13->required = NO;
  opt13->answer = "10";
  opt13->description = _("Street width correction W3 (3-15)");

  opt14 = G_define_option();
  opt14->key = "frequency_correction";
  opt14->type = TYPE_DOUBLE;
  opt14->required = NO;
  opt14->answer = "10";
  opt14->description = _("Frequency correction W4 (3-25)");

  opt15 = G_define_option();
  opt15->key = "building_height_correction";
  opt15->type = TYPE_DOUBLE;
  opt15->required = NO;
  opt15->answer = "20";
  opt15->description = _("Building Height Correction W5 (10-30)");

  opt5 = G_define_option();
  opt5->key = "street_width";
  opt5->type = TYPE_DOUBLE;
  opt5->required = NO;
  opt5->answer = "15";
  opt5->description = _("Street width W6 [m]");

  opt6 = G_define_option();
  opt6->key = "distance_between_buildings";
  opt6->type = TYPE_DOUBLE;
  opt6->required = NO;
  opt6->answer = "30";
  opt6->description = _("Distance between buildings W7 [m]");


  opt7 = G_define_option();
  opt7->key = "building_height";
  opt7->type = TYPE_DOUBLE;
  opt7->required = NO;
  opt7->answer = "12";
  opt7->description = _("Building height W8 [m]");

  opt9 = G_define_option();
  opt9->key = "phi_street";
  opt9->type = TYPE_DOUBLE;
  opt9->required = NO;
  opt9->answer = "90";
  opt9->description = _("Street orientation [deg]"); 


  opt8 = G_define_option();
  opt8->key = "area_type";
  opt8->type = TYPE_STRING;
  opt8->required = NO;
  opt8->description = _("Area type");
  opt8->options = "metropolitan,medium_cities";
  opt8->answer = "medium_cities";
  /*------------------------------------*/

  /* options and flags parser */
  if (G_parser(argc, argv))
  exit(EXIT_FAILURE);

  /* stores options and flags to variables */
  name = input->answer;
  result = output->answer;
  verbose = (!flag1->answer);
  G_scan_easting(opt1->answers[0], &east, G_projection());
  G_scan_northing(opt1->answers[1], &north, G_projection());
  sscanf(opt2->answer, "%lf", &ant_height);
  sscanf(opt4->answer, "%lf", &radius);
  sscanf(opt3->answer, "%lf", &frequency); 

  sscanf(opt10->answer, "%lf", &W0_main);
  sscanf(opt11->answer, "%lf", &W1_main);
  sscanf(opt12->answer, "%lf", &W2_main);
  sscanf(opt13->answer, "%lf", &W3_main);
  sscanf(opt14->answer, "%lf", &W4_main);
  sscanf(opt15->answer, "%lf", &W5_main);
  sscanf(opt5->answer, "%lf", &W6_main);
  sscanf(opt6->answer, "%lf", &W7_main);
  sscanf(opt7->answer, "%lf", &W8_main);
  sscanf(opt9->answer, "%lf", &PHI_Street_main);
        
  /* returns NULL if the map was not found in any mapset, mapset name otherwise */
  mapset = G_find_raster(name, "");
  if (mapset == NULL)
  G_fatal_error(_("Raster map <%s> not found"), name);
   
  if (G_legal_filename(result) < 0)
  G_fatal_error(_("<%s> is an illegal file name"), result);


  /* Rast_open_old - returns file destriptor (>0) */
  if ((infd = Rast_open_old(name, mapset)) < 0)
  G_fatal_error(_("Unable to open raster map <%s>"), name);

  /* open input raster */   
  Rast_get_cellhd(name, mapset, &cellhd);

  G_debug(3, "number of rows %d", cellhd.rows);

  G_get_window(&window);

  /* Allocate input buffer */
  inrast = Rast_allocate_buf(FCELL_TYPE);

  /* Allocate output buffer, use input map data_type */
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf(FCELL_TYPE);

  //G_message(_("nrows %d and ncols %d"),nrows,ncols);

  /* controlling, if we can write the raster */
  if ((outfd = Rast_open_new(result, FCELL_TYPE)) < 0)
  G_fatal_error(_("Unable to create raster map <%s>"), result);

  /* check if specified transmitter location inside window   */
  if (east < window.west || east > window.east
      || north > window.north || north < window.south)
    G_fatal_error(_("Specified base station  coordinates are outside current region bounds."));  

  /* map array coordinates for transmitter */
  tr_row = (window.north - north) / window.ns_res;
  tr_col = (east - window.west) / window.ew_res;

  /* total height of transmitter */
  FCELL trans_elev;
  Rast_get_row(infd, inrast, tr_row, FCELL_TYPE);
  trans_elev = ((FCELL *) inrast)[tr_col];

  // check if transmitter is on DEM 
  if ( isnan((double)trans_elev))                         
  {
    G_fatal_error(_("Transmitter outside raster DEM map."));
  }

  /*--- define structure variables----*/   
  double BSAntHeight = ant_height;     
  double MSAntHeight = rec_height;     
  int xN = window.rows;           
  int yN = window.cols;   
  double scale = window.ew_res;       
  double freq = frequency;

  double W0 = W0_main;    
  double W1 = W1_main;    
  double W2 = W2_main;    
  double W3 = W3_main;    
  double W4 = W4_main;    
  double W5 = W5_main;    
  double W6 = W6_main;        
  double W7 = W7_main;        
  double W8 = W8_main;        
  double PHIStreet  = PHI_Street_main;            

  double ResDist = 1;

//  (I.O. 2-oct-2013)
//  double BSyIndex = (east-window.west)/scale+0.5; 
//  double BSxIndex = (window.north-north)/scale+0.5;
  double BSyIndex = ( east - window.west) / scale - 0.5;
  double BSxIndex = ( window.north - north) / scale - 0.5;

  double radi = radius;

  struct StructWaIk IniWaIk = {BSxIndex,BSyIndex, BSAntHeight, MSAntHeight, xN, yN, scale, freq, W0, W1, W2, W3, W4, W5, W6, W7, W8, PHIStreet, ResDist, radi};
  /*---------------------------------*/    
    
/* do WA-IK */

  /* allocate the memory to contain the whole file */
  /*RASTER*/
  double **m_rast;
  int i;
  m_rast = (double **)G_calloc(nrows, sizeof(double *));
  //m_rast [0]= (double *)G_calloc(nrows * ncols, sizeof(double));
  double *tmp_rast = (double *)G_calloc(nrows * ncols, sizeof(double));
  memset (tmp_rast, 0, nrows * ncols * sizeof(double));
  for (i=0; i<nrows;i++)  
  {
    m_rast [i] = tmp_rast + i*ncols;
  }
  /*PATH LOSS*/
  double **m_loss;
  int k;
  m_loss = (double **)G_calloc(nrows, sizeof(double *));
  //m_rast [0]= (double *)G_calloc(nrows * ncols, sizeof(double));
  double *tmp_loss = (double *)G_calloc(nrows * ncols, sizeof(double));
  memset (tmp_loss, 0, nrows * ncols * sizeof(double));
  for (k=0; k<nrows;k++)  
  {
    m_loss[k] = tmp_loss + k*ncols;
  }

  /* Write file (raster) to array - for each row */

  for (row = 0; row < nrows; row++) 
  {   
    if (verbose)
      G_percent(row, nrows, 2);
    FCELL f_in;     
    /* read input map */
    Rast_get_row(infd, inrast, row, FCELL_TYPE);
    /* process the data */
    for (col = 0; col < ncols; col++) 
    { 
      f_in = ((FCELL *) inrast)[col];
      m_rast[row][col] = (double)f_in;
    }
  }
    
  WaIkPathLossSub (m_rast, m_loss, IniWaIk, opt8->answer); 
    
  double path_loss_num;
  FCELL  null_f_out;

  for (row = 0; row < nrows; row++)
  {
    G_percent(row, nrows, 2);
    for (col = 0; col < ncols; col++) 
    {
      path_loss_num = m_loss[row][col];
      if (path_loss_num == 0)
      {
        Rast_set_f_null_value(&null_f_out, 1);   
        ((FCELL *) outrast)[col] =null_f_out;
      }   
      else
      {
        ((FCELL *) outrast)[col]  = (FCELL)path_loss_num;
      } 
    }
    /* write raster row to output raster map */
    Rast_put_row(outfd, outrast, FCELL_TYPE);
  }


   // G_message(_("END"));

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


int WaIkPathLossSub(double** Raster, double** PathLoss, struct StructWaIk IniWaIk, char *area_type)
/*************************************************************************************************
 *
 *  Function WaIkPathLossSub calculates PathLoss in dB using Walfish-Ikegami path loss formula
 *    **PathLoss: array of path loss in dB
 *    **Raster:       input DEM file
 *  
 *    T.Javornik, Jan. 2010
 *
 *************************************************************************************************/
{
  // Walfish-Ikegami model constants and variables
  double BSxIndex = IniWaIk.BSxIndex;       //  normalized position of BS -> UTMx/resolution 
  double BSyIndex = IniWaIk.BSyIndex;       //  normalized position of BS -> UTMy/resolution
  double AntHeightBS = IniWaIk.BSAntHeight; //  Antenna height of BS [m]
  double AntHeightMS = IniWaIk.MSAntHeight; //  Antenna height of MS [m]
  int xN = IniWaIk.xN;                      //  dimension of the input(Raster) and output (PathLoss)
  int yN = IniWaIk.yN;                      //  dimension of the input(Raster) and output (PathLoss)
  double scale = IniWaIk.scale;             //  Resolution Walfish-Ikegami model

  double W0 = IniWaIk.W0;                 // Free space loss correction
  double W1 = IniWaIk.W1;                 // Reduced base antenna height correction
  double W2 = IniWaIk.W2;                 // Range correction
  double W3 = IniWaIk.W3;                 // Street width correction
  double W4 = IniWaIk.W4;                 // Frequency correction
  double W5 = IniWaIk.W5;                 // Building height correction
  double W6 = IniWaIk.W6;                 //  Street width [m]
  double W7 = IniWaIk.W7;                 //  Distance between buildings [m] 
  double W8  = IniWaIk.W8;                //  Building height [m] 
  double PHI_Street  = IniWaIk.PHIStreet; //  Street orientation
    
  double freq  = IniWaIk.freq;      //  carrier frequency
  double ResDist = IniWaIk.ResDist; //  distance BS - MS sampling rate [normalized with scale]
  double radi = IniWaIk.radi;       // radius of calculation

  //--------------------------------------------------------------------------------------------------------------------------------------------  
  double ZoBS;                    
  double ZoTransBS,ZoTransMS;              // BS and MS height about the sea level
  double ZoTransBS_delta, ZoTransMS_delta; // BS and MS height above/below roof
  double FreeSpacePathLoss = 0;            //  Free space path loss    
  double PathLoss_RTS;                     //  Roof-top-to-street diffraction and scatter loss
  double PathLoss_MSD;                     //  Multi-screen loss
  double PathLoss_Street;                  //  Street orientation loss 
    
  double PathLossTmp = 0; // tmp path loss

  double tiltBS2MS; // (ZoBS-ZoMS)/distBS2MSNorm    
    
  int ix; int iy; 
  double DiffX, DiffY; // Difference in X and Y direction
    
  double DistBS2MSNorm, DistBS2MSKm; // distance between MS and BS in Km sqrt(x2+y2+z2) * scale / 1000
                                     // normalized distance between MS and BS in xy plan sqrt(x2+y2)
  double ZObs2LOS = 0;
  double DistObs2BS = 0;

//  (I.O. 2-oct-2013)
//  ZoBS = Raster[(int)BSxIndex][(int)BSyIndex];    // BS height above the sea level calculated from raster DEM file
  ZoBS = Raster[(int)( BSxIndex + 0.5)][(int)( BSyIndex + 0.5)];

  ZoTransBS = ZoBS + AntHeightBS;      // BS transmitter height above the sea level
  ZoTransBS_delta =  AntHeightBS - W8; // BS transmitter height above the roof

  // Path loss MSD factos
  double PathLoss_MSD1 = 0;
  double ka =0;
  double kd =0;
  double kf =0;
    
  // PathLoss_MSD1 -------------------------------------------------------------------------------
  if (AntHeightBS > W8)
    PathLoss_MSD1 = -18*log10(1+ZoTransBS_delta);
  else if (AntHeightBS <= W8)
    PathLoss_MSD1 = 0;
  // kf
  if (strcmp(area_type, "metropolitan") == 0)
    kf = -4 + 1.5*((freq/925)-1);
  else if (strcmp(area_type, "medium_cities") == 0)
    kf = -4 + 0.7*((freq/925)-1);       
  else
    G_fatal_error(_("Unknown area type: [%s]."), area_type);
  //------------------------------------------------------------------------------------------------------
    
  for (ix = 0; ix < xN; ix++)
  {
    G_percent(ix, xN, 2);
    for (iy = 0; iy < yN; iy++)
    {
      DiffX = (BSxIndex-ix); DiffY = (BSyIndex-iy);
      ZoTransMS = Raster[ix][iy]+AntHeightMS; // ZoMS

      ZoTransMS_delta = W8 - AntHeightMS; // MS receiver hight below roof

      DistBS2MSKm = sqrt(DiffX*DiffX+DiffY*DiffY)*scale/1000;
      DistBS2MSNorm = sqrt(DiffX*DiffX+DiffY*DiffY);
            
      if (DistBS2MSKm < 0.01)
      {
        DistBS2MSKm = 0.01;
      }
      if ((DistBS2MSKm) > radi)
      {    
        continue;
      }

      // Calc position of the height and position of the highest obstacle
      tiltBS2MS= ZoTransBS - ZoTransMS; 

      if (DistBS2MSNorm > 0)
      {
        tiltBS2MS = -tiltBS2MS/(DistBS2MSNorm);
      }       
      else
      {
        tiltBS2MS = 0;
      }
      DoProfile( &ZObs2LOS, &DistObs2BS, ResDist, Raster, BSxIndex, BSyIndex, ZoTransBS, ix, iy, tiltBS2MS);        
      //PathLossTmp=0;

      if (ZObs2LOS < 0)   // **** LOS ****
      {
        PathLossTmp = 42.6 + 26*log10(DistBS2MSKm) + 20*log10(freq);
      }
      else if (ZObs2LOS >= 0) // **** NLOS ****
      {
        // Free space loss  ----------
        FreeSpacePathLoss = W0 + 20*log10(freq) + 20*log10(DistBS2MSKm);            
                
        // Rooftop-to-street difraction and scatter loss - PathLoss_RTS  ----------
        if (PHI_Street>=0 && PHI_Street < 35)
          PathLoss_Street = -10 + 0.354*PHI_Street;
        else if (PHI_Street>=35 && PHI_Street < 55)
          PathLoss_Street = 2.5 - 0.075*(PHI_Street - 35);
        else if (PHI_Street>=55 && PHI_Street < 91 )
          PathLoss_Street = 4 - 0.114*(PHI_Street - 55);
        else
        {
          PathLoss_Street = 0;  // Tx location (PHI_Street= nan)
        }   
        if (W8>AntHeightMS)
        {
          PathLoss_RTS = -8.2 - W3*log10(W6) + W4*log10(freq) + W5*log10(ZoTransMS_delta) + PathLoss_Street;  
        }
        else
        {
          PathLoss_RTS= 0;
        }

        // Multi-screen loss - PathLoss_MSD  ----------
        // ka
        if (DistBS2MSKm >=0.5 && AntHeightBS <= W8) 
          ka = W1 - 0.8*ZoTransBS_delta;
        else if (DistBS2MSKm < 0.5 && AntHeightBS <= W8)    
          ka = W1 - 0.8*ZoTransBS_delta*(DistBS2MSKm/0.5);
        else if (AntHeightBS > W8)  
          ka = W1;
        // kd
        if (AntHeightBS > W8)
          kd = W2; 
        else if (AntHeightBS <= W8)
          kd = W2 - 15*(ZoTransBS_delta/W8);
                    
        PathLoss_MSD = PathLoss_MSD1 + ka + kd*log10(DistBS2MSKm) + kf*log10(freq) - 9*log10(W7);
        if (PathLoss_MSD < 0)
        {
          PathLoss_MSD = 0;       
        }
        // ------------------------------------------------------------------------------------------------------------------------------
        PathLossTmp=FreeSpacePathLoss + PathLoss_RTS + PathLoss_MSD;
      }                      
      // write data to pathloss --------------------------------------------------------------------------------------------------
      PathLoss[ix][iy] = PathLossTmp; 
    } // end irow
  } // end icol
    return 0;
}

