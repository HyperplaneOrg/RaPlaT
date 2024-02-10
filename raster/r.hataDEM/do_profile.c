/****************************************************************************
 *
 * MODULE:       do_profile function for r.hataDEM
 * AUTHOR(S):    Andrej Vilhar, Andrej Hrovat, Jozef Stefan Institute                
 *               Igor Ozimek (modifications & corrections), Jozef Stefan Institute
 *
 * PURPOSE:      Calculates radio coverage from a single base station 
 *               according to model 9999
 *
 * COPYRIGHT:    (C) 2009-2014 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/* History:
30-may-2013 (I.O.)
 - minor modifications
*/


#include "local_proto.h"


int DoProfile( double *ZObs2LOS, double *DistObs2BS, double ResDist, double** Raster, double xBS, double yBS, double ZoTransBS,
               double xMS, double yMS, double tiltBS2MS)
/***********************************************************************************************
 *
 *  Function finds the highest obstacle between [xBS,yBS,ZoTransBS] and calculates:
 *      *DistObs2BS:    nomalized distance to from BS to obstacle [normalized]
 *                      real distance (*DistObs2BS) * resolution 
 *      *ZObs2LOS:      distance between LOS line and obstacle [m]
 *                      *ZObs2LOS < 0:  obstacle is below LOS line
 *                      *ZObs2LOS > 0:  obstacle intersect LOS line
 *
 *      ResDist:        step for sample DEM [Normalized with resolution]
 *      Raster:         DEM array;
 *      xBS, yBS:       BS location [Normalized with distance]
 *      ixMS, iyMS:     Mobile station loction [Normalized with distance]
 *      ZoTransBS:      The transmitter elevation (height) above sea level DEM + AntHeightBS
 *      tiltBS2MS:      Tilt Between LOS line BS-MS and horizontal line: 
 *                      (ZoTransBS-ZoTransMS)/(DistBS2MS) 
 *
 *      T.Javornik, Jan. 2010
 *
 *************************************************************************************************/
{
  double DistBS2MS;
  double cols, rows;
  double AZI, X, Y;
  double DistTmp = 0;
  double HeightTmp = 0;
  double yTmp, xTmp;
  int iX, iY;

  // output variables initialisation
  *DistObs2BS = 0;
  *ZObs2LOS = -10000.0;

//  iX = (int)xBS; iY = (int)yBS;
  cols = xBS - xMS;
  rows = yBS - yMS;
  DistBS2MS = sqrt( cols * cols + rows * rows);

  AZI = atan( cols / rows);
  X = ResDist * sin( AZI); xTmp = xBS;
  Y = ResDist * cos( AZI); yTmp = yBS;
  if ( Y < 0 ) Y = Y * -1.0;
  if ( X < 0 ) X = X * -1.0;
  if ( rows < 0  && cols < 0 ) { X = X;  Y = Y; }  // NE quadrant
  if ( rows >= 0 && cols < 0 ) { X = X;  Y = -Y;}  // SE quadrant
  if ( rows >= 0 && cols >= 0) { X = -X; Y = -Y;}  // SW quadrant
  if ( rows < 0  && cols >= 0) { X = -X; Y = Y; }  // NW quadrant

  DistBS2MS = DistBS2MS - ResDist;
  while ( DistTmp < DistBS2MS)
  {
    DistTmp = DistTmp + ResDist; xTmp = xTmp + X; yTmp = yTmp + Y;   // calculate position in raster
    iX = (int)( xTmp + 0.5); iY = (int)( yTmp + 0.5);                // index of raster
    HeightTmp = -DistTmp * tiltBS2MS - ZoTransBS + Raster[ iX][ iY]; // calculate relative height above the LOS

    if ( HeightTmp > *ZObs2LOS)
    {
      // estimate max height above the LOS
      *ZObs2LOS = HeightTmp;
      *DistObs2BS = DistTmp;
    }
  }  // end while

  return 0;
}  // end doProfile

