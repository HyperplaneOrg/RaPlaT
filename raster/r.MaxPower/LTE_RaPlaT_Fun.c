/*******************************************************************
 *
 * LTE_RaPlaT_Fun: program for testing RaPlat LTE routnes
 *
 * Tomaz Javornik (13.5.2014)
 * modifications and r.MaxPower integration: Igor Ozimek (avg. 2014)
 *
 * v.1.0 21.08.2014
 *       - added option -o i (suborutine output is interference
 *       . change input signal sigMax: signal from best server
 *                             sigSum: sum of signals from all BS
 * v.1.1 16.09.2014
 *       - added impact of receiver noise figure NF = 7 dB (LTE_NF)
 *       - added impact of interference margin = 3 dB (LTE_INTERFERENCE_MARGIN)
 *
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *******************************************************************/

#include "stdafx.h"


/***********************************************************************/
int BwMHz2nRB( double *BwMHz) 
/***********************************************************************
*
*  Function returns the number of Resource blocks in LTE systems
*    BwMHz: LTE bandwidth in MHz
*    nRB: number of resource blocks
*    if BwMHz is not LTE bandwidth, default value of RB is set 
*
*************************************************************************/
{
  int nRB = 0;
  int iBwMHz;
  iBwMHz = (int)( *BwMHz) * 10;
  switch ( iBwMHz)
  {
  case 14:  nRB = 6;   break;
  case 30:  nRB = 15;  break;
  case 50:  nRB = 25;  break;
  case 100: nRB = 50;  break;
  case 150: nRB = 75;  break;
  case 200: nRB = 100; break;
  default:
    nRB = LTE_NRB;
    *BwMHz = LTE_BWMHZ;
    if ( LTE_WARNING_FLAG)
    {
      LTE_WARNING;
      printf("   Error in number of resource blocks! This is not a LTE Bandwidth! \n");
      printf("   Default number of RB is set: %5i \n", nRB); 
    }
    break;
  }
  return nRB;
};


/***********************************************************************/
double LTEOverHead( double *BwMHz, int *nPDCCH, int *AntennaNum, char *cpf)
/*******************************************************************
*
*  Function returns overhead of LTE system in downlink
*    BwMHz: LTE bandwidth in MHz
*    CPF: cyclic prefix type 'n' = normal, 'e' = extended
*    nPDCCH: number of PDCCH (1,2,3, or 2,3,4 for BwMHz = 1.4) 
*    AntennaNum: number of antennas
*
*    LTE in Bullets, Table 34 and Table 59
*    Tomaz Javornik (26.11.2013)           
*
*******************************************************************/
{
  double OverHead = -1.0;
  double OverHeadadd = -1.0;

  if ( ( *cpf != 'n') && ( *cpf != 'e'))
  {
    if ( LTE_WARNING_FLAG) LTE_WARNING;
    if ( LTE_WARNING_FLAG) printf("    Error in Cyclic prefix! The normal Cyclic prefix is used! \n");
    *cpf = 'n';
  }

  // calculates the overhead due to PCHFICH, PHICH and PDCCH channles
  if ( *nPDCCH == 4)
  {
    if ( (int)( *BwMHz * 10) != 14)
    {
      *BwMHz = 1.4;
      if ( LTE_WARNING_FLAG)
      {
        LTE_WARNING;
        printf("    Error in LTEOverHead! \n");
        printf("    Channel bandwidth or number of PDCCH are not correct! \n");
        printf("    Setting bandwidth [MHz] = %f \n", *BwMHz);
      }
    }
    OverHead = 0.264;
    if (*cpf == 'e') OverHead = 0.278;
  }
    
  switch (*nPDCCH)
  {
    case 1:
      OverHead = 0.048;
      if ( *cpf == 'e') OverHead = 0.056; 
      break;       
    case 2:
      OverHead = 0.119;
      if ( *cpf == 'e') OverHead = 0.139; 
      break;       
    case 3:
      OverHead = 0.190;
      if ( *cpf == 'e') OverHead = 0.222; 
      break;
    case 4:
      break;
    default:
      *nPDCCH = LTE_nPDCCH;
      OverHead  = LTE_OVERHEAD;
      if ( LTE_WARNING_FLAG)
      {
        LTE_WARNING;
        printf("    Error in LTEOverHead! \n");
        printf("    Number of PDCCH is not correct! Default value is set: %i \n", *nPDCCH);
      }
      break;
  }

  switch (*AntennaNum)
  { // Overhead due to reference signals
    case 1:
      OverHeadadd = 0.048;
      if ( *cpf == 'e') OverHeadadd = 0.056;
      break;
    case 2:
      OverHeadadd = 0.095; 
      if ( *cpf == 'e') OverHeadadd = 0.111; 
      break;     
    case 4:
      OverHeadadd = 0.143;
      if ( *cpf == 'e') OverHeadadd = 0.167; 
      break;         
    default:
      *AntennaNum = 1;
      OverHeadadd = 0.048;
      if ( *cpf == 'e') OverHeadadd = 0.056;
      if ( LTE_WARNING_FLAG)
      {
        LTE_WARNING;
        printf("    Error in LTEOverHead! \n");
        printf("    Number of antennas is not correct! Using values for 1 Antenna! \n");
      }
      break;
  }

  OverHead = OverHead + OverHeadadd;
  return OverHead;
}


/***********************************************************************/
double LTEOverHeadApp( double *BwMHz, int *nPDCCH, int *AntennaNum, char *cpf)
/*******************************************************************
*
*  Function returns overhead of LTE system in downlink
*    BwMHz: LTE bandwidth in MHz
*    CPF: cyclic prefix type 'n' = normal, 'e' = extended
*    nPDCCH: number of PDCCH (1,2,3, or 2,3,4 for BwMHz = 1.4) 
*    AntennaNum: number of antennas
*
*    LTE in Bullets, Table 61
*    Tomaz Javornik (28.07.2014)           
*
*******************************************************************/
{
  double OverHead = -1.0;
  int iOverHead = -1;
  int i = 0;
  int j = 0; 

  int tA[4][6] = { {1000, 870, 879, 886, 888, 890}, {765, 799, 808, 815, 817, 818}, {694, 728, 737, 743, 746, 747}, {623, 1000, 1000, 1000, 1000, 1000}};
  int tB[4][6] = { {1000, 849, 860, 867, 870, 871}, {728, 766, 776, 784, 787, 788}, {644, 683, 693, 701, 703, 705}, {575, 1000, 1000, 1000, 1000, 1000}};
  int tC[4][6] = { {1000, 835, 844, 851, 853, 854}, {731, 764, 773, 779, 781, 782}, {660, 692, 701, 708, 710, 711}, {588, 1000, 1000, 1000, 1000, 1000}};
  int tD[4][6] = { {1000, 809, 819, 826, 828, 830}, {689, 726, 735, 743, 745, 746}, {606, 642, 652, 659, 662, 663}, {550, 1000, 1000, 1000, 1000, 1000}};


  if ( ( *cpf != 'n') && ( *cpf != 'e'))
  {
    *cpf = 'n';
    if ( LTE_WARNING_FLAG)
    {
      LTE_WARNING;
      printf("    Error in Cyclic prefix! The normal Cyclic prefix is used! \n");
    }
  }

  if ( ( *nPDCCH > 4) || ( *nPDCCH < 1))
  {
    *nPDCCH = 2;
    if ( LTE_WARNING_FLAG)
    {
      LTE_WARNING;
      printf("    Error in number of PDDCH symbols! The number of PDCCH symbols is set to 2! \n");
    }
  }

  if ( (*AntennaNum > 2) || ( *AntennaNum < 1))
  {
    *AntennaNum = 1;
    if ( LTE_WARNING_FLAG)
    {
      LTE_WARNING;
      if(LTE_WARNING_FLAG) printf("    Error in number Tx and Rx antenna! The number of Tx and Rx antennas is set to 1! \n");
    }
  }

  // calculate the PDSCH overhead for applications
  if ( *nPDCCH == 4)
  {
    if( (int)( *BwMHz * 10) != 14)
    {
      *BwMHz = 1.4;
      if ( LTE_WARNING_FLAG)
      {
        LTE_WARNING;
        printf("    Error in LTEOverHead! \n");
        printf("    Channel bandwidth or number of PDCCH are not correct! \n");
        printf("    Setting bandwidth [MHz] = %f \n", *BwMHz);
      }
    }
  }
    
  switch ( (int)( *BwMHz * 10))
  {
    case 14:  i = 0; break;
    case 30:  i = 1; break;
    case 50:  i = 2; break;
    case 100: i = 3; break;
    case 150: i = 4; break;
    case 200: i = 5; break;
    default:  i = 2; break;
  };

  j = *nPDCCH - 1;
  switch ( *AntennaNum)
  { // Overhead due to reference signals
    case 2:
      iOverHead = tC[j][i];
      if ( *cpf == 'e') iOverHead = tD[j][i];
      break;
    default:
      iOverHead = tA[j][i];
      if ( *cpf == 'e') iOverHead = tB[j][i];
      break;
  };

  if ( iOverHead == 1000) iOverHead = 879;
  OverHead = (double)( 1000 - iOverHead) / 1000.; 
  return OverHead;
}


/***********************************************************************/
int PdBm2LteThroughput( int nRows, int nCols, float *sigMax, float *sigOut, float *sigSum, char ChanType, char OutputFlag, double *BwMHz, int *nPDCCH, int *nAntenna, char *cpf)
/*******************************************************************
*
* Function returns:
*   RSPR - Reference signal received power 
*   max. spectral efficiencty assuming no interference
*   max. throughput [bit/s]
* 
* Parameters:
*   nRows:         number of Rows in sigMax & sigOut
*   nCols:         number of Colomns in sigMax & sigOut
*   sigMax:        received signal power from best serving Base Station [dBm]
*   sigSum:        sum of power from all base stations [dBm]; 
                   RSSI without noise (RSSI - No)
*   sigOut:        output matrix
*   ChanType:      channel type 'g' Gaussian, 'r' Rayleigh
*   OutputFlag:    determines the output of the calculation
*     'p' --> RSRP (received signal representative power)
*     'r' --> RSSI (received signal strenght)
*     'q' --> RSRQ (received signal representative quality)
*     'c' --> CINR (max CINR, interference free)
*     's' --> max. spectral efficiency considering only AWGN
*     't' --> max. throughput
*     'i' --> inteference in (dBm) 
*   BwMHz:         bandwidth in MHz
*   nPDCCH:        number of physical downlink control channel
*   cpf:           cyclix prefix type 'n' normal 'e' extended
*   nAntenna:      number of transmit antennas
*
*                 10% overhead due to retransmission
*                 5% additional overhead
*
* Tomaz Javornik (13.5.2014)
*
*******************************************************************/
{
  int ReturnVal = -1;                      // procedure output 
  int nRB = 0;                             // number of resoruce blocks
  int nRE = 84;                            // number of resoruce elements 
  double SpecEff2ThroughPut = 1;           // factor which converts Spectraefficency per bin into throughput
  double NodBm = -132.07;                  // AWGN noise power in dBm for 15 kHz BW in one RE
                                           // 10*log10(k*B*T) + 30 = 10*log10(1.380e-23*300*15e3) + 30
  double NoBwmW;                           // noise in entire Bw in mW
  double InterfmW;                         // interference in mW 
  double OverHead;
  double tmpSNR;                           // signal to noise (interference ratio)
  double log12nRB;                         // factor 10.0 * log10( nRB * 12.0) 
 
  double tmpdouble;
  double tmpRSSI;                          // temporal value of RSSI
  int i, j, k, tmpInt;

  nRB = BwMHz2nRB( BwMHz);                // number of resource blocks calculation
  OverHead = LTEOverHeadApp( BwMHz, nPDCCH, nAntenna, cpf);  // overhead calculations 
  if ( *cpf == 'e') nRE = 72;
  nRE = (int)( nRE / 0.5e-3);
  log12nRB = 10.0 * log10( nRB * 12.0);

  NoBwmW = 12.0 * nRB * pow( 10, 0.1 * (NodBm + LTE_NF));  // noise in mW entire bandwdith & receiver generated noise

  // constant which convert the spectral efficiency to throughput in Mbit/s assuming 10% retransmission and 5% additonal overhead
  SpecEff2ThroughPut = SpecEff2ThroughPut * nRB * 180.0e3 * ( 1.0 - OverHead) / ( 1.0e6 *  1.10 * 1.05); 

  long ij;
  for ( i = 0; i < nCols; i++)
  {
    for ( j = 0; j < nRows; j++)
    {
      ij = i * nRows + j;
      
      // default option (-p)	
      // RSRP = RSSI (without noise and intererence) - 10 * log10( 12.0 * nRB) 
      // RSRP = sigMax[ij] - 10 * log10( 12.0 * nRB) = sigMax[ij] - log12nRB
      // RSRP = -44:-140 dBm
      sigOut[ij] = sigMax[ij] - log12nRB;  // RSRP is default output
      tmpRSSI = pow( 10.0, 0.1 * sigSum[ij]) + NoBwmW;  // add noise power
      tmpRSSI = 10.0 * log10(tmpRSSI);                  // covert to dBm
      tmpSNR = sigOut[ij] + log12nRB - 10.0 * log10( NoBwmW);
      InterfmW = pow( 10.0, 0.1 * sigSum[ij]) - pow( 10.0, 0.1 * sigMax[ij]);
      
      if( fabs( sigMax[ij] - sigSum[ij]) < 0.0001) InterfmW = FLT_MIN;		
      if( InterfmW <= 0.0) InterfmW = FLT_MIN;
      if( OutputFlag == 'p')   // RSRP (received signal representative quality)               
      {  
        if( sigOut[ij] < -140.0) sigOut[ij] = -140.0;
        if( sigOut[ij] > -44.0) sigOut[ij] = -44.0;
      }
		

      // (-r) RSSI = sum of power od all signals + noise power		
      if( OutputFlag == 'r')   // -r RSSI (received signal strenght)   
      {  
        sigOut[ij] = tmpRSSI;
      }

      // RSRQ (-q) Interference
      // RSRQ = 10 * log10( nRB) + RSRP - RSSI; nRB is number of used resource blocks
      // RSRQ = 10.0 * log10( nRB) + sigOut[ij] - tmpRSSI
      // RSRQ = -19.5 : -3 dB      	
      if( OutputFlag == 'q')   // RSRQ (received signal representative quality)               
      {  
        sigOut[ij] = 10.0 * log10( nRB) + sigOut[ij] - tmpRSSI;
        if( sigOut[ij] < -19.5) sigOut[ij] = -19.5;
        if( sigOut[ij] > -3.0) sigOut[ij] = -3.0;
      }
  
      // value of interfering signla in dBm (i)
      if( OutputFlag == 'i')   // RSRQ (received signal representative quality)               
      {
        if ( InterfmW == FLT_MIN) sigOut[ij] = DB_MIN_VAL;
        else
        {
          sigOut[ij] = 10.0 * log10( InterfmW);
          if ( sigOut[ij] < DB_MIN_VAL) sigOut[ij] = DB_MIN_VAL;
        }
      }

      // CINR = RSRP + 10 * log10( 12.0 * nRB) - 10 * log10( Interf + Noise) (over all bandwidth));
      if( OutputFlag == 'c')   // max CINR / no iterference 
      {  
        sigOut[ij] = tmpSNR; 
      }

      
      if( OutputFlag == 's' || OutputFlag == 't' ) // max. spectral efficiency or throughput
      {  
        k = -1;
        tmpdouble = (tmpSNR - LTE_INTERFERENCE_MARGIN) / LTE_TABLEFACT;
        tmpInt = (int)floorf( tmpdouble );
        while( tmpInt > LTE_CINRTABLE[ k + 1])
        {
          k = k +1;
          if ( k > LTE_CINRSIZE) { k = LTE_CINRSIZE - 1; break;}
        }
        sigOut[ij] = 0;
        if( k > -1)
        { 
          sigOut[ij] = LTE_EFFICTABLE[k] * LTE_TABLEFACT;  // max. spectral efficiency - considering only AWGN
          if ( OutputFlag == 't')
          {  // switch 't'; 
            sigOut[ij] = sigOut[ij] * SpecEff2ThroughPut;  // max. throughput in Mbit/s
          }
        }
      }
    }
  }

  ReturnVal = 0;
  return ReturnVal;
}


/***********************************************************************
*
* Run example: ./Run.exe -o t -i -80 -b 5.0
*
***********************************************************************/

#if 0

int main( int argc, char *argv[] )
{
  printf("-----------------------------------------------------\n");
  printf("|                                                   |\n");
  printf("|      Program for testing LTE RaPlaT routines.     |\n");
  printf("|                                                   |\n");
  printf("-----------------------------------------------------\n");


  float *sigMax;                          // received signal from best server [dBm]
  float *sigOut;                          // results of calculation [dBm]
  float *sigSum;                          // sum of signals from all base station in [dBm]
  int nRows = 3;                          // Number of rows  
  int nCols = 5;                          // Number of columns
  char ChanType = 'g';                    // Channel Type
  char OutputFlag = 'p';                  // Type of output 's' spectral efficieny, 't' maximal throughput per bin
  double BwMHz = 5.0;                     // Bandwidth in MHz
  int nAntenna = 1;                       // Number of transmit antennas
  char cpf = 'n';                         // cyclic prefix, 'n' normal, 'e' extended
  int nPDCCH = 1;                         // number of physical downlink control channels  
  float sigMaxIni = -112.0;
  float sigInter =  -112.0;         

  int i,j,ij;
  float tmp;
     
  i = 1;  
  while (i < argc) 
  {
    if(strcmp(argv[i],"-h") == 0) 
    {                      
      printf(" Input arguments help: \n");
      printf(" -h help \n");
      printf(" -b channel bandwidth in MHz (default = -5 MHz) \n");
      printf(" -o ouput type:  c = CINR, s = Spec. eff., t = Max.throughput, \n");
      printf("                 r = RSSI, p =  RSRP, q = RSRQ \n");
      printf(" -i initial RSSI (default = -112 dBm) \n");
      return -1;    
    }
    if(strcmp(argv[i],"-b") == 0) 
    {                      
      i = i + 1;
      BwMHz = atof(argv[i]);  
    }
    if(strcmp(argv[i],"-o") == 0) 
    {                      
      i = i + 1;
      OutputFlag = argv[i][0];  
    }
    if(strcmp(argv[i],"-i") == 0) 
    {                      
      i = i + 1;
        sigMaxIni = atof(argv[i]);  
    }
    i = i + 1;                          
  }


  // Memory allocation
  sigMax = (float *)malloc(nCols*nRows*sizeof(float));
  for (i = 0; i < nCols*nRows; i++)   sigMax[i] = -999.0;

  sigOut = (float *)malloc(nCols*nRows*sizeof(float));  
  for (i = 0; i < nCols*nRows; i++) sigOut[i] = -999.0;

  sigMax[0] = sigMaxIni;
  sigMax[1] = sigMaxIni + 2;
  sigMax[2] = sigMaxIni + 4;
  sigMax[3] = sigMaxIni + 6;
  sigMax[4] = sigMaxIni + 8;
  sigMax[5] = sigMaxIni + 10;
  sigMax[6] = sigMaxIni + 12;
  sigMax[7] = sigMaxIni + 14;
  sigMax[8] = sigMaxIni + 16;
  sigMax[9] = sigMaxIni + 18;
  sigMax[10] = sigMaxIni + 20;
  sigMax[11] = sigMaxIni + 22;
  sigMax[12] = sigMaxIni + 24;
  sigMax[13] = sigMaxIni + 26;
  sigMax[14] = sigMaxIni + 28;
  
  sigSum = (float *)malloc(nCols*nRows*sizeof(float));
  tmp = pow( 10.0, 0.1 * (sigInter));
  for (i = 0; i < nCols*nRows; i++)
  {
     sigSum[i] = pow( 10.0, 0.1 * sigMax[i])  + tmp;
     sigSum[i] = 10 * log10( sigSum[i]); 
  }

  PdBm2LteThroughput(nRows,nCols, sigMax, sigOut, sigSum, ChanType, OutputFlag, &BwMHz, &nPDCCH, &nAntenna, &cpf);

  // Test print
  printf("\n");
  printf("-----------------------------------------------------\n");
  printf("|       Test Output!                                |\n");
  printf("-----------------------------------------------------\n");
  printf("\n");
  printf(" Inteference [dBm] = %f \n", sigInter);
  printf("\n");
  printf("----- Best server =  RSSI - (No + Interference) -----\n");
  for(j=0;j<nRows;j++) 
  {
    for(i=0;i<nCols;i++) 
    {
      ij = i * nRows + j;
      printf(" %3.1f ", sigMax[ij]); 
    }
    printf(" \n");
  }
  printf("\n");

  switch (OutputFlag) 
  {
    case 'c':
      printf("-------------  CINR [dB]   ----------\n");
      break;

    case 's':
      printf("-------  Spec. Eff. [bits/s/Hz] -----\n");
      break;

    case 't':
      printf("------ Max throughput  [Mbit/s] -----\n");
      break;

    case 'q':
      printf("---------------  RSRQ [-19.5 : -3] [dB]  ---------\n");
      printf("----------  RSRQ [-10.8] [dB]  noise limited -----\n");
      break;

    case 'r':
      printf("--------------  RSSI [dBm] ----------\n");
      break;

   case 'i':
      printf("--------------  Interference [dBm] ----------\n");
      break;


    default:
      printf("--------------  RSRP [-44:-140][dBm] ----------\n");
      break;
  }

  for( j = 0 ; j < nRows ; j++) 
  {
    for( i = 0 ; i < nCols; i++) 
    {
      ij = i * nRows + j;
      printf(" %6.3f ", sigOut[ij]); }
      printf(" \n");
    }
    printf("\n");   
  }

#endif

