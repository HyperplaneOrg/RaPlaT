/*******************************************************************
 *
 *  LTE_RaPlaT_Fun.h: header file for  LTE_RaPlaT_Fun.c
 *
 *  Tomaz Javornik (13.5.2014)
 *  - (17.9.2014) add LTE_NF and  LTE_INTERFERENCE_MARGIN (TJ)
 *
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *******************************************************************/

#define LTE_CINRSIZE 15
#define LTE_TABLEFACT 0.0001
#define LTE_BWMHZ 10.0       // default bandwidth in MHz
#define LTE_NRB 50           // default number of resource blocks
#define LTE_nPDCCH  1        // number of Physical downlink control channels
#define LTE_CPF 'n'          // normal cyclic prefix is default
#define LTE_OVERHEAD 0.048   // LTE overhead due to  PCHFICH, PHICH and PDCCH channles
#define LTE_WARNING_FLAG 1   // LTE worning flag, write warning on console
#define LTE_WARNING printf(" LTE warning! Default LTE values applied for calculations! \n")  // LTE ouput warning
#define LTE_NF 7.0                   // Noise figure of the receiver in dB
#define LTE_INTERFERENCE_MARGIN 3.0  // Interference margin in dB	


// spectral efficiency table in [bits/s/Hz] for LTE system
static int LTE_EFFICTABLE[LTE_CINRSIZE] = { 1523,  2344, 3770, 
                                            6016,  8770, 11758,
                                           14766, 19141, 24063,
                                           27305, 33223, 39023,
                                           45234, 51152, 55547};

// CINR for gassian channel in [dB]
static int LTE_CINRTABLE[LTE_CINRSIZE] = { -70000, -50714, -31429,
                                           -12143,   7143,  26429,
                                            45714,  65000,  84286,
                                           103571, 122857, 142143,
                                           161429, 180714, 200000};

