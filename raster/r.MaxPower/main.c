/****************************************************************************
 *
 * MODULE:       r.MaxPower
 * AUTHOR(S):    Tine Celcer, Jozef Stefan Institute
 *               Igor Ozimek, JSI (a major rewrite to improve performance - speed and storage space)
 *
 * PURPOSE:      Sorts the recieved power levels from different cells in descending order,
 *               writes data in the form of a GRASS map and as a database table
 *             
 * COPYRIGHT:    (C) 2009-2018 Jozef Stefan Institute
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with RaPlaT
 *               for details.
 *
 *****************************************************************************/

/*
GRASS GIS 7 module, tested on Ubuntu 18.04

History (I.O.):

30-nov-2011 - a major rewrite, mostly of fill_database for faster execution of database inserts
 - explicit column names not used for inserts (SQL)
 - multirow inserts (SQL) - additional parameter dbperf
 - START TRANSACTION / COMMIT (speed-up for transactional data tables PostgreSQL, MySQL with InnoDB)
 - fast insert using intermediate CSV file and LOAD DATA LOCAL INFILE (MySQL) or COPY (PostgreSQL)
 - rearranged (unified) C style (wrong idents - tab/space confusion)
 - skipping null values (conditional compilation - SKIPNULL constant defiend below)
 - rearrange output map closing (before data table writing)
 - region confusion resolved: current region should be used (and not the region of the first sector map)
 - various minor corrections & modifications
 - included funcionality of db.GenerateTable (creation of a new empty table)
 - workaround for detecting already exisiting pg database (db_table_exists() doesn't really work for PostgreSQL)
 - GRASS' OVEWRITE environment variable used to allow database overwrite, instead of special -o flag 

09-dec-2011
Changes related to the new antenna ID field
 - changed format of the 'cell_input' file - antenna ID parameter inserted at the second place (after cell name)
 - changed data table format - 'id' column added (descriptive 'cell' and 'model' columns retained for now)
   ('model' column too short to fit all the parameters of some new models - will be cut short)
Other changes:
 - ';;' before '\n' not required any more in the intermediate (temporary) radio cell/sector table file
 - workaround for PostgreSQL (see above) not needed any more... (a transient problem ???)

17-dec-2011
 - maximum sizes of "cell" an "model" data table columns (strings) enlarged to 32 and 128 respectivelly (were 10 and 35)
 - size of "ID" data table enlarged from smallint to integer

25-dec-2011
 - PostgreSQL workaround is back... :(

21-feb-2012
  - output data table's columns idN and EcN0 are now stored as REAL (not INTEGER x100) - 4 bytes required as before
  - the first line of the input radio cell/sector table file (containing the number of the following lines) is removed
    (not needed any more)
  - some other modifications and code cleaning

17-oct-12:
 - additional input parameter rx_threshold (areas with received dBm less then that are treated as without coverage)
   (it is an optional parameter - backward compatibility is maintained)
 - flag -1 added: Rx (dBm) values in output map replaced by 1.0 when above rx_threshold

02-oct-2013
 - modified Description (includes version)
 - pixel coordinates in database now alligned to the centers of pixel squares (instead their WN corners)

4-jul-2014
 - added support for CSV-format output file instead of database output, in the form of a quasi database engine 'csv';
   the file written is actually the one normally created as a temporary intermidiate file for fast mode MySQl/PostgreSQL
   data tabel creation (when dbperf=99), it contains no header line

29-jul-2014
 - added switch -s: the output raster map contains the sum of all signal instead of the strongest received signal at each
   raster point
 - some minor corrections

25-aug-2014
 - large arrays for input/output maps (raster images) and related data are reorganized: input map, output strongest signals
   (one or more maps, from one or more transmiiters, the first one is "best server" masp), best server indeces map, sum of
   all signals map);
   changes (to remove the mess - 2D arrays have been used for both 2D and 3D data, to simplify/speed-up some processing, and
   to use less memory space):
    - individual maps are stored in 1-dimensional arrays, multiple output strongest signals are stored in 2-dimensional arrays
      with the slow-running dimension for the map index
    - all maps arrays are 'float' type now (sum of signals was 'double' before)
 - added additional LTE-related computations (implemented in the new LTE_RaPlaT_Fun.c module)
 - added r.MaxPower parameter 'generate' for specifying computations (old and the new additional LTE-related ones)
   'generate=rss-max' (default) generates the old output map,
   'generate=coverage' generates the old '-1' (flag) map (single color map)
   'generate=rss-sum' generates the old '-s' map (sum of signals map)
 - flags '-1' and '-s' removed (superseded by 'generate=coverage' and 'generate=rss-sum', respectivelly)
 - added r.MaxPower parameter 'bandwidth' (needed for the new LTE-related computations)

23-sep-2014
 - updated version of LTE_RaPlaT_Fun.c (v_1.1)

21-jul-2017 (I.O.)
 - modified for GRASS GIS 7 (was GRASS GIS 6)
     GRASS API (function names, ...)
 - minor corrections (to remove compiler warnings, redundant code, etc.)

24-aug-2017 (I.O.)
 - added support for sqlite database

7-dec-2018
 - minor changes to avoid C compiler warnings

23-mar-2021
 - removal of the leftovers of the -1 flag
*/


#define VERSION "v24mar2021"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>
#include <grass/dbmi.h>


#include "common.h" // defines DB_MIN_VAL

#define GENERATE_BASIC "rss-max,coverage,rss-sum,rss-maxix"
#define GENERATE_LTE "lte-rssi,lte-rsrp,lte-rsrq,lte-cinr,lte-maxspecteff,lte-maxthrput,lte-interfere"
#define GENERATE_default "rss-max"

#define SKIPNULL 1 /* 0: no skip, 1: skip (do not invlude null points in the output data table) */


int create_table ( char *drv_name, char *db_name, char *tbl_name, int cell_num, int ovr)
{
  dbConnection conn;
  dbDriver *driver; 
  dbTable *table;
  int i,ncols_tab;
  char temp_str[100];  
  char temp_str1[100];

  ncols_tab=4*cell_num+4;   


  if (db_table_exists(drv_name, db_name, tbl_name) == 1)
  { 
    if (ovr)
    {
      if ( db_delete_table(drv_name, db_name,tbl_name) != DB_OK)
        G_fatal_error(_("Cannot delete/overwrite existing table <%s>"), tbl_name);
    }
    else G_fatal_error(_("Table <%s> already exists"), tbl_name);
  }


  /* set connection */
  db_get_connection(&conn);  /* read current */
  if (drv_name) conn.driverName = drv_name;
  if (db_name) conn.databaseName = db_name;
  db_set_connection(&conn);  


  driver = db_start_driver_open_database(drv_name, db_name);
  

#if 1  // PostgreSQL workaround  

  // workaround for PostgreSQl - db_table_exists() (above) simply doesn't find the existing table

  dbString sql[100]; // for PostgreSQL workaround
//  dbString table_name[100]; // for PostgreSQL workaround

  if (strcmp("pg",drv_name) == 0)
  {
    if (ovr)
    {
      G_message(_("\nDeleting table <%s> (if exists)"), tbl_name); 
      db_set_string( sql, "DROP TABLE IF EXISTS ");
      db_append_string( sql, tbl_name);
      if (db_execute_immediate( driver, sql) != DB_OK)
        G_fatal_error(_("Failed dropping table (PostgreSQL)! "));
    }
/*
    else
    {
      db_init_string(&table_name);
      db_set_string(&table_name, tbl_name);    
      if( db_describe_table(driver,&table_name,&table) == DB_OK)
        G_fatal_error(_("Table <%s> already exists"), tbl_name);
    }
*/
  }

#endif  // PostgreSQL workaround

  table = db_alloc_table(ncols_tab);
  db_set_table_name(table, tbl_name);


  db_set_column_name(&table->columns[0],"x");
  db_set_column_sqltype(&table->columns[0],DB_SQL_TYPE_INTEGER);
  db_set_column_length(&table->columns[0],6);

  db_set_column_name(&table->columns[1],"y");
  db_set_column_sqltype(&table->columns[1],DB_SQL_TYPE_INTEGER);
  db_set_column_length(&table->columns[1],6);
 
  db_set_column_name(&table->columns[2],"resolution");
  db_set_column_sqltype(&table->columns[2],DB_SQL_TYPE_INTEGER);
  db_set_column_length(&table->columns[2],4);

  int c_ix = 3;
  for (i = 1; i < cell_num + 1; i++)
  {
    strcpy(temp_str,"cell");
    sprintf(temp_str1,"%d",i);
    strcat(temp_str,temp_str1);
    db_set_column_name(&table->columns[c_ix],temp_str);
    db_set_column_sqltype(&table->columns[c_ix],DB_SQL_TYPE_CHARACTER);
    db_set_column_length(&table->columns[c_ix++],32);

    strcpy(temp_str,"id");
    sprintf(temp_str1,"%d",i);
    strcat(temp_str,temp_str1);
    db_set_column_name(&table->columns[c_ix],temp_str);
    db_set_column_sqltype(&table->columns[c_ix],DB_SQL_TYPE_INTEGER);
    db_set_column_length(&table->columns[c_ix++],6);

    strcpy(temp_str,"Pr");
    strcat(temp_str,temp_str1);
    db_set_column_name(&table->columns[c_ix],temp_str);
    db_set_column_sqltype(&table->columns[c_ix],DB_SQL_TYPE_REAL);
    db_set_column_length(&table->columns[c_ix++],6);

    strcpy(temp_str,"model");
    strcat(temp_str,temp_str1);
    db_set_column_name(&table->columns[c_ix],temp_str);
    db_set_column_sqltype(&table->columns[c_ix],DB_SQL_TYPE_CHARACTER);
    db_set_column_length(&table->columns[c_ix++],128);
  }

  db_set_column_name(&table->columns[ncols_tab-1],"EcN0");
  db_set_column_sqltype(&table->columns[ncols_tab-1],DB_SQL_TYPE_REAL);
  db_set_column_length(&table->columns[ncols_tab-1],6);   
       
  if ( db_create_table(driver,table) != DB_OK)
    G_fatal_error(_("Cannot create table <%s>! "), tbl_name);
  
  db_close_database(driver);
  db_shutdown_driver(driver);

  return DB_OK;
}

//------------------------------------------------------------------------------


int fill_database(char *drv_name, char *db_name, char *tbl_name, int db_perf, int ncols, int nrows, int x_start, int y_start, int res, float **arr_power, int **arr_index, char *cell_name[], int *antenna_id, char *model_name[], int ncells, float *arr_EcNo, int verbose, int ovr)
{
  dbDriver *driver = NULL;
  int col,row;
  int cell;
  int x,y;
  char buf[500];
  dbString sql[500];

  int csv2db;
  int SQLPacketSize;
  char csv_filename[20]="/tmp/fileXXXXXX"; //filename template for the temporaray CSV file
  FILE *csv_fp;
  int csv_fd;

  int csv_only; // =1 if we only want to create csv file without creating data table; =0 otherwise


  if ( strcmp( "csv", drv_name) == 0) csv_only = 1;
  else csv_only = 0;


  // only db_perf = 1 is allowed for dbf and sqlite database drivers
  if ( strcmp( "dbf", drv_name) == 0 || strcmp( "sqlite", drv_name) == 0) db_perf = 1;


  //open database, start driver, open table
  if ( !csv_only)
  {
    driver = db_start_driver_open_database(drv_name, db_name);

    if (db_perf >= 0 && db_perf < 99) 
    {
      csv2db = 0;
      SQLPacketSize = db_perf;
    }
    else
    {
      csv2db = 1;
      SQLPacketSize = 1;
    }
  }
  else csv2db = 1;


  db_init_string(sql);


  if (!csv2db)
  {
    // START TRANSACTION
//    db_set_string(sql,"START TRANSACTION");
//    if (db_execute_immediate(driver,sql) != DB_OK)
    if (db_begin_transaction (driver) != DB_OK) 
      G_fatal_error(_("Failed starting transaction! "));
  }
  else
  {
    if ( !csv_only)
    {
      // Temporary file for the intermidiate CSV file (to be read by the database manager MySQl or PostgreSQL)
      csv_fd = mkstemp(csv_filename);
      // PostgreSQL requires read access for other users
      fchmod(csv_fd, 0644); // PostgresSQL requires "other users" read access
      csv_fp = fdopen(csv_fd, "wb+");
      G_message(_("Creating indermediate temporary CSV file '%s'..."), csv_filename);
    }
    else
    {
      csv_fp = fopen( tbl_name, "r");
      if ( csv_fp)
      {
        fclose( csv_fp);
        if ( !ovr) G_fatal_error(_("Output csv file already exists! "));
      }
      csv_fp = fopen( tbl_name, "wb+");
      if ( !csv_fp)
        G_fatal_error(_("Cannot open output csv file for writing! "));
    }
  }


  long arr_ix;
  for (row = 0; row < nrows; row++)
  {
    y = y_start - row*res; 
    if (verbose) G_percent(row, nrows, 2);
     
    int SQLrowcnt = 0; // counts the rows while building a multirow SQL INSERT query
    for (col = 0; col < ncols; col++) 
    {

//      if (verbose) G_percent(row*ncols + col, nrows*ncols, 2);

      arr_ix = row * ncols + col;

      if( !SKIPNULL || arr_power[ 0][ arr_ix] != DB_MIN_VAL) // skip processing of NULL points
      {
        x = x_start + col * res;         

        if (!csv2db)
        {
          if (SQLrowcnt == 0)
          {
            db_set_string(sql,"INSERT INTO ");
            db_append_string(sql, tbl_name);
            db_append_string(sql, " VALUES (");
          }
          else db_append_string(sql, " ,(");
          SQLrowcnt++;
        }
        else db_set_string(sql, "");

        sprintf(buf, "%d,%d,%d",x,y,res);
        db_append_string(sql, buf);

        for (cell = 1; cell <= ncells; cell++) 
        {
          int arr_ix2 = arr_index[ cell-1][ arr_ix];
          sprintf(buf, ",'%s',%d,%.2f,'%s'", cell_name[ arr_ix2], antenna_id[ arr_ix2], arr_power[ cell-1][ arr_ix], model_name[ arr_ix2]);
          db_append_string(sql, buf);
        }

        sprintf(buf, ",%.2f",arr_EcNo[ arr_ix]);
        db_append_string(sql, buf);

        if (!csv2db) db_append_string(sql, ")");
        else fprintf(csv_fp, "%s\n", db_get_string(sql));
 
        G_debug(3, "%s", db_get_string(sql));
        //G_message(_("sql string: %s"), sql[0]);
      }
      else G_debug(3, "%s", "Skipping null value");

      if (!csv2db && SQLrowcnt > 0)
      {
        if (SQLrowcnt == SQLPacketSize || col+1 == ncols)
        {
          if (db_execute_immediate(driver,sql) != DB_OK) G_fatal_error(_("Failed writting data in the database! "));
          SQLrowcnt = 0;
        }
      }

    } //cols loop
  } //rows loop   

  if (!csv2db)
  {
    // COMMIT
//    db_set_string(sql,"COMMIT");
//    if (db_execute_immediate(driver,sql) != DB_OK)
    if (db_commit_transaction(driver) != DB_OK)
      G_fatal_error(_("Failed commiting transaction! "));
  }
  else
  {
    if ( !csv_only)
    {
      fflush(csv_fp);
      rewind(csv_fp);

      G_message(_("Converting intermediate temporary CSV file to data table..."));
      G_message(_("... please wait ..."));
      db_init_string(sql);
      if (strcmp("mysql", drv_name) == 0)
      {
        // LOAD DATA LOCAL INFILE ‘datoteka.csv' INTO TABLE <tabela> FIELDS TERMINATED BY ';' ENCLOSED BY "'";
        // LOAD DATA LOCAL INFILE ‘datoteka.csv' INTO TABLE <tabela> FIELDS TERMINATED BY ',' ENCLOSED BY '"';
        db_set_string(sql,"LOAD DATA LOCAL INFILE '");
        db_append_string(sql, csv_filename);
        db_append_string(sql,"' INTO TABLE ");
        db_append_string(sql, tbl_name);
        db_append_string(sql, " FIELDS TERMINATED BY ',' ENCLOSED BY \"'\"");
      }
      else if (strcmp("pg", drv_name) == 0)
      {
        // COPY <tabela> FROM ‘datoteka.csv' WITH DELIMITER ';' CSV QUOTE '''';
        // COPY <tabela> FROM ‘datoteka.csv' WITH CSV;
        db_set_string(sql,"COPY ");
        db_append_string(sql, tbl_name);
        db_append_string(sql," FROM '");
        db_append_string(sql, csv_filename);
        db_append_string(sql,"' CSV QUOTE ''''");
      }
      else G_fatal_error(_("dbperf=99 not supported for the database driver '%s' ! "), drv_name);
//      G_message(_("SQL: '%s'"), db_get_string(sql));
      if (db_execute_immediate(driver,sql) != DB_OK)
      G_fatal_error(_("Converting CSV to data table failed! "));
    }

    fclose(csv_fp);
  }
  //***close database and shut down database driver***
  if ( !csv_only)
  {
    db_close_database(driver);
    db_shutdown_driver(driver);
  }

  return DB_OK;
}

//------------------------------------------------------------------------------

  // check if a string is included in string list (delimited by commas and/or spaces)
  int str_in_strlist(char *strlist, char *str)
  {
    char *str1, *str2;

    str1 = strstr( strlist, str);
    str2 = str1 + strlen( str); // points past the found substring
    if ( !str1) return 0;
    if ( ( str1 == strlist                    || *(str1-1) == ' ' || *(str1-1) == ',') &&
         ( str2 == strlist + strlen( strlist) || *str2     == ' ' || *str2     == ','))
      return 1;
    return 0;
  }

//------------------------------------------------------------------------------

 int PdBm2LteThroughput( int, int, float*, float *, float *, char, char, double *, int *, int *, char *);

/*
 * main function
 */
int main(int argc, char *argv[])
{
       
  struct Cell_head window; /* database window */
  struct Cell_head cellhd; /* it stores region information, and header information of rasters */     
  struct History history; /* holds meta-data (title, comments,..) */
  struct GModule *module; /* GRASS module for parsing arguments */
  struct Option *output,*table_name,*file_name, *driver_name,*database, *dbperf, *cell_number, *rx_threshold; /* options*/
  struct Option *opt_generate, *opt_channel_type, *opt_bandwidth; /* more oprtions */
  struct Flag *flag_q; /* flags */

  char *outraster_name; /* output raster name */
  char *tbl_name; /* dbf table name */
  char *drv_name;
  char *db_name;
  char *in_file; /* input cell data file name */
  int db_perf;
  int verbose;
  int ovr;
  void *inrast; /* input buffer */
  unsigned char *outrast; /* output buffer */
  int infd,outfd; /* input and output file descriptor */

  int map_number; /* number of input raster maps*/
  int nrows, ncols;
  int row, col, map, n, i;
//  int x_start,y_start,res_x,res_y;
  int x_start,y_start,res_x;
  int cell_num; /* number of cells in dbf table*/

  double rx_thresh;

#define CELL_INPUT_PARAMS_MAX 20
  char *cell_input_param[CELL_INPUT_PARAMS_MAX + 1];
  int cell_input_params_num;
  int model_name_len;
  char *sep;
  char buffer[250];
  FILE *cell_input;

  FCELL f_in, f_out;
  int temp_index;
  float temp_float;

  //-----------------------------------------------------------------------------------------

  /* initialize GIS environment */
  G_gisinit(argv[0]); /* reads grass env, stores program name to G_program_name() */

  /* initialize module */
  module = G_define_module();
  G_add_keyword(_("raster"));
  G_add_keyword(_("MaxPower"));
//  module->description = _("RaPlaT - MaxPower module");
  char mod_desc[100];
  *mod_desc = 0;
  strncat( mod_desc, _("RaPlaT - MaxPower module"), 80);
  strcat( mod_desc, " (");
  strncat( mod_desc, VERSION, 15);
  strcat( mod_desc, ")");
  module->description = mod_desc;

  file_name = G_define_option();
  file_name->key = "cell_input";
  file_name->type = TYPE_STRING;
  file_name->required = YES;
  file_name->gisprompt = "old_file,file,input";
  file_name->description = _("Cells data table");

  opt_generate = G_define_option();
  opt_generate->required = NO;
  opt_generate->key = "generate";
  opt_generate->type = TYPE_STRING;
  opt_generate->description = _("Selection of the generated output contents");
  opt_generate->options = GENERATE_BASIC","GENERATE_LTE;
  opt_generate->answer = GENERATE_default;

  rx_threshold = G_define_option();
  rx_threshold->key = "rx_threshold";
  rx_threshold->type = TYPE_DOUBLE;
  rx_threshold->required = NO;
  rx_threshold->description = _("Minimum received power [dBm] for radio signal coverage");
//  rx_threshold->answer = "-999"; // Should be equal to DB_MIN_VAL !
  char DB_MIN_VAL_str[100];
  snprintf( DB_MIN_VAL_str, sizeof( DB_MIN_VAL_str) - 1, "%i", DB_MIN_VAL);
  rx_threshold->answer = DB_MIN_VAL_str;

  opt_channel_type = G_define_option();
  opt_channel_type->required = NO;
  opt_channel_type->key = "chan_type";
  opt_channel_type->type = TYPE_STRING;
  opt_channel_type->description = _("Channel type - Gaussian or Rayleigh (currently only Gaussian)");
  opt_channel_type->options = "gaussian";
  opt_channel_type->answer = "gaussian";

  opt_bandwidth = G_define_option();
  opt_bandwidth->required = NO;
  opt_bandwidth->key = "bandwidth";
  opt_bandwidth->type = TYPE_DOUBLE;
  opt_bandwidth->description = _("Bandwidth [MHz] (required for LTE computations)");
  opt_bandwidth->answer = "5";

  output = G_define_standard_option(G_OPT_R_OUTPUT);

  table_name = G_define_option();
  table_name->required = NO;
  table_name->key = "table";
  table_name->type = TYPE_STRING;
  table_name->description = _("Table name");

  driver_name = G_define_option();
  driver_name->required = NO;
  driver_name->key = "driver";
  driver_name->type = TYPE_STRING;
  driver_name->description = _("Driver name");
  char driver_name_options[100];
  strncpy( driver_name_options, db_list_drivers(), 90);
  strcat(driver_name_options, ",none,csv");
  driver_name->options = driver_name_options;
  driver_name->answer = "none";

  database = G_define_option();
  database->required = NO;
  database->key = "database";
  database->type = TYPE_STRING;
  database->description = _("Database name");
  database->answer = "$GISDBASE/$LOCATION_NAME/$MAPSET/dbf/";

  cell_number = G_define_option();
  cell_number->key = "cell_num";
  cell_number->type = TYPE_INTEGER;
  cell_number->required = NO;
  cell_number->description = _("Number of successive path loss values to be written in the table");
  cell_number->answer = "5";

  dbperf = G_define_option();
  dbperf->required = NO;
  dbperf->key = "dbperf";
  dbperf->type = TYPE_INTEGER;
  dbperf->description = _("Database insert performance (rows/INSERT; 99: special fast mode via CSV)");
  dbperf->options = "1-99";
  dbperf->answer = "20";


  /* define flags */
  flag_q = G_define_flag();
  flag_q->key = 'q';
  flag_q->description = _("Quiet");   


  /* options and flags parser */
  if (G_parser(argc, argv)) exit(EXIT_FAILURE);

  /* stores options and flags to variables */    
  outraster_name = output->answer;
  tbl_name = table_name->answer;
  drv_name = driver_name->answer;
  db_name = database->answer;
  verbose = !flag_q->answer;
  in_file = file_name->answer;
  sscanf( dbperf->answer, "%i", &db_perf);
  sscanf( cell_number->answer, "%d", &cell_num);
  sscanf( rx_threshold->answer, "%lf", &rx_thresh);

  char *generate = opt_generate->answer;
  char *channel_type = opt_channel_type->answer;
  double bandwidth;
  sscanf( opt_bandwidth->answer, "%lf", &bandwidth);

  ovr = 0;
  char *envGOverwriteStr = getenv( "GRASS_OVERWRITE");
  if  (envGOverwriteStr)
    if ( strcmp( envGOverwriteStr, "1") == 0) ovr = 1;


  //*****extract data from measurement file*****
  if( (cell_input = fopen(in_file,"r")) == NULL )
    G_fatal_error(_("Unable to open file <%s>"), in_file);
 

  // find the number of input cells (=lines of text) in the 'cell_input' stream
  map_number=0;
  while (fgets(buffer, 250, cell_input))
    map_number++;
  if (ferror( cell_input))
    G_fatal_error(_("Error reading cell_input file"));   
  if (map_number == 0)
    G_fatal_error(_("Empty cell_input file"));   
  G_message(_("Processing %d cells..."), map_number);
  rewind (cell_input);


  //***** create/overwrite data table (and tweak cell_num)
  if ( cell_num > map_number) cell_num = map_number;
  if ( strcmp( "none", drv_name) != 0 && strcmp( "csv", drv_name) != 0)
  { 
    if ( create_table( drv_name, db_name, tbl_name, cell_num, ovr) != DB_OK)
      G_fatal_error(_("Error creating/overwriting data table! ")); 
  }
  else if ( strcmp( "none", drv_name) == 0) cell_num = 1; //only highest signal for maxpower output raster

  char *name[map_number];  /* input raster file name */
  char *cell_name[map_number];  /* input cell name */
  int antenna_id[map_number];  /* antenna ID number */
  char *model_name[map_number];
  const char *mapset[map_number];  /* mapset name */
  double Pt[map_number];


  // *****store cell data****** 
  for (n = 0; n<map_number; n++)
  {
    if ( strcmp( fgets (buffer, 250, cell_input), buffer))
      G_fatal_error(_("Error reading cell_input file")); 
    else if ( buffer[ strlen(buffer)-1] != '\n')
      G_fatal_error(_("Error in cell_input file - line too long (\\n not found)"));
    buffer[ strlen(buffer)-1] = ';'; // replace \n (after the last parameter) with ;

    //parsing the output string
    cell_input_params_num = 0;
    cell_input_param[0] = buffer;
    for (i = 0; i < CELL_INPUT_PARAMS_MAX; i++)
    {
      sep = strchr( cell_input_param[i], ';');
      if ( !sep) break;
      *sep = '\0';
      cell_input_params_num++;
      cell_input_param[i + 1] = sep + 1;
    }

    if ( cell_input_params_num < 5)
      G_fatal_error(_("Number of parameters in a line of the cell_input file is too small (should be >=5)"));

    cell_name[n] = malloc( strlen( cell_input_param[0]) + 1); 
    strcpy( cell_name[n], cell_input_param[0]);  //get cell name from file
     
    antenna_id[n] = atoi( cell_input_param[1]);  //get antenna ID from file
     
    name[n] = malloc( strlen( cell_input_param[2]) + 1);
    strcpy( name[n], cell_input_param[2]);  //get input raster from file

    Pt[n] = atof( cell_input_param[3]);  //get cell transmit power from file
  

    model_name_len = 0;
    for( i = 4; i < cell_input_params_num; i++)
      model_name_len += strlen( cell_input_param[i]) + 1; // 1 for '_' (but no final '_')

    model_name[n] = malloc( model_name_len);
    strcpy( model_name[n], cell_input_param[4]);  //get the name of the channel model

    for (i = 5; i < cell_input_params_num; i++)
    {
      strcat( model_name[n], "_");
      strcat( model_name[n], cell_input_param[i]);
    }


    /* G_find_cell2 returns NULL if the map was not found in any mapset, mapset name otherwise */
    mapset[n] = G_find_raster(name[n], "");
    if (mapset == NULL) G_fatal_error(_("Raster map <%s> not found"), name[n]);
  }


  G_get_set_window(&window);

  if (G_legal_filename(outraster_name) < 0)
    G_fatal_error(_("<%s> is an illegal file name"), outraster_name);
   
  /* Allocate output buffer, use input map data_type */
  nrows = Rast_window_rows();  
  ncols = Rast_window_cols();
  outrast = Rast_allocate_buf(FCELL_TYPE);

  /* controlling, if we can write the raster */
  if ((outfd = Rast_open_new(outraster_name, FCELL_TYPE)) < 0)
    G_fatal_error(_("Unable to create raster map <%s>"), outraster_name);

  int num_points = nrows*ncols;

  /*POWER ARRAY*/
  float *tmp_arrpower = (float *)G_calloc( num_points * cell_num, sizeof(float));
  memset ( tmp_arrpower, 0, cell_num * num_points * sizeof(float));
  float **arr_power = (float **)G_calloc( cell_num, sizeof(float *));
  for ( i = 0 ; i < cell_num ; i++) arr_power [i] = tmp_arrpower + i * num_points;
   
  /*INDEX ARRAY*/
  int *tmp_arrindex = (int *)G_calloc( num_points * cell_num, sizeof(int));
  memset ( tmp_arrindex, 0, cell_num * num_points * sizeof(int));
  int **arr_index = (int **)G_calloc( cell_num, sizeof(int *));
  for ( i = 0; i < cell_num;i++) arr_index [i] = tmp_arrindex + i * num_points;
     
  /*ECNO ARRAY*/
  float *arr_EcNo = (float *)G_calloc( num_points, sizeof(float));
  memset ( arr_EcNo, 0, num_points * sizeof(float));
    
  /*SUM POWER ARRAY*/
  float *arr_sumpower = (float *)G_calloc( num_points, sizeof(double));
  memset ( arr_sumpower, 0, num_points * sizeof(double));

  /* Write rasters to array and sort power values and indexes - for each point */
  int count_mem=0;
  //G_message(_("\n...check_progress..., numn points = %d"),num_points);


  G_message(_("Sorting receive power values")); 

  long arr_ix;
  for (map = 0; map<map_number; map++)
  {
            
    if (verbose) G_percent(map+1, map_number, 2);

    /* G_open_cell_old - returns file descriptor (>0) */
    if ((infd = Rast_open_old(name[map], mapset[map])) < 0)
      G_fatal_error(_("Unable to open raster map <%s>"), name[map]);

    /* open input raster */
    Rast_get_cellhd(name[map], mapset[map], &cellhd);
     
    /* Allocate input buffer */
    inrast = Rast_allocate_buf(FCELL_TYPE);

    /* for each row */
    for (row = 0; row < nrows; row++)
    {
      Rast_get_row(infd, inrast, row, FCELL_TYPE);

      /* process the data */
      for (col = 0; col < ncols; col++)
      { 
        arr_ix = row * ncols + col;
        
        if ( map == 0) arr_sumpower[ arr_ix]=0;

        if ( Rast_is_f_null_value( &((FCELL *) inrast)[col])) f_in = DB_MIN_VAL;
        else f_in = (Pt[map] - ((FCELL *) inrast)[col]); //calculate receive power in dBm (Pr)

        //if (f_in > 0) f_in=DB_MIN_VAL; //in case of void (null) inrast value - to avoid case Pr=Pt

        if ( f_in > DB_MIN_VAL)
          arr_sumpower[ arr_ix] = arr_sumpower[ arr_ix] + pow( 10.0, f_in / 10.0); //sum power in mW 

        if ( map < cell_num)
        {  
          arr_power[ map][ arr_ix] = (float)f_in;
          arr_index[ map][ arr_ix] = map;
          count_mem = map;
        }
        else
        {
          if ((float)f_in < arr_power[ cell_num-1][ arr_ix]) continue;
          else
          {
            arr_power[ cell_num-1][ arr_ix] = (float)f_in;
            arr_index[ cell_num-1][ arr_ix] = map;
          }
        }

        //sort receive power values and cell indexes
        if(map) //do not sort for the first cell
        {
          for ( i = count_mem ; i > 0; i--)
          {
            if( arr_power[ i][ arr_ix] > arr_power[ i-1][ arr_ix])
            {
              temp_float = arr_power[i][ arr_ix];
              arr_power[   i][ arr_ix] = arr_power[ i-1][ arr_ix];
              arr_power[ i-1][ arr_ix] = temp_float;

              temp_index = arr_index[ i][ arr_ix];
              arr_index[   i][ arr_ix] = arr_index[ i-1][ arr_ix];
              arr_index[ i-1][ arr_ix] = temp_index;
            } 
            else break;    
          }
        }
      
      } //cols loop
    } //rows loop

    Rast_close(infd);
    G_free(inrast);      

  } //map (raster) loop

  G_message(_("Finished sorting receive power values")); 


  // At this point we have:
  //  - n strongest received signals at each raster point (sorted by value), in dBm
  //  - sum of all recevied signals, in mW

  // convert sumpower [mW] -> [dBm]
  for (row = 0; row < nrows; row++)
  {
    for (col = 0; col < ncols; col++)
    {
      arr_ix = row * ncols + col;
      if ( arr_sumpower[ arr_ix] == 0) arr_sumpower[ arr_ix] = DB_MIN_VAL;
      else
      {
        arr_sumpower[ arr_ix] = 10.0 * log10( arr_sumpower[ arr_ix]);
        if ( arr_sumpower[ arr_ix] < DB_MIN_VAL) arr_sumpower[ arr_ix] = DB_MIN_VAL;
      }
    }
  }

  // Prepare what is required for non-default output contents (default is GENERATE_default (rss-max))
  float *out_raster = NULL;
  int *out_raster_int = NULL;
  float *arr_out = NULL;

  int coverage_f = 0;

  if      ( strstr( "rss-max", generate))    out_raster = arr_power[0];
  else if ( strstr( "coverage", generate)) { out_raster = arr_power[0]; coverage_f = !0; }
  else if ( strstr( "rss-sum", generate))    out_raster = arr_sumpower;
  else if ( strstr( "rss-maxix", generate))  out_raster_int = arr_index[0];

  else if ( str_in_strlist( GENERATE_LTE, generate))
  {
    // PdBm2LteThroughput will be used to compute the output raster 
    float *arr_out = (float *)G_calloc( num_points, sizeof(double));
    out_raster = arr_out;

    char ChanType;
    if ( strcmp( channel_type, "gaussian") == 0) ChanType = 'g';
    else G_fatal_error(_("Wrong ChanType (internal error)"));

    char OutputFlag;
    if      ( strcmp( generate, "lte-rsrp") == 0)        OutputFlag = 'p';
    else if ( strcmp( generate, "lte-rssi") == 0)        OutputFlag = 'r';
    else if ( strcmp( generate, "lte-rsrq") == 0)        OutputFlag = 'q';
    else if ( strcmp( generate, "lte-cinr") == 0)        OutputFlag = 'c';
    else if ( strcmp( generate, "lte-maxspecteff") == 0) OutputFlag = 's';
    else if ( strcmp( generate, "lte-maxthrput") == 0)   OutputFlag = 't';
    else if ( strcmp( generate, "lte-interfere") == 0)    OutputFlag = 'i';
    else G_fatal_error(_("Wrong OutpuFlag (internal error)"));

    int nPDCCH = 2;   // number of physical downlink control channels
    char cpf = 'n';   // can be 'n' (normal) or 'e' (extended)
    int nAntenna = 1; // Number of transmit antennas 

    PdBm2LteThroughput( (int) nrows, (int) ncols, (float *) arr_power[0], (float *) arr_out, (float *) arr_sumpower,
                        (char) ChanType, (char) OutputFlag,
                        (double *) &bandwidth, (int *) &nPDCCH, (int *) &nAntenna, (char *) &cpf);
  }

  else G_fatal_error(_("'Generate' type not supported <%s> (internal error)"), generate);  


  //****write output raster and calculate Ec/No****
  for (row = 0; row < nrows; row++)
  {
    for (col = 0; col < ncols; col++)
    {
      arr_ix = row * ncols + col;

      // compute EcNo
      arr_EcNo[ arr_ix] = ( arr_power[ 0][ arr_ix] - arr_sumpower[ arr_ix]);

      // to the output raster map
      if ( out_raster) f_out = (CELL) out_raster [ arr_ix];
      else             f_out = (CELL) out_raster_int [ arr_ix];
      // change value DB_MIN_VAL or lower to 'null' (undefined, transparent)
      if ( f_out <= DB_MIN_VAL) Rast_set_f_null_value( &f_out, 1);

      // treshold processing - no output for points with the strongest signal below the treshold
      if ( arr_power[ 0][arr_ix] <= rx_thresh) Rast_set_f_null_value( &f_out, 1);
      else if ( coverage_f) f_out = 1.0;

      ((FCELL *) outrast)[ col] = f_out;

      //G_message(_("sum_power_dBm = %f"), arr_sumpower[ arr_ix]); 
      //G_message(_("EcN0 = %f"), arr_EcNo[ arr_ix]);
    }
         
    Rast_put_row( outfd, outrast, FCELL_TYPE);
  }
  Rast_close( outfd);
  G_free( outrast);


  //****write values in data table****

  if (strcmp("none",drv_name))
  {  
    x_start = (int)round( window.west + window.ew_res / 2.0);
    y_start = (int)round( window.north - window.ns_res / 2.0);
    res_x = (int)round( window.ew_res);
//    res_y = (int)round( window.ns_res);

    if (strcmp("csv",drv_name))
    {
      G_message(_("Writing MaxPower data in table '%s'..."), tbl_name);

      if ( fill_database(drv_name, db_name, tbl_name, db_perf, ncols, nrows, x_start, y_start, res_x, arr_power, arr_index,
                         cell_name, antenna_id, model_name, cell_num, arr_EcNo, verbose, ovr) != DB_OK)
        G_fatal_error(_("Error writing data in database! ")); 

      G_message(_("Finished writing MaxPower data in table '%s'..."), tbl_name);    
    }
    else
    {
      G_message(_("Writing MaxPower data in csv file '%s'..."), tbl_name);

      if ( fill_database(drv_name, db_name, tbl_name, db_perf, ncols, nrows, x_start, y_start, res_x, arr_power, arr_index,
                         cell_name, antenna_id, model_name, cell_num, arr_EcNo, verbose, ovr) != DB_OK)
        G_fatal_error(_("Error writing data in csv file! ")); 

      G_message(_("Finished writing MaxPower data in csv file '%s'..."), tbl_name);    
    }
  }    


  /* memory cleanup */
  G_free( arr_power);
  G_free( tmp_arrpower);
  G_free( arr_index);
  G_free( tmp_arrindex);
  G_free( arr_EcNo);
  G_free( arr_sumpower);

  if ( arr_out) G_free( arr_out);

  for ( n = 0; n < map_number; n++)
  {
    free( cell_name[n]);
    free( name[n]);
    free( model_name[n]);
  }


  /* add command line incantation to history file */
  Rast_short_history( outraster_name, "raster", &history);
  Rast_command_history( &history);
  Rast_write_history( outraster_name, &history);     

  exit( EXIT_SUCCESS);
}

