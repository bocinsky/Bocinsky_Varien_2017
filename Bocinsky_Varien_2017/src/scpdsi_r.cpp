#include <Rcpp.h>
using namespace Rcpp;

#include <stdlib.h>
using namespace std;

#define _USE_MATH_DEFINES
#include <math.h>

#include <cstring>
#include <stdio.h>
#include <ctype.h>
// #include <direct.h>  // unccoment for windows machines.


//=============================================================================
//pdsi.cpp              University of Nebraska - Lincoln            Jul 15 2003
//=============================================================================
//
// Current Version:  2.0
//
//Calculates the Palmer Drought Severity Index for a given station.  
//
//Recent modifications:
//  - Implemented user-defined calibration interval (Jun 5 2006)--goddard
//  - A more intuitive format for the parameter file is accepted:
//       AWC Latitude
//    The '-n' flag was added to facilitate the old format:
//       AWC TLA
//    (Jul 15 2003)
//  - Calculates the Weekly CMI (Jul 1 2003)
//  - Works in with southern latitude (Apr 1 2003)
//  - Works with the metric flag (Mar 12 2003)
//
//Self-calibrating, multiple time scale version.  Based upon weekly_pdsi.cpp
//which is a self-calibrating, weekly PDSI program developed by Nathan Wells.
//This version is capable of calculating a self-calibrating weekly PDSI, 
//a self-calibrating monthly PDSI, and the original monthly PDSI.
//
//Most recently translated to C++ from FORTRAN by Rob Reed and Nate Wells,
//University of Nebraska - Lincoln, advised by Dr. Steve Goddard - July 2001. 
//
//Methodology based on Research Paper No. 45; Meteorological Drought; by
//Wayne C. Palmer for the U.S. Weather Bureau, February 1965.
//
//Based mostly on the source code of pdsi.f, a FORTRAN program that calculates 
//monthly PDSI values.  The source code came from NCDC, originally written by 
//Tom Karl, and revised by Roger Winchell, Alan McNab, and Richard Heim.  
//
//Slight modifications in the algorithm were made based on the source of 
//palmcode.f, a FORTRAN program that calculates weekly PDSI values, received 
//from Tom Heddinghaus at NCEP, who is also the original author of that 
//particular code.
//
//Additional modifications were made to adapt the algorithm to a weekly time
//scale based on recalculations of several constants as described in Palmer's
//paper.  
//
//Additional modifications were made to attempt to make this program 
//completely independent of emperically derived formulas.  This will allow the
//program to perform accurately in any enviroment.  Most changes came in the 
//addition of the Calibrate() function.  These changes were made in order to 
//make comparisions between stations more accurate.  --August 2001
//
//The incorporation of a self-calibrating monthly and the oringal monthly PDSI
// -- May 2002
//
//Changed the calibration method to calibrate the climatic characteristic based
//on the quartiles instead of the max and min.  -- Jun 2002.
//-----------------------------------------------------------------------------
//
// 4 input files for weekly PDSI calculations:
//
//weekly_T and weekly_P:
//  Weekly temperature and precipitation files; each line starts with the year
//  and is followed by 52 data entries.  It is important to note that only 52
//  weeks are on each line, even though 52 weeks is only 364 days.  This data 
//  must be gathered in such a way that the kth week of the year always 
//  represents the same calendar interval.  For example, the first week should
//  always represent Jan 1 through Jan 7.  
//
//wk_T_normal:
//  The average weekly temperature for each week of the year.  One line, 52
//  data entries.
//
//parameter:
//  contains two numbers specific to each station: the underlying soil water
//  capacity - Su (in inches), and the negative of the tangent of the latitude
//  of the station - TLA.  
//
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
// 4 input files for monthly PDSI calculations:
//
//monthly_T and monthly_P:
//  Monthly temperature and precipitation files; each line starts with the year
//  and is followed by 12 data entries.  
//
//mon_T_normal:
//  The average monthly temperature for each week of the year.  One line, 12 
//  data entries.
//
//parameter:
//  same as above.
//
//-----------------------------------------------------------------------------
//Extra notes on input files:
//  The format of these files matches the format of the original FORTRAN 
//  program.  There is no precise need for this format to be used.
//  
//  The program is able to calculate the weekly PDSI, the original monthly PDSI
//  and the monthly self-calibrating PDSI if all the data is in the same 
//  directory.
//
//  It is possible to use the filename T_normal in place of either mon_T_normal
//  or wk_T_normal.  For example, the program will try to open mon_T_normal 
//  first, and then T_normal.  If T_normal is opened, it will check to make 
//  sure it is the right format by counting the number of entries in the file.
//  This was done to allow the program to work on the exact same input data as
//  the original FORTRAN program, allowing comparisons.
//-----------------------------------------------------------------------------
//
//------ Output Files:
//
//There are two formats of output, table and column, which are selected by 
//command line options.  The table output files end with .tbl and the column
//output files end with .clm.  The table files list a whole year's worth of 
//resultsing values on each line while the column files list the year and week 
//followed by the resulting value on each line.
//
//PDSI.tbl and/or PDSI.clm:
//  The Palmer Drought Severity Index values
//
//PHDI.tbl and/or PHDI.clm:
//  The Palmer Hydrological Drought Index values.
//
//WPLM.tbl and/or WPLM.clm:
//  The "Weighted" Palmer Index.  An average of either X1 and X3, or X1 and X2
//  weighted by the probability of the current spell ending.  For more 
//  information, see how the WPLM is calculated in the pdsi::write() function.
//
//ZIND.tbl and/or ZIND.clm:
//  The Z Index, or the moisture anomaly
//
//------ Other possible output files, depending on certain flags:
//
//WB.tbl
//  The water ballance coefficients (Alpha, Beta, Gamma, & Delta) for each
//  week/month.
//
//bigTable.tbl
//  Z, % Prob(end), X1, X2, and X3 for every week/month.
//
//potentials
//  P, PE, PR, PRO, PL, and P - PE for every week/month.
//
//dvalue
//  The moisture departure, d, for every week/month.
//
//=============================================================================
//           end of introductory comments
//=============================================================================

// This defines the type number as a double.  This is used to easily change
// the PDSI's variable types.
typedef double number;
typedef int flag;
#define min(a,b) ((a) < (b) ? (a) : (b));
#define MISSING -99.00

//-----------------------------------------------------------------------------
//**********   START OF STRUCTURE DEFINITION FOR STRUCT:  node        *********
//-----------------------------------------------------------------------------
// The node struct is used specifically in the linked list class and is not    
// relevant to the actual PDSI.                                                
//-----------------------------------------------------------------------------
struct node {            // A structure for a node
public:
  number key;            // Where the data is stored
  struct node *next;     // Where the next node is
  struct node *previous; // Where the previous node is
};

//-----------------------------------------------------------------------------
//**********   CLOSE OF STRUCTURE DEFINITION FOR STRUCT:  node        *********
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//**********   START OF CLASS DEFINITIONS FOR THE CLASS:  llist       *********
//-----------------------------------------------------------------------------
// The llist class is a dynamic storage class.  It is used in the PDSI to 
// eliminate problems with filling static arrays.  
//-----------------------------------------------------------------------------
class llist {            // A linked list class
private:
  node *head;            // A pointer to the head of the linked-list
  int size;
  number kthLargest(int k);  //returns the kth largest number in the list
  number long_select(int k); //returns the kth largest number in the list
  //without using an array or sorting.  Used when
  //there is not enough memory to place the entire
  //list in a new array.
  number percentile(double percentage); //returns the specified percentile
  
public:
  llist();               // The constructor
  ~llist();              // The destructor
  // The insert function takes an argument of type number and places it on
  // the head of the linked list.
  void insert(number x); 
  // The remove functions remove from either the head (head_remove) or the
  // tail (tail_remove) of the linked list.
  number head_remove();  // remove the first node and returns its key
  number tail_remove();  // remove the last node and returns its key
  // These are other useful functions used in dealing with linked lists
  int is_empty();// Returns 1 if the llist is empty 0 otherwise
  int get_size();
  number sumlist();  // Sums the items in list
  void sumlist(number &prev_sum, int sign);//sums items in list into prev_sum
  number maxlist();
  number minlist();
  number quartile(int q);        //returns the qth quartile
  number safe_percentile(double percentage); //safe version 
  
  number* returnArray();
  
  // The set_node function sets the key of the node pointed to by set.  It
  // checks the linked list to make sure set is in the node to prevent 
  // runtime errors from occuring.  It was written specifically for the PDSI
  // program and its backtracking function.
  node *set_node(node *set=NULL, number x=0);
  friend void copy(llist &L1,const llist &L2); // Copies L2 to L1
};
//-----------------------------------------------------------------------------
//**********   CLOSE OF CLASS DEFINITIONS FOR THE CLASS:  llist       *********
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//**********   START OF CLASS DEFINITIONS FOR THE CLASS:  pdsi        *********
//-----------------------------------------------------------------------------
// The pdsi class contains and deals with all necessary pdsi variables and
// calculations.  The paramater values and starting/ending years are set upon
// creation, but the actual pdsi values are not calculated until Calcpdsi is
// called.
//-----------------------------------------------------------------------------
class pdsi {
public:
  // The pdsi constructor takes in an array of arguments (argv[]) and an
  // integer number of arguments (argc).  These arguments should contain
  // the various flags desired to custimize the performance of the pdsi
  pdsi();
  ~pdsi();
  
  // The set_flags function takes in an array of flags/arguments (flags[])
  // and an integer number of flags/arguments (num_flags).  It then sets
  // the various pdsi options accordingly.  If incorrect or invalid flags
  // are input, the set_flags function terminates the program.
  void set_flags(int num_flags, char *flags[]);
  
  number getPDSI(int period, int year);
  number getZIND(int period, int year);
  number getPHDI(int period, int year);
  number getWPLM(int period, int year);
  
  number* getPDSIArray(int &size);
  number* getPDSIArray(int start_per, int start_yr, 
                       int end_per, int end_yr, int &size);
  number* getZINDArray(int &size);
  number* getZINDArray(int start_per, int start_yr,
                       int end_per, int end_yr, int &size);
  number* getWPLMArray(int &size);
  number* getWPLMArray(int start_per, int start_yr, 
                       int end_per, int end_yr, int &size);
  number* getPHDIArray(int &size);
  number* getPHDIArray(int start_per, int start_yr,
                       int end_per, int end_yr, int &size);
  
  number* getYearArray(int &size);
  number* getPerArray(int &size);
  
  NumericVector getPDSIArray_R();
  void set_flags_R();
  void MonthlyPDSI_R(NumericVector monthly_T,
                     NumericVector monthly_P,
                     NumericVector mon_T_normal,
                     double awc,
                     double lat,
                     int start_year,
                     int end_year);
  
private:
  //these variables keep track of what type of PDSI is being calculated.
  bool Weekly;
  bool Monthly;
  bool SCMonthly;
  
  // The variables for storing the starting year of calculation and the
  // total number of years to calculate
  int startyear;
  int endyear;
  int totalyears;
  
  int period_length;        //set to 1 for monthly, otherwise, legth of period
  int num_of_periods;       //number of periods of period_length in a year.
  
  // The variables used as flags to the pdsi class
  flag bug;
  flag output_mode;
  flag verbose;
  flag s_year;
  flag e_year;
  flag extra;
  flag metric;
  flag south;
  flag nadss;
  
  /* added on 9/21/05 to allow for user-define calibration start year (jdokulil) */
  flag setCalibrationStartYear;
  flag setCalibrationEndYear;
  int calibrationStartYear;
  int calibrationEndYear;
  /* end addition */
  
  /* SG: Steve Goddard modifications */  
  /* SG 6/5/06: add variables to allow user-defined calibration intervals */
  int currentCalibrationStartYear;
  int currentCalibrationEndYear;
  int nStartYearsToSkip;
  int nEndYearsToSkip;
  int nCalibrationYears;
  int nStartPeriodsToSkip;
  int nEndPeriodsToSkip;
  int nCalibrationPeriods;
  /* SG 6/5/06: End adding variables to allow user-defined calibration intervals */  
  
  
  // Various constants used in calculations
  number TLA; // The negative tangent of latitude is used in calculating PE
  number AWC; // The soils water capacity
  number I;   // Thornthwaites heat index
  number A;   // Thornthwaites exponent
  number tolerance; // The tolerance for various comparisons
  
  // The arrays used to read in the normal temp data, a years worth of 
  // actual temp data, and a years worth of precipitation data
  number TNorm[52];
  number T[52];
  number P[52];
  
  // These variables are used in calculation to store the current period's
  // potential and actual water balance variables as well as the soil
  // moisture levels
  number ET;            // Actual evapotranspiration
  number R;             // Actual soil recharge 
  number L;             // Actual loss
  number RO;            // Actual runoff
  number PE;            // Potential evapotranspiration
  number PR;            // Potential soil recharge
  number PL;            // Potential Loss
  number PRO;           // Potential runoff
  number Su;            // Underlying soil moisture
  number Ss;            // Surface soil moisture
  
  // These arrays are used to store the monthly or weekly sums of the 8 key 
  // water balance variables and the precipitation
  number ETSum[52];
  number RSum[52];
  number LSum[52];
  number ROSum[52];
  number PESum[52];
  number PRSum[52];
  number PLSum[52];
  number PROSum[52];
  number PSum[52];
  
  // These arrays store the monthly or weekly water balance coefficients 
  number Alpha[52];
  number Beta[52];
  number Gamma[52];
  number Delta[52];
  
  // The CAFEC percipitation
  number Phat;
  
  // These variables are used in calculating the z index
  number d;     // Departure from normal for a period
  number D[52]; // Sum of the absolute value of all d values by period
  number k[52]; // Palmer's k' constant by period
  number K;     // The final K value for a period
  number Z;     // The z index for a period (Z=K*d)
  
  // These variables are used in calculating the PDSI from the Z
  // index.  They determine how much of an effect the z value has on 
  // the PDSI based on the climate of the region.  
  // They are calculated using CalcDurFact()
  number drym;
  number dryb;
  number wetm;
  number wetb;
  
  //these two variables weight the climate characteristic in the 
  //calibration process
  number dry_ratio; 
  number wet_ratio;
  
  // The X variables are used in book keeping for the computation of
  // the pdsi
  number X1;    // Wet index for a month/week
  number X2;    // Dry index for a month/week
  number X3;    // Index for an established wet or dry spell
  number X;     // Current period's pdsi value before backtracking
  
  // These variables are used in calculating the probability of a wet
  // or dry spell ending
  number Prob;  // Prob=V/Q*100
  number V;     // Sumation of effective wetness/dryness
  number Q;     // Z needed for an end plus last period's V
  
  // These variables are statistical variables computed and output in 
  // verbose mode
  number DSSqr[52];
  number DEPSum[52];
  number DKSum;
  number SD;
  number SD2;
  
  // linked lists to store X values for backtracking when computing X
  llist Xlist;//final list of PDSI values
  llist altX1;//list of X1 values
  llist altX2;//list of X2 values
  
  // These linked lists store the Z, Prob, and 3 X values for
  // outputing the Z index, Hydro Palmer, and Weighted Palmer
  llist XL1;
  llist XL2;
  llist XL3;
  llist ProbL;
  llist ZIND;
  llist PeriodList;
  llist YearList;
  
  // Class Functions
  
  // These funcitons initialize some variables and linked lists
  // and check to make sure all the necessary data is there.
  int initialize();
  
  // These functions calculate the Potentials needed
  // Calculates Potential Evapotranspiration from Thornthwaite
  void CalcMonPE(int month, int year);
  void CalcPR();// Calculates Potential Recharge
  void CalcPL();// Calculates Potential Loss
  void CalcPRO();// Calculates Potential Runoff
  
  // This function calculates the actual values from the potentials
  void CalcActual(int per);// Calculates Actual values
  
  // and the X values.  Used for uncalibrated PDSI.
  void CalcZ();     // Calculates the Z-index
  void CalcX();     // Calculates the PDSI and X1, X2, and X3
  
  // This function backtracks through the X values 
  // and replaces them with the appropriate value of X1 or X2
  // when necessary
  void Backtrack(number X1, number X2);
  void ChooseX(number& newX, number& newX1, number& newX2,
               number& newX3, int bug);
  
  // These functions calculate Thornthwaite coefficients
  number CalcThornA(number I);// Calculates Thornthwaite Exponent
  
  // This function calculates the constants used in calculating the PDSI value
  // from the Z index.  These constants affect how much influence the Z index
  // has on the PDSI value.
  void CalcDurFact(number &slope, number &intercept, int sign);
  number get_Z_sum(int length, int sign);
  void LeastSquares(int *x, number *y, int n, int sign, number &slope, number &intercept); 
  
  number getValue(llist &List, int period, int year);
  
  number* getSubArray(llist &List, int start_per, int start_yr, 
                      int end_per, int end_yr, int &size);
  
  // These are simple functions to determine if the characters of a string 
  // form an integer or a floating point number.
  inline int is_int(char *string,int length);
  inline int is_flt(char *string,int length);
  
  
  number CalcMonThornI_R();
  void SumAll_R(); // Creates sums of Actual and Potential values
  void Calcd_R();
  void CalcWBCoef_R();
  void CalcOrigK_R();
  void CalcOneX_R(int period_number, int year);
  // void CalcZ_R();
  
  NumericVector monthly_temp;
  NumericVector monthly_precip;
  NumericVector mon_temp_normal;
  IntegerVector all_years;
  
  NumericVector potentials_P;
  NumericVector potentials_PE;
  NumericVector potentials_PR;
  NumericVector potentials_PRO;
  NumericVector potentials_PL;
  NumericVector potentials_PPE;
  
  NumericVector dvalue_d;
};
//-----------------------------------------------------------------------------
//**********   CLOSE OF CLASS DEFINITIONS FOR THE CLASS:  pdsi        *********
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//These three functions, partition(), select(), and exch() are used to select
//the kth largest number in an array.
//-----------------------------------------------------------------------------
void select(number a[], int l, int r, int k);
int partition(number a[], int left, int right);
void exch(number &x, number &y);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//**********                       R PROGRAM                       *********
//-----------------------------------------------------------------------------
// The main program takes in command line arguments and passes them to the 
// constructor of a pdsi variable tester.  It then calls Calcpdsi to calculate
// the pdsi values.  Finally it calls write to output these values to file.
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector scpdsi(NumericVector monthly_T,
                     NumericVector monthly_P,
                     NumericVector mon_T_normal,
                     double awc,
                     double lat,
                     int start_year,
                     int end_year) {
  pdsi PDSI;
  
  PDSI.set_flags_R(); // Sets the flags of PDSI
  
  PDSI.MonthlyPDSI_R(monthly_T,
                     monthly_P,
                     mon_T_normal,
                     awc,
                     lat,
                     start_year,
                     end_year);
  
  return PDSI.getPDSIArray_R();
}
//-----------------------------------------------------------------------------
//**********                      R PROGRAM END                       *********
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//**********   START OF FUNCTION DEFINITIONS FOR CLASS:  llist        *********
//-----------------------------------------------------------------------------
// The pdsi constructor sets all flags to default, scans the temperature file 
// to get the starting and ending years of the data, and reads in the values 
// from the parameter file.
//-----------------------------------------------------------------------------
pdsi::pdsi() {
  
  //set several parameters to their defaults
  period_length = 1;
  num_of_periods = 52;
  verbose=0;
  bug=0;
  output_mode=0;
  tolerance=0.00001;
  metric=0;
  nadss=0;
  setCalibrationStartYear=0;
  setCalibrationEndYear=0;
}
//-----------------------------------------------------------------------------
// The destructor deletes the temporary files used in storing various items
//-----------------------------------------------------------------------------
pdsi::~pdsi() {
  
  if(verbose > 0)
    printf("Calculations Complete\n");
}

//-----------------------------------------------------------------------------
// The initialize function clears all arrays and linked lists.
// it also checks to make sure all the necessary input files exist.
// returns 1 if they do, and -1 if they do not.
// returns 0 if there is too much missing data.  
//  -- there should be 25 years worth of data. (300 months or 1300 weeks)
//-----------------------------------------------------------------------------
int pdsi::initialize() {
  while(!Xlist.is_empty())
    Xlist.tail_remove();
  while(!altX1.is_empty())
    altX1.tail_remove();
  while(!altX2.is_empty())
    altX2.tail_remove();
  while(!XL1.is_empty())
    XL1.tail_remove();
  while(!XL2.is_empty())
    XL2.tail_remove();
  while(!XL3.is_empty())
    XL3.tail_remove();
  while(!ProbL.is_empty())
    ProbL.tail_remove();
  while(!ZIND.is_empty())
    ZIND.tail_remove();
  while(!PeriodList.is_empty())
    PeriodList.tail_remove();
  while(!YearList.is_empty())
    YearList.tail_remove();
  // monthly_temp.erase(0,monthly_temp.length()-1);
  // monthly_precip.erase(0,monthly_precip.length()-1);
  // mon_temp_normal.erase(0,mon_temp_normal.length()-1);
  // all_years.erase(0,all_years.length()-1);
  // potentials_P.erase(0,potentials_P.length()-1);
  // potentials_PE.erase(0,potentials_PE.length()-1);
  // potentials_PR.erase(0,potentials_PR.length()-1);
  // potentials_PRO.erase(0,potentials_PRO.length()-1);
  // potentials_PL.erase(0,potentials_PL.length()-1);
  // potentials_PPE.erase(0,potentials_PPE.length()-1);
  // dvalue_d.erase(0,dvalue_d.length()-1);
  
  return 1;
}

//-----------------------------------------------------------------------------
// The set_flags function takes an argument list (flags) and a number of 
// arguments (num_flags).  It then sets the various pdsi flags accordingly
//-----------------------------------------------------------------------------
void pdsi::set_flags_R() {
  int read_from =1;
  unsigned int n;
  
  // These flags are used for two argument flags
  flag week=0, year=0, both=0, out=0, e_flag=0, t_year=0;
  flag in_dir = -1, out_dir = -1;
  // Initializes the output years flags to 0
  s_year = 0;
  e_year = 0;
  extra = 0;
  
  initialize();
  
}// End set_flags

//-----------------------------------------------------------------------------
// MonthlyPDSI calls all functions necessary to compute the final PDSI.  
// This is done using the 4 functions SumAll, CalcWBCoef, Calcd, and CalcK.
// Produces a the original PDSI
//-----------------------------------------------------------------------------
void pdsi::MonthlyPDSI_R(NumericVector monthly_T,
                         NumericVector monthly_P,
                         NumericVector mon_T_normal,
                         double awc,
                         double lat,
                         int start_year,
                         int end_year) {
  
  monthly_temp = monthly_T;
  monthly_precip = monthly_P;
  mon_temp_normal = mon_T_normal;
  
  Monthly = true;
  Weekly = false;
  SCMonthly = false;
  
  period_length = 1;
  num_of_periods = 12;
  
  startyear = start_year;
  endyear = end_year;
  totalyears = endyear - startyear + 1;
  all_years = Rcpp::IntegerVector(totalyears);
  
  potentials_P = Rcpp::DoubleVector(totalyears * num_of_periods);
  potentials_PE = Rcpp::DoubleVector(totalyears * num_of_periods);
  potentials_PL = Rcpp::DoubleVector(totalyears * num_of_periods);
  potentials_PPE = Rcpp::DoubleVector(totalyears * num_of_periods);
  potentials_PR = Rcpp::DoubleVector(totalyears * num_of_periods);
  potentials_PRO = Rcpp::DoubleVector(totalyears * num_of_periods);
  
  dvalue_d = Rcpp::DoubleVector(totalyears * num_of_periods);
  
  int syear = startyear;
  
  for(int year = 1; year <= totalyears; year++) {
    all_years[year-1] = syear;
    // cout << syear;
    // cout << "\n";
    syear++;
  }
  
  /* SG 6/5/06: The original Monthly PDSI should not support a calibration interval, so clear the vars */
  currentCalibrationStartYear = startyear;
  currentCalibrationEndYear = endyear;
  nEndYearsToSkip = 0;
  nStartYearsToSkip = 0;
  nCalibrationYears = totalyears;
  nStartPeriodsToSkip = 0;
  nEndPeriodsToSkip = 0;
  nCalibrationPeriods = nCalibrationYears * num_of_periods;
  /* SG 6/5/06: end addition */
  
  // This block opens the parameter file and sets the initial Su and TLA values
  // must be called after the variable period_length is determined in the 
  // set_flags function
  
  
  AWC = awc;
  if(AWC <= 0) {
    printf("Invalid value for AWC: %f\n",Su);
    exit(0);
  }
  Ss = 1.0;   //assume the top soil can hold 1 inch
  if(AWC < Ss){
    //always assume the top layer of soil can 
    //hold at least the Ss value of 1 inch.
    AWC = Ss;
  }
  Su = AWC - Ss;
  if(Su < 0)
    Su = 0;
  
  TLA = -tan(M_PI*lat/180);
  if(lat >= 0){
    if(verbose>1)
      printf("TLA is positive, assuming location is in Southern Hemisphere. TLA: %f\n",TLA);
    south = 0;
  }
  else{
    south = 1;
    TLA = -TLA;
  }
  
  I = CalcMonThornI_R();
  
  A = CalcThornA(I);
  if(verbose>1) {
    printf ("AWC = % .9f  TLA = % .8f\n", AWC, TLA);
    printf ("HEAT INDEX, THORNTHWAITE A: %14.5f %11.5f\n", I, A);
  }
  
  // Output seen only in maximum verbose mode
  if(verbose>1)
    printf ("processing station 1\n");
  // SumAll is called to compute the sums for the 8 water balance variables
  SumAll_R();
  // This outputs those sums to the screen
  if(verbose>1) {
    printf ("STATION = %5d %18c SUMMATION OF MONTHLY VALUES OVER %4d YEARS\n", 0, ' ', totalyears);
    printf ("%36c CALIBRATION YEARS:\n", ' ');
    printf ("%4s %7s %8s %8s %8s %8s %8s %8s %8s %8s %10s", "PER", "P", "S", "PR", "PE", "PL", "ET", "R", "L", "RO", "DEP\n\n");
  }
  for (int i=0;i<num_of_periods;i++) {
    DEPSum[i] = ETSum[i] + RSum[i] - PESum[i] + ROSum[i];
    if(verbose>1) {
      printf ("%4d", (period_length*i)+1);
      printf ("%9.2f %8.2f %8.2f %8.2f %8.2f", PSum[i], PROSum[i], PRSum[i], PESum[i], PLSum[i]);
      printf ("%9.2f %8.2f %8.2f %8.2f %8.2f", ETSum[i], RSum[i], LSum[i], ROSum[i], DEPSum[i]);
      printf ("\n");
    }
    DSSqr[i] = 0;
  }
  
  // CalcWBCoef is then called to calculate alpha, beta, gamma, and delta
  CalcWBCoef_R();
  // Next Calcd is called to calculate the monthly departures from normal
  Calcd_R();
  // Finally CalcK is called to compute the K and Z values.  CalcX is called
  // within CalcK.
  CalcOrigK_R();
  
}


void pdsi::CalcMonPE(int month, int year) {
  number Phi[]={-.3865982,-.2316132,-.0378180,.1715539,.3458803,.4308320,.3916645,.2452467,.0535511,-.15583436,-.3340551,-.4310691}; 
  //these values of Phi[] come directly from the fortran program.
  int Days[]={31,28,31,30,31,30,31,31,30,31,30,31};
  number Dum, Dk;
  int offset;
  if(south)
    offset = 6;
  else
    offset = 0;
  
  if (T[month] <= 32)
    PE=0;
  else { 
    Dum=Phi[(month+offset)%12]*TLA;
    Dk = atan(sqrt(1-Dum*Dum)/Dum);
    if (Dk < 0)
      Dk += 3.141593;
    Dk = (Dk +.0157)/1.57;
    if (T[month]>=80)
      PE=(sin(T[month]/57.3-.166)-.76)*Dk;
    else {
      Dum=log(T[month]-32);
      PE=(exp(-3.863233+A*1.715598-A*log(I)+A*Dum))*Dk;
    }
  }
  
  // This calculation of leap year follows the FORTRAN program 
  // It does not take into account factors of 100 or 400
  /*
   if (year%4==0 && month==1)
   PE=PE*29;
   else
   PE=PE*Days[month];
   */
  //this calculation has been updated to accurately follow leap years
  if(month == 1){
    if(year%400 == 0)
      PE=PE*29;
    else if(year%4 == 0 && year%100 != 0)
      PE=PE*29;
    else
      PE=PE*28;
  }
  else
    PE=PE*Days[month];
  
}

//-----------------------------------------------------------------------------
// This function calculates the Thornthwaite heat index I.  This is done by 
// reading in the weekly normal temperature from file.  Any above freezing 
// temperatures are adjusted and then added to the index.  The equations have
// been modified to handle temperature in degrees Fahrenheit
//This function calculates the Thornthwaite heat index I for montly PDSI's.
//-----------------------------------------------------------------------------
number pdsi::CalcMonThornI_R() {
  number I=0;
  int i=0,j=0;
  float t[13];

  // The monthly temperatures are read in to a temparary array.
  // This was done because the fscanf function was unable to handle an array
  // of type double with the %f.  There might be something that could be used
  // in place of the %f to get a double to work.
  
  for(int i = 0; i < 12; i++){
    t[i] = mon_temp_normal[i];
  }
  
  //check to make sure file only had 12 entries

  // Then we move the temperatures to the TNorm array and calclulate I
  for(i=0;i<12;i++) {
    if(metric)
      TNorm[i] = t[i]*(9.0/5.0)+32;
    else
      TNorm[i]=t[i];
    // Prints the normal temperatures to the screen
    if(verbose>1)
      printf (" %.2f", TNorm[i]);
    // Adds the modified temp to heat if the temp is above freezing
    if (TNorm[i]>32)
      I=I+pow((TNorm[i]-32)/9,1.514);
  }
  // Prints a newline to the screen and closes the input file
  if(verbose>1)
    printf ("\n");

  return I;
}


//-----------------------------------------------------------------------------
// CalcThornA calculates the Thornthwaite exponent a based on the heat index I.
//-----------------------------------------------------------------------------
number pdsi::CalcThornA(number I) {
  number A;
  A=6.75*(pow(I,3))/10000000 - 7.71*(pow(I,2))/100000 + 0.0179*I + 0.49;
  return A;
}

//-----------------------------------------------------------------------------
// This function calculates d, the moisture departure, for each period of every
// year examined.
//
// d = P - Phat
//
// Phat = (alpha * PE) + (beta * PR) + (gamma * PRO) - (delta * PL)
//
// where PE, PR, PRO, and PL are read in from a file named "potentials"
// where each line has year, month/week, P, PE, PR, PRO, PL and P-PE
//
// this function also calculates the mean of the absolute values of d
// over all the years examined, called D, which is used to  
// calculate k and K
//
// Note:  Phat is P with a ^ on it in the mathematical equations explaining
// the Palmer Drought Index
//-----------------------------------------------------------------------------
void pdsi::Calcd_R() {
  // FILE *fin;        // File pointer for the temp. input file potentials
  // FILE *fout;       // File pointer for the temp. output file dvalue
  int per;           // The period in question
  int yr;           // The year in question
  number p;         // The precip for that period
  // float scn1, scn2, scn3, scn4, scn5, scn6; // Temp. variables for fscanf
  char letter = ' ';
  int i = 0;
  // These variables are used in calculating terminal outputs and are not 
  // important to the final PDSI
  number D_sum[52];
  number DSAct[52];
  number SPhat[52];
  
  for(i=0;i<52;i++){
    D_sum[i] = 0.0;
    DSAct[i] = 0.0;
    SPhat[i] = 0.0;
  }
  
  // The potentials file is opened for reading in the previously stored values
  // if((fin=fopen("potentials","r"))==NULL) {
  //   if(verbose>0)
  //     printf ("Error opening temp file with all the potential values.\n");
  //   exit(1);
  // }
  // 
  // if((fout=fopen("dvalue","w")) == NULL) { 
  //   if(verbose>0)
  //     printf ("Error opening temp file for d values.\n");
  //   exit(1);
  // }
  
  // This reads in the first line, which contains column headers.
  // while(letter != '\n')
  //   letter = fgetc(fin);
  
  // This reads in the values previously stored in "potentials"
  for(int i = 0; i < potentials_P.length(); i++) {
    // while(fscanf(fin,"%d %d %f %f %f %f %f %f", &yr, &per, &scn1, &scn2, &scn3, &scn4, &scn5, &scn6) != EOF) {
    per = i % 12;   //adjust the period # for use in arrays.
    
    // cout << per;
    // cout << "\n";
    
    p=potentials_P[i];
    PE=potentials_PE[i];
    PR=potentials_PR[i];
    PRO=potentials_PRO[i];
    PL=potentials_PL[i];
    //scn6 is P - PE, which can be ignored for calculations.
    
    if(p!=MISSING&&PE!=MISSING&&PR!=MISSING&&PRO!=MISSING&&PL!=MISSING){
      // Then the calculations for Phat and d are done
      Phat=(Alpha[per]*PE)+(Beta[per]*PR)+(Gamma[per]*PRO)-(Delta[per]*PL);
      d=p - Phat;
      
      dvalue_d[i] = d;
      // The values of d are output to a temp file for later use.
      // fprintf (fout, "%d %d %f\n", yr, (period_length*per)+1, d);
      
      /* SG 6/5/06: Need to only update statistical values when in 
       **            user defined calibration interval. When not used
       **            nCalibrationYears==totalyears; hence no change then
       */
      if (yr >= currentCalibrationStartYear && yr <= currentCalibrationEndYear) {
        // D_sum is the sum of the absolute values of d
        // and is used to find D
        // D_sum[per] += abs(d);
        if(d < 0.0)
          D_sum[per] += -(d);
        else
          D_sum[per] += d;
        
        
        // The statistical values are updated
        DSAct[per] += d;
        DSSqr[per] += d*d;
        SPhat[per] += Phat;
      }
    }
    else {
      d = MISSING;
      dvalue_d[i] = d;
      // fprintf (fout, "%d %d %f\n", yr, (period_length*per)+1, d);
    }
  }
  
  // If the user specifies, the various sums are output to the screen
  if(verbose>1) {
    printf ("\n%34c CHECK SUMS OF ESTIMATED VARIABLES\n\n", ' ');
    printf ("%4s %8s %8s %8s", "PER", "SCET", "SCR", "SCRO");
    printf ("%9s %8s %8s\n\n", "SCL", "SCP", "SCD");
  }
  for(i = 0; i < num_of_periods; i++) {
    if(verbose>1) {
      printf ("%4d%9.2f%9.2f", (period_length*i)+1, Alpha[i]*PESum[i], Beta[i]*PRSum[i]);
      printf ("%9.2f%9.2f", Gamma[i]*PROSum[i], Delta[i]*PLSum[i]); 
      printf ("%9.2f%9.2f\n", SPhat[i], DSAct[i]); 
    }
    // D becomes the mean of D_sum
    /* SG 6/5/06: changed totalyears to nCalibrationYears to support
     **            user defined calibration intervals. When not used
     **            nCalibrationYears==totalyears; hence no change then
     */
    D[i] = D_sum[i] / nCalibrationYears;
  }
  // The files are closed 
  // fclose(fin);
  // fclose(fout);
}// end of Calcd()

//-----------------------------------------------------------------------------
// This function uses previously calculated sums to find K, which is the 
// original weighting factor used in the Palmer Index calculation.
//-----------------------------------------------------------------------------
void pdsi::CalcOrigK_R() {
  int month, year;
  number sums;        //used to calc k
  float dtemp;
  DKSum = 0;
  
  // FILE * inputd; // File pointer for input file dvalue
  // The dvalue file is open for reading. 
  // if((inputd=fopen("dvalue","r"))==NULL) {
  //   if(verbose>0)
  //     printf ("Error reading the file with d values.\n");
  //   exit(1);
  // }
  // Calculate k, which is K', or Palmer's second approximation of K
  for(int per = 0; per < num_of_periods; per++){
    if(PSum[per] + LSum[per] == 0)
      sums = 0;//prevent div by 0
    else
      sums = (PESum[per] + RSum[per] + ROSum[per]) / (PSum[per] + LSum[per]);
    
    if(D[per] == 0)
      k[per] = 0.5;//prevent div by 0
    else
      k[per] = (1.5) * log10((sums + 2.8) / D[per]) + 0.5;
    
    DKSum += D[per]*k[per];
  } 
  
  
  // Set duration factors to Palmer's original duration factors
  drym = .309;
  dryb = 2.691;
  
  wetm = drym;
  wetb = dryb;
  
  
  // Initializes the book keeping indices used in finding the PDSI
  Prob = 0.0;
  X1 = 0.0;
  X2 = 0.0;
  X3 = 0.0;
  X = 0.0;
  V = 0.0;
  Q = 0.0;
  
  // open file point to bigTable.tbl if necessary
  // FILE* table;
  // if(extra == 2 || extra == 9){
  //   table = fopen("bigTable.tbl","w");
  //   if(table == NULL){
  //     if(verbose > 0) 
  //       printf("Error opening file \"bigTable.tbl\"\n");
  //   }
  //   else {
  //     //write column headers 
  // 
  //       // fprintf(table, "YEAR  MONTH     Z     %Prob     "); 
  //       // fprintf(table, "X1       X2      X3\n"); 
  //     
  //   }
  // }
  // else
  // table = NULL;
  
  
  // Reads in all previously calclulated d values and calculates Z
  // then calls CalcX to compute the corresponding PDSI value
  for(int i = 0; i < dvalue_d.length(); i++){
    // while((fscanf(inputd,"%d %d %f", &year, &month, &dtemp))!=EOF) {
    month = (i % num_of_periods) + 1;
    year = 1900;
    PeriodList.insert(month);
    YearList.insert(year);
    d = dvalue_d[i];
    month--;
    K = (17.67/DKSum) * k[month];
    if(d != MISSING)
      Z = d*K;
    else
      Z = MISSING;
    
    ZIND.insert(Z);
    CalcOneX_R(month,year);
  }
  // fclose(inputd);
  // if(table)
  //   fclose(table);
  // Now that all calculations have been done they can be output to the screen
  if(verbose>1) {
    int i;
    if(Weekly)   
      printf ("STATION = %5d %24c PARAMETERS AND MEANS OF WEEKLY VALUE FOR %d YEARS\n\n", 0, ' ', totalyears);
    else
      printf ("STATION = %5d %24c PARAMETERS AND MEANS OF MONTHLY VALUE FOR %d YEARS\n\n", 0, ' ', totalyears);      
    printf ("%4s %8s %8s %8s %8s %8s %7s %8s", "MO", "ALPHA", "BETA", "GAMMA", "DELTA", "K", "P", "S");
    printf ("%9s %8s %8s %8s %8s %8s %8s\n\n", "PR", "PE", "PL", "ET", "R", "L", "RO");
    for (i=0;i<num_of_periods;i++) {
      printf ("%4d %8.4f %8.4f %8.4f %8.4f", (period_length*i)+1, Alpha[i], Beta[i], Gamma[i], Delta[i]);
      printf ("%9.3f %8.2f %8.2f %8.2f", 17.67/DKSum*k[i], PSum[i]/totalyears, PROSum[i]/totalyears, PRSum[i]/totalyears);
      printf ("%9.2f %8.2f %8.2f %8.2f", PESum[i]/totalyears, PLSum[i]/totalyears, ETSum[i]/totalyears, RSum[i]/totalyears);
      printf ("%9.2f %8.2f\n", LSum[i]/totalyears, ROSum[i]/totalyears);
    }
    printf ("\n\n\n%4s %8s %8s %8s %8s %8s\n\n", "PER", "D-ABS", "SIG-D", "DEP", "S-DEP", "SIG-S");
    for (i=0;i<num_of_periods;i++) {
      printf ("%4d %8.3f %8.2f %8.2f ", (period_length*i)+1, D[i], sqrt(DSSqr[i]/(totalyears-1)), DEPSum[i]/totalyears);
      if (i==7) {
        number E, DE;
        E=SD/totalyears;
        DE=sqrt((SD2-E*SD)/(totalyears-1));
        printf ("%8.2f %8.2f", E, DE);
      }
      printf ("\n");
    }
  }
}//end of CalcOrigK_R()

//-----------------------------------------------------------------------------
// This function calculates X, X1, X2, and X3
//
// X1 = severity index of a wet spell that is becoming "established"
// X2 = severity index of a dry spell that is becoming "established"
// X3 = severity index of any spell that is already "established"
//
// newX is the name given to the pdsi value for the current week.
// newX will be one of X1, X2 and X3 depending on what the current 
// spell is, or if there is an established spell at all.
//-----------------------------------------------------------------------------
void pdsi::CalcOneX_R(int period_number, int year) {
  number newV;    //These variables represent the values for 
  number newProb; //corresponding variables for the current period.
  number newX;    //They are kept seperate because many calculations
  number newX1;   //depend on last period's values.  
  number newX2;
  number newX3;
  number ZE;      //ZE is the Z value needed to end an established spell
  
  number m, b, c;
  
  flag wd;        //wd is a sign changing flag.  It allows for use of the same
  //equations during both a wet or dry spell by adjusting the
  //appropriate signs.
  
  if(X3>=0){
    m = wetm;
    b = wetb;
  }
  else{
    m = drym;
    b = dryb;
  }
  c = 1 - (m / (m + b));
  
  
  if(Z != MISSING){
    // This sets the wd flag by looking at X3
    if(X3>=0) wd=1;
    else wd=-1;
    // If X3 is 0 then there is no reason to calculate Q or ZE, V and Prob
    // are reset to 0;
    if(X3==0) {
      newX3=0;
      newV=0;
      newProb=0;
      ChooseX(newX, newX1, newX2, newX3, bug);
    }
    // Otherwise all calculations are needed.
    else {
      newX3 = (c * X3 + Z/(m+b));
      ZE = (m+b)*(wd*0.5 - c*X3);
      Q=ZE+V;  
      newV = Z - wd*(m*0.5) + wd*min(wd*V+tolerance,0);
      
      if((wd*newV)>0) {
        newV=0;
        newProb=0;
        newX1=0;
        newX2=0;
        newX=newX3;
        while(!altX1.is_empty())
          altX1.head_remove();
        while(!altX2.is_empty())
          altX2.head_remove();
      }
      else {
        newProb=(newV/Q)*100;
        if(newProb>=100-tolerance) {
          newX3=0;
          newV=0;
          newProb=100;
        }
        ChooseX(newX, newX1, newX2, newX3, bug);
      }
    }
    
    // if(table != NULL){
    //   //output stuff to a table
    //   //year, period, z, newProb, newX1, newX2, newX3
    //   fprintf(table, "%5d %5d %7.2f %7.2f ",year,period_number,Z,newProb);
    //   fprintf(table, "%7.2f %7.2f %7.2f\n",newX1, newX2, newX3);
    // }
    
    //update variables for next month:
    V = newV;
    Prob = newProb;
    X1 = newX1;
    X2 = newX2;
    X3 = newX3;
    
    //add newX to the list of pdsi values
    Xlist.insert(newX);
    XL1.insert(X1);
    XL2.insert(X2);
    XL3.insert(X3);
    ProbL.insert(Prob);
  }
  else{
    //This month's data is missing, so output MISSING as PDSI.  
    //All variables used in calculating the PDSI are kept from 
    //the previous month.  Only the linked lists are changed to make
    //sure that if backtracking occurs, a MISSING value is kept 
    //as the PDSI for this month.
    
    // if(table != NULL){
    //   //output stuff to a table
    //   //year, period, z, newProb, newX1, newX2, newX3
    //   fprintf(table, "%5d %5d %7.2f %7.2f ",year,period_number,Z,MISSING);
    //   fprintf(table, "%7.2f %7.2f %7.2f\n",MISSING, MISSING, MISSING);
    // }
    Xlist.insert(MISSING);
    XL1.insert(MISSING);
    XL2.insert(MISSING);
    XL3.insert(MISSING);
    ProbL.insert(MISSING);
  }
  
}//end of CalcOneX_R

//-----------------------------------------------------------------------------
// This function calculates the sum of all actual and potential values for each
// indivual period over the given period of years.  To do this it must calculate
// these values for each period in the period.  The potential values are then
// stored for future use.
//-----------------------------------------------------------------------------
void pdsi::SumAll_R() {
  // FILE * fout;
  // char Temp[150], Precip[150];
  int actyear;
  number DEP=0;
  SD=0;
  SD2=0;
  int nCalibrationPeriodsLeft = nCalibrationPeriods; /* init periods left */
  
  // Initializes the sums to 0;
  for(int i = 0; i < 52; i++) {
    ETSum[i] = 0;
    RSum[i] = 0;
    LSum[i] = 0;
    ROSum[i] = 0;
    PSum[i] = 0;
    PESum[i] = 0;
    PRSum[i] = 0;
    PLSum[i] = 0;
    PROSum[i] = 0;
  }
  
  // Opens the temperature and precipitation files for reading the input and
  // potentials file for temparary storage.
  // fout = fopen("potentials","w");
  // 
  // // strcpy(Temp,"monthly_T");
  // // strcpy(Precip,"monthly_P");
  // 
  // // input_temp = fopen(Temp,"r");
  // // input_prec = fopen(Precip,"r");
  // 
  // //print column headers:
  // fprintf(fout," Year  MONTH       P         PE         PR         PRO");
  // fprintf(fout,"        PL        P-PE \n");
  
  
  // This loop runs to read in and calculate the values for all years
  int this_y = 0;
  int this_y_start;
  int year = 1;
  for(IntegerVector::iterator it = all_years.begin(); it != all_years.end(); ++it) {
    actyear = *it;
    year++;
    this_y_start = this_y;
    // cout << actyear;
    // cout << "\n";
    
    for(int per = 0; per < num_of_periods; per++) {
      P[per] = monthly_precip[this_y];
      T[per] = monthly_temp[this_y];
      this_y++;
    }
    
    // This loop runs for each per in the year
    for(int per = 0; per < num_of_periods; per++) {
      if(P[per] >= 0 && T[per] != MISSING){
        // calculate the Potential Evapotranspiration first
        // because it's needed in later calculations
        
        // cout << "Year: ";
        // cout << actyear;
        // cout << ", ";
        // 
        // cout << "Month: ";
        // cout << per;
        // cout << ", ";
        // 
        // cout << "Precipitation: ";
        // cout << P[per];
        // cout << ", ";
        // 
        // cout << "Temperature: ";
        // cout << T[per];
        // cout << "\n";
        // 
        CalcMonPE(per,actyear);
        
        CalcPR();         // calculate Potential Recharge, Potential Runoff,
        CalcPRO();        // and Potential Loss
        CalcPL();
        CalcActual(per);  // Calculate Evapotranspiration, Recharge, Runoff,
        // and Loss
        
        // Calculates some statistical variables for output 
        // to the screen in the most verbose mode (verbose > 1)
        
        if (per > 4 && per < 8) {
          DEP = DEP + P[per] + L - PE;
          if (per == 7) {
            SD=SD+DEP;
            SD2=SD2+DEP*DEP;
            DEP=0;
          }
        }
        
        if (year > nStartYearsToSkip) {  
          if (nCalibrationPeriodsLeft > 0) { /* Within the calibration interval, so update Sums */
            // Reduce number of calibration years left by one for each year actually summed 
            nCalibrationPeriodsLeft--; 
            
            // Update the sums by adding the current water balance values
            ETSum[per] += ET;
            RSum[per] += R;
            ROSum[per] += RO;
            LSum[per] += L;
            PSum[per] += P[per];
            PESum[per] += PE;
            PRSum[per] += PR;
            PROSum[per] += PRO;
            PLSum[per] += PL;
          }
        }
        
        //P, PE, PR, PRO, and PL will be used later, so 
        //these variables need to be stored to an outside file
        //where they can be accessed at a later time
        // fprintf(fout,"%5d %5d %10.6f ",actyear,(period_length*per)+1,P[per]);
        // fprintf(fout,"%10.6f %10.6f %10.6f %10.6f ",PE,PR,PRO,PL);
        // fprintf(fout,"%10.6f\n",P[per]-PE);
        
        potentials_P[this_y_start] = P[per];
        potentials_PE[this_y_start] = PE;
        potentials_PR[this_y_start] = PR;
        potentials_PRO[this_y_start] = PRO;
        potentials_PL[this_y_start] = PL;
        potentials_PPE[this_y_start] = P[per]-PE;
        
      }//matches if(P[per]>= 0 && T[per] != MISSING)
      else {
        // fprintf(fout,"%5d %5d %f ",actyear, (period_length*per)+1,MISSING);
        // fprintf(fout,"%10.6f %10.6f %10.6f ", MISSING, MISSING, MISSING);
        // fprintf(fout,"%10.6f %10.6f\n",MISSING, MISSING);
        
        potentials_P[this_y_start] = MISSING;
        potentials_PE[this_y_start] = MISSING;
        potentials_PR[this_y_start] = MISSING;
        potentials_PRO[this_y_start] = MISSING;
        potentials_PL[this_y_start] = MISSING;
        potentials_PPE[this_y_start] = MISSING;
        
      }
      
      this_y_start++;
    }//end of period loop
  }//end of year loop
  
  // We are done with these files for now so close them
  // fclose(fout);
  // fclose(input_temp);
  // fclose(input_prec);
}

//-----------------------------------------------------------------------------
// CalcPR calculates the Potential Recharge of the soil for one period of the 
// year being examined.  PR = Soils Max Capacity - Soils Current Capacity or
// AWC - (SU + Ss)
//-----------------------------------------------------------------------------
void pdsi::CalcPR() {
  PR = AWC - (Su + Ss);
}
//-----------------------------------------------------------------------------
// CalcPRO calculates the Potential Runoff for a given period of the year being
// examined.  PRO = Potential Precip - PR. Palmer arbitrarily set the Potential
// Precip to the AWC making PRO = AWC - (AWC - (Su + Ss)). This then simplifies
// to PRO = Su + Ss
//-----------------------------------------------------------------------------
void pdsi::CalcPRO() {
  PRO = Ss + Su;
}
//-----------------------------------------------------------------------------
// CalcPL calculates the Potential Loss of moisture in the soil for a period of
// one period of the year being examined. If the Ss capacity is enough to
// handle all PE, PL is simple PE.  Otherwise, potential loss from Su occurs at
// the rate of (PE-Ss)/AWC*Su.  This means PL = Su*(PE - Ss)/AWC + Ss
//-----------------------------------------------------------------------------
void pdsi::CalcPL() {
  if(Ss >= PE)
    PL = PE;
  else {
    PL = ((PE - Ss) * Su) / (AWC) + Ss;
    if(PL > PRO)  // If PL>PRO then PL>water in the soil.  This isn't
      PL = PRO;   // possible so PL is set to the water in the soil
  }
}
//-----------------------------------------------------------------------------
// CalcActual calculates the actual values of evapotranspiration,soil recharge,
// runoff, and soil moisture loss.  It also updates the soil moisture in both
// layers for the next period depending on current weather conditions.
//-----------------------------------------------------------------------------
void pdsi::CalcActual(int per) {
  number R_surface = 0.0;   // recharge of the surface layer
  number R_under = 0.0;    // recharge of the underlying layer
  number surface_L = 0.0;   // loss from surface layer
  number under_L = 0.0;    // loss from underlying layer
  number new_Su, new_Ss;    // new soil moisture values
  
  
  if(P[per] >= PE) {
    // The precipitation exceeded the maximum possible evapotranspiration
    // (excess moisture)
    ET = PE;   // Enough moisture for all potential evapotranspiration to occur
    L = 0.0;   // with no actual loss of soil moisture
    
    if((P[per] - PE) > (1.0 - Ss)) {
      // The excess precip will recharge both layers. Note: (1.0 - SS) is the 
      // amount of water needed to saturate the top layer of soil assuming it
      // can only hold 1 in. of water.
      R_surface = 1.0 - Ss;
      new_Ss = 1.0;
      
      if((P[per] - PE - R_surface) < ((AWC - 1.0) - Su)) {
        // The entire amount of precip can be absorbed by the soil (no runoff)
        // and the underlying layer will receive whats left after the top layer
        // Note: (AWC - 1.0) is the amount able to be stored in lower layer
        R_under = (P[per] - PE - R_surface);
        RO = 0.0;
      }
      else {
        // The underlying layer is fully recharged and some runoff will occur
        R_under = (AWC - 1.0) - Su;
        RO = P[per] - PE - (R_surface + R_under);
      }
      new_Su = Su + R_under;
      R = R_surface + R_under;//total recharge
    }
    else {
      // There is only enough moisture to recharge some of the top layer.
      R = P[per] - PE;
      new_Ss = Ss + R;
      new_Su = Su;
      RO = 0.0;
    }
  }// End of if(P[per] >= PE)
  else {
    // The evapotranspiration is greater than the precipitation received.  This
    // means some moisture loss will occur from the soil.
    if(Ss > (PE - P[per])) {
      // The moisture from the top layer is enough to meet the remaining PE so 
      // only the top layer losses moisture.
      surface_L = PE - P[per];
      under_L = 0.0;
      new_Ss = Ss - surface_L;
      new_Su = Su;
    }
    else {
      // The top layer is drained, so the underlying layer loses moisture also.
      surface_L = Ss;
      under_L = (PE - P[per] - surface_L) * Su / AWC;
      if(Su < under_L)
        under_L = Su;
      new_Ss = 0.0;
      new_Su = Su - under_L;
    }
    R = 0;// No recharge occurs
    L = under_L + surface_L;// Total loss
    RO = 0.0;// No extra moisture so no runoff
    ET = P[per] + L;// Total evapotranspiration
  }
  Ss = new_Ss;//update soil moisture values
  Su = new_Su;
}//end of CalcActual(int per)

//-----------------------------------------------------------------------------
// This function calculates Alpha, Beta, Gamma, and Delta, the normalizing
// climate coefficients in the water balance equation.
// If the user desires, the results are output to the screen and a file.
//-----------------------------------------------------------------------------
void pdsi::CalcWBCoef_R() {
  
  // FILE *wb;
  
  // The coefficients are calculated by per
  for (int per=0; per < num_of_periods; per++) {
    
    //calculate alpha:
    if(PESum[per] != 0.0)
      Alpha[per] = ETSum[per] / PESum[per];
    else if(ETSum[per] == 0.0)
      Alpha[per] = 1.0;
    else
      Alpha[per] = 0.0;
    
    //calculate beta:
    if(PRSum[per] != 0.0)
      Beta[per] = RSum[per] / PRSum[per];
    else if(RSum[per] == 0.0)
      Beta[per] = 1.0;
    else
      Beta[per] = 0.0;
    
    //calculate gamma:
    if(PROSum[per] != 0.0)
      Gamma[per] = ROSum[per] / PROSum[per];
    else if(ROSum[per] == 0.0)
      Gamma[per] = 1.0;
    else
      Gamma[per] = 0.0;
    
    //calculate delta:
    if(PLSum[per] != 0.0)
      Delta[per] = LSum[per] / PLSum[per];
    else 
      Delta[per] = 0.0;
  }
  
  // if(extra==1 || extra == 9){
  //   //output water balance coefficients
  //   wb = fopen("WB.tbl", "w");
  //   fprintf(wb, "PERIOD   ALPHA     BETA    GAMMA    DELTA\n");
  //   if(verbose>1)
  //     printf ("\nPERIOD   ALPHA     BETA    GAMMA    DELTA\n");
  //   for(int i = 0; i < num_of_periods; i++){
  //     fprintf(wb, "%3d %10.4f %8.4f %8.4f %8.4f \n", (period_length*i)+1, Alpha[i], Beta[i], Gamma[i], Delta[i]);
  //     if(verbose>1)
  //       printf ("%3d %10.4f %8.4f %8.4f %8.4f \n", (period_length*i)+1, Alpha[i], Beta[i], Gamma[i], Delta[i]);
  //   }
  //   fclose(wb);
  // }
}//end CalcWBCoef()


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void pdsi::Backtrack(number X1, number X2) {
  number num1,num2;
  node * ptr=NULL;
  num1=X1;
  while (!altX1.is_empty() && !altX2.is_empty()) {
    if (num1>0) {
      num1=altX1.head_remove();
      num2=altX2.head_remove();
    }
    else {
      num1=altX2.head_remove();
      num2=altX1.head_remove();
    }
    if (-tolerance<=num1 && num1<=tolerance) num1=num2;
    ptr=Xlist.set_node(ptr,num1);
  }
}//end of backtrack()
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void pdsi::ChooseX(number& newX, number& newX1, number& newX2, number& newX3, int bug)
{
  number m, b;
  number wetc, dryc;
  
  if(X3>=0){
    m = wetm;
    b = wetb;
  }
  else{
    m = drym;
    b = dryb;
  }
  
  wetc = 1 - (wetm / (wetm+wetb));
  dryc = 1 - (drym / (drym+wetb));
  
  newX1 = (wetc*X1 + Z/(wetm+wetb));
  if(newX1 < 0)
    newX1 = 0;
  newX2 = X2;
  
  if(bug==0){
    newX2 = (dryc*X2 + Z/(drym+dryb));
    if(newX2 > 0)
      newX2 = 0;
  }
  
  if((newX1 >= 0.5)&&(newX3 == 0)){
    Backtrack(newX1, newX2);
    newX = newX1;
    newX3 = newX1;
    newX1 = 0;
  }
  else{
    newX2 = (dryc*X2 + Z/(drym+dryb));
    if(newX2 > 0)
      newX2 = 0;
    
    if((newX2 <= -0.5)&&(newX3 == 0)){
      Backtrack(newX2, newX1);
      newX = newX2;
      newX3 = newX2;
      newX2 = 0;
    }
    else if(newX3 == 0) {
      if(newX1 == 0){
        Backtrack(newX2, newX1);
        newX = newX2;
      }
      else if(newX2 == 0){
        Backtrack(newX1, newX2);
        newX = newX1;
      }
      else{
        altX1.insert(newX1);
        altX2.insert(newX2);
        newX = newX3;
      }
    }
    
    else{
      //store X1 and X2 in their linked lists for possible use later
      altX1.insert(newX1);
      altX2.insert(newX2);
      newX = newX3;
    }
  }
}//end of chooseX
void pdsi::CalcDurFact(number &slope, number &intercept, int sign){
  //calculates m and b, which are used to calculated X(i) 
  //based on the Z index.  These constants will determine the 
  //weight that the previous PDSI value and the current Z index 
  //will have on the current PDSI value.  This is done by finding 
  //several of the driest periods at this station and assuming that
  //those periods represents an extreme drought.  Then a linear 
  //regression is done to determine the relationship between length
  //of a dry (or wet) spell and the accumulated Z index during that
  //same period.  
  //
  //it appears that there needs to be a different weight given to 
  //negative and positive Z values, so the variable 'sign' will 
  //determine whether the driest or wettest periods are looked at.
  
  int num_list = 10;
  number sum[10];
  int length[10];
  int i;
  
  if(Weekly){
    if(period_length==1){
      length[0]=13;
      length[1]=26;
      length[2]=39;
      length[3]=52;
      length[4]=78;
      length[5]=104;
      length[6]=130;
      length[7]=156;
      length[8]=182;
      length[9]=208;
    }
    else if(period_length==2){
      length[0]=6;
      length[1]=13;
      length[2]=19;
      length[3]=26;
      length[4]=39;
      length[5]=52;
      length[6]=65;
      length[7]=78;
      length[8]=91;
      length[9]=104;
    }
    else if(period_length==4){
      length[0]=3;
      length[1]=6;
      length[2]=10;
      length[3]=13;
      length[4]=20;
      length[5]=26;
      length[6]=33;
      length[7]=39;
      length[8]=46;
      length[9]=52;
    }
    else if(period_length==13){
      length[0]=2;
      length[1]=3;
      length[2]=4;
      length[3]=5;
      length[4]=6;
      length[5]=8;
      length[6]=10;
      length[7]=12;
      length[8]=14;
      length[9]=16;
    }
  }
  else{
    length[0]=3;
    length[1]=6;
    length[2]=9;
    length[3]=12;
    length[4]=18;
    length[5]=24;
    length[6]=30;
    length[7]=36;
    length[8]=42;
    length[9]=48;
  }
  
  
  for(i = 0; i < num_list; i++){
    sum[i] = get_Z_sum(length[i],sign);
    //printf("sum[%d] = %f\n",i,sum[i]);
  }
  if(verbose > 1){
    printf("Points used in linear regression for Duration Factors:\n");
    for(i=0;i<num_list;i++)
      printf("%7d  ",length[i]);
    printf("\n");
    for(i=0;i<num_list;i++)
      printf("%7.2f  ",sum[i]);
    printf("\n");
  }
  
  LeastSquares(length, sum, num_list, sign, slope, intercept);
  
  //printf("the line is: y = %f * x + %f\n",slope,intercept);
  
  //now divide m and b by 4 or -4 becuase that line represents
  //pdsi of either 4.0 or -4.0
  
  slope = slope / (sign*4);
  intercept = intercept / (sign*4);
}//end of CalcDurFact()

number pdsi::get_Z_sum(int length, int sign) {
  number sum, max_sum, z;
  llist tempZ, list_to_sum, list_of_sums;
  
  number highest_reasonable;
  number percentile;
  number reasonable_tol = 1.25;
  /* SG 6/5/06: Add variable to implement user-defined calibration interval */
  int nCalibrationPeriodsLeft;
  
  copy(tempZ,ZIND); 
  sum = 0;
  
  /* SG 6/5/06: Remove the periods from the list until we get to the
   **            start of the calibration interval
   */
  for (int i=0; (i < nStartPeriodsToSkip) && (!tempZ.is_empty()) ; i++) {
    tempZ.tail_remove(); /* remove periods before the start of the interval */
  }
  /* SG 6/5/06: We now have a list that begins at the calibration interval.
   **            However, if the list has more periods than the length of the
   **            calibration interval, we must be sure to not go past the
   **            calibration interval length
   */
  nCalibrationPeriodsLeft = nCalibrationPeriods; /* init periods left */
  //first fill the list to be summed  
  for(int i = 0; i < length; i++){
    if(tempZ.is_empty()){
      printf("Error: tempZ is empty.\n");
      i = length;
    }
    else {
      z = tempZ.tail_remove();
      nCalibrationPeriodsLeft--; /* reduce by one period for each remove */
  /* assumes that nCalibrationPeriods is >= length, reasonable 
   ** This is a reasonable assumption and does not hurt if
   ** anything if false--just calibrate over a slightly longer
   ** interval, which is already way too short in that case */
  if(z != MISSING){
    sum += z;
    list_to_sum.insert(z);
  }
  else{
    i--;
  }
    }
    //printf("i = %d  z= %f  sum = %f\n",i,z,sum);
  }
  
  
  //now for each remaining Z value,
  //recalculate the sum based on last value in the
  //list to sum and the next Z value
  max_sum = sum;
  list_of_sums.insert(sum);
  while(!tempZ.is_empty() && nCalibrationPeriodsLeft > 0){
    z = tempZ.tail_remove();
    nCalibrationPeriodsLeft--; /* reduce by one period for each remove */
  if(z != MISSING){
    sum -= list_to_sum.tail_remove();
    sum += z;
    list_to_sum.insert(z);
    list_of_sums.insert(sum);
  }
  if(sign * sum > sign * max_sum)
    max_sum = sum;
  }
  
  //highest reasonable is the highest (or lowest)
  //value that is not due to some freak anomaly in the
  //data.
  //"freak anomaly" is defined as a value that is either
  //   1) 25% higher than the 98th percentile
  //   2) 25% lower than the 2nd percentile
  //
  highest_reasonable = 0; 
  if(sign == 1)
    percentile = list_of_sums.safe_percentile(.98);
  if(sign == -1)
    percentile = list_of_sums.safe_percentile(.02);
  
  while(!list_of_sums.is_empty()){
    sum = list_of_sums.tail_remove();
    if(sign * sum > 0 ){
      if( (sum / percentile) < reasonable_tol ) {
        if(sign * sum > sign * highest_reasonable )
          highest_reasonable = sum;
      }
    }
  }
  
  if(sign == -1)
    return max_sum;
  else if(sign == 1)
    //return max_sum;
    return highest_reasonable;
  else
    return MISSING;
}//end of get_Z_sum()

void pdsi::LeastSquares(int *x, number *y, int n, int sign, number &slope, number &intercept) {
  number sumX, sumX2, sumY, sumY2, sumXY;
  number SSX, SSY, SSXY;
  number xbar, ybar;
  
  number correlation = 0;
  number c_tol = 0.85;
  
  number max = 0;
  number max_diff = 0;
  int max_i = 0;
  
  number this_x, this_y;
  int i;
  
  sumX = 0; sumY = 0; sumX2 = 0; sumY2 = 0; sumXY = 0;
  for(i = 0; i < n; i++){
    this_x = x[i];
    this_y = y[i];
    
    sumX += this_x;
    sumY += this_y;
    sumX2 += this_x * this_x;
    sumY2 += this_y * this_y;
    sumXY += this_x * this_y;
  }
  
  xbar = sumX/n;
  ybar = sumY/n;
  
  SSX = sumX2 - (sumX * sumX)/n;
  SSY = sumY2 - (sumY * sumY)/n;
  SSXY = sumXY - (sumX * sumY)/n;
  
  correlation = SSXY / (sqrt(SSX) * sqrt(SSY));
  
  if(verbose > 1 && (sign*correlation) < c_tol ){
    printf("original correlation = %.4f \n",correlation);
  }
  
  i = n - 1;
  while((sign*correlation) < c_tol && i > 3){
    //when the correlation is off, it appears better to 
    //take the earlier sums rather than the later ones.
    this_x = x[i];
    this_y = y[i];
    
    sumX -= this_x;
    sumY -= this_y;
    sumX2 -= this_x * this_x;
    sumY2 -= this_y * this_y;
    sumXY -= this_x * this_y;
    
    SSX = sumX2 - (sumX * sumX)/i;
    SSY = sumY2 - (sumY * sumY)/i;
    SSXY = sumXY - (sumX * sumY)/i;
    
    xbar = sumX/i;
    ybar = sumY/i;
    
    correlation = SSXY / (sqrt(SSX) * sqrt(SSY));
    i--;
  }
  
  if(verbose > 1){
    printf("final correlation =  %.4f\n\n",correlation);
  }
  slope = SSXY / SSX;
  
  n = i+1;
  for(int i = 0; i < n; i++){
    if(sign*(y[i] - slope * x[i]) > sign*max_diff){
      max_diff = y[i] - slope * x[i];
      max_i = i;
      max = y[i];
    }
  }
  intercept = max - slope*x[max_i];
}//end of LeastSquares()

number pdsi::getPDSI(int period, int year) {
  return getValue(Xlist, period, year);
}
number pdsi::getZIND(int period, int year) {
  return getValue(ZIND, period, year);
}
number pdsi::getWPLM(int period, int year) {
  
  number x1,x2,x3,p,wp;
  
  x1 = getValue(XL1, period, year);
  x2 = getValue(XL2, period, year);
  x3 = getValue(XL3, period, year);
  p = getValue(ProbL, period, year);
  if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING)  
    wp = MISSING;  
  else{  
    p = p / 100;
    if (x3==0) {
      // There is not an established wet or dry spell so PHDI = PDSI (ph=x) 
      // and the WPLM value is the maximum absolute value of X1 or X2
      wp=x1;
      if (-x2>(x1+tolerance))
        wp=x2;
    }
    else if (p>(0+tolerance/100) && p<(1-tolerance/100)) {
      // There is an established spell but there is a possibility it has or is
      // ending.  The WPLM is then a weighted average between X3 and X1 or X2
      if (x3 < 0) 
        // X3 is negative so WPLM is weighted average of X3 and X1
        wp=(1-p)*x3 + p*x1;
      else
        // X3 is positive so WPLM is weighted average of X3 and X2
        wp=(1-p)*x3 + p*x2;
    }
    else
      // There is an established spell without possibility of end meaning the
      // WPLM is simply X3
      wp=x3;
  }
  return wp;
}

number pdsi::getPHDI(int period, int year) {
  number x, x3;
  
  x = getValue(Xlist, period, year);
  x3 = getValue(XL3, period, year);
  if(x == MISSING || x3 == MISSING) 
    return MISSING;
  if(x3==0)
    return x;
  else
    return x3;
}

number pdsi::getValue(llist &List, int period, int year) {
  llist tempPer, tempYear, tempList;
  number per, yr, val;
  bool loop_exit = false;
  
  copy(tempList, List);
  copy(tempPer, PeriodList);
  copy(tempYear, YearList);
  
  while(! loop_exit) {
    if(tempList.is_empty())
      loop_exit = true;
    if(tempPer.is_empty())
      loop_exit = true;
    if(YearList.is_empty())
      loop_exit = true;
    
    val = tempList.head_remove();
    per = tempPer.head_remove();
    yr = tempYear.head_remove();
    
    if(yr == year && per == period)
      return val;
  }
  return MISSING;
}

number* pdsi::getYearArray(int &size) {
  size = YearList.get_size();
  return YearList.returnArray();
}
number* pdsi::getPerArray(int &size) {
  size = PeriodList.get_size();
  return PeriodList.returnArray();
}

number* pdsi::getPDSIArray(int &size) {
  size = Xlist.get_size();
  return Xlist.returnArray();
}

NumericVector pdsi::getPDSIArray_R() {
  int size = Xlist.get_size();
  number* A = Xlist.returnArray();
  
  Rcpp::NumericVector vectorOutput(size);
  for(int i=0;i<size;i++)
  {
    vectorOutput(i) = A[i];
  }
  
  return vectorOutput;
}

number* pdsi::getPDSIArray(int start_per, int start_yr, 
                           int end_per, int end_yr, int &size) {
  return getSubArray(Xlist, start_per, start_yr, end_per, end_yr, size);
}
number* pdsi::getZINDArray(int &size){
  size = ZIND.get_size();
  return ZIND.returnArray();
}
number* pdsi::getZINDArray(int start_per, int start_yr, 
                           int end_per, int end_yr, int &size) {
  return getSubArray(ZIND, start_per, start_yr, end_per, end_yr, size);
}
number* pdsi::getPHDIArray(int &size) {
  number *x, *x3;
  x = Xlist.returnArray();
  x3 = XL3.returnArray();
  size = Xlist.get_size();
  number *A = new number[size];
  if(A == NULL){
    size = 0; 
    return A;
  }
  for(int i = 0; i < size; i++){
    if(x[i] != MISSING){
      if(x3[i]==0)
        A[i] = x[i];
      else
        A[i] = x3[i];
    }
    else
      A[i] = MISSING;
  }
  delete [] x;
  delete [] x3;
  return A;
}
number* pdsi::getPHDIArray(int start_per, int start_yr, 
                           int end_per, int end_yr, int &size) {
  number *x, *x3; 
  int tempsize;
  x = getSubArray(Xlist, start_per, start_yr, end_per, end_yr, size);
  x3 = getSubArray(XL3, start_per, start_yr, end_per, end_yr, tempsize); 
  number *A = new number[size];
  if(A == NULL){
    size = 0;
    return A;
  }
  for(int i = 0; i < size; i++){ 
    if(x[i] != MISSING){ 
      if(x3[i]==0) 
        A[i] = x[i]; 
      else 
        A[i] = x3[i]; 
    } 
    else 
      A[i] = MISSING; 
  }
  delete [] x;
  delete [] x3;
  return A;
}
number* pdsi::getWPLMArray(int &size) {
  number *A;
  number *x1Array,*x2Array,*x3Array,*pArray; 
  number x1, x2, x3, p, wp;
  
  x1Array = XL1.returnArray();
  x2Array = XL2.returnArray();
  x3Array = XL3.returnArray();
  pArray = ProbL.returnArray();
  
  size = XL1.get_size();
  
  A = new number[size];
  if(A == NULL){
    size = 0;
    return A;
  }
  
  for(int i = 0; i < size; i++){
    x1 = x1Array[i];
    x2 = x2Array[i];
    x3 = x3Array[i];
    p = pArray[i];
    
    if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING)
      wp = MISSING;
    else{
      p = p / 100; 
      if (x3==0) { 
        // There is not an established wet or dry spell so PHDI = PDSI (ph=x)  
        // and the WPLM value is the maximum absolute value of X1 or X2 
        wp=x1; 
        if (-x2>(x1+tolerance)) 
          wp=x2; 
      } 
      else if (p>(0+tolerance/100) && p<(1-tolerance/100)) { 
        // There is an established spell but there is a possibility it has or is 
        // ending.  The WPLM is then a weighted average between X3 and X1 or X2 
        if (x3 < 0)  
          // X3 is negative so WPLM is weighted average of X3 and X1 
          wp=(1-p)*x3 + p*x1; 
        else 
          // X3 is positive so WPLM is weighted average of X3 and X2 
          wp=(1-p)*x3 + p*x2; 
      } 
      else 
        // There is an established spell without possibility of end meaning the 
        // WPLM is simply X3 
        wp=x3; 
    }
    A[i] = wp;
  }
  delete [] x1Array;
  delete [] x2Array;
  delete [] x3Array;
  delete [] pArray;
  return A;
}
number* pdsi::getWPLMArray(int start_per, int start_yr, 
                           int end_per, int end_yr, int &size) {
  number *A; 
  number *x1Array,*x2Array,*x3Array,*pArray;  
  number x1, x2, x3, p, wp; 
  int tempsize;
  
  x1Array = getSubArray(XL1, start_per, start_yr, end_per, end_yr, size);
  x2Array = getSubArray(XL2, start_per, start_yr, end_per, end_yr, tempsize);
  x3Array = getSubArray(XL3, start_per, start_yr, end_per, end_yr, tempsize);
  pArray = getSubArray(ProbL, start_per, start_yr, end_per, end_yr, tempsize);
  
  A = new number[size]; 
  if(A == NULL){ 
    size = 0; 
    return A; 
  } 
  
  for(int i = 0; i < size; i++){ 
    x1 = x1Array[i]; 
    x2 = x2Array[i]; 
    x3 = x3Array[i]; 
    p = pArray[i];
    
    if(x1 == MISSING || x2 == MISSING || x3 == MISSING || p == MISSING) 
      wp = MISSING; 
    else{ 
      p = p / 100; 
      if (x3==0) {  
        // There is not an established wet or dry spell so PHDI = PDSI (ph=x)   
        // and the WPLM value is the maximum absolute value of X1 or X2  
        wp=x1;  
        if (-x2>(x1+tolerance))  
          wp=x2;  
      }  
      else if (p>(0+tolerance/100) && p<(1-tolerance/100)) {  
        // There is an established spell but there is a possibility it has or is  
        // ending.  The WPLM is then a weighted average between X3 and X1 or X2  
        if (x3 < 0)   
          // X3 is negative so WPLM is weighted average of X3 and X1  
          wp=(1-p)*x3 + p*x1;  
        else  
          // X3 is positive so WPLM is weighted average of X3 and X2  
          wp=(1-p)*x3 + p*x2;  
      }  
      else  
        // There is an established spell without possibility of end meaning the  
        // WPLM is simply X3  
        wp=x3;  
    }
    A[i] = wp; 
  } 
  delete [] x1Array;
  delete [] x2Array;
  delete [] x3Array;
  delete [] pArray;
  return A; 
}
number* pdsi::getSubArray(llist &List, int start_per, int start_yr, 
                          int end_per, int end_yr, int &size) {
  
  llist temp;
  number *Array, *year, *period;
  int i,j;
  int cur_per, cur_yr;
  int per_len=0;
  int num_missing;
  
  Array = List.returnArray();
  year = YearList.returnArray();
  period = PeriodList.returnArray();
  
  for(j = 0; j < PeriodList.get_size(); j++){
    if(period[j] > per_len)
      per_len = (int)period[j];
  }   
  printf("per_len is: %d\n",per_len);
  printf("size of list is: %d\n",PeriodList.get_size() );
  
  if( (start_yr > year[0]) ||
      ( (start_yr == year[0]) && (start_per > period[0]) ) ) { 
    i = 0;
    while( ( ( year[i] < start_yr ) || 
           ( year[i] == start_yr && period[i] < start_per ) ) &&
           ( i < List.get_size() ) ) {
      i++;
    }
    while( ( ( year[i] < end_yr ) ||
           ( year[i] == end_yr && period[i] <= end_per ) ) &&
           ( i < List.get_size() ) ) {
      temp.insert(Array[i]);
      i++;
    }
    if(i == List.get_size()){
      cur_yr = (int)year[i-1];
      cur_per = (int)period[i-1];
      if((cur_per%per_len) == 0){
        cur_per = 1;
        cur_yr++;
      }
      else
        cur_per++;   
      while( (cur_yr < end_yr) ||
             ( (cur_yr == end_yr) && (cur_per <= end_per)) ) {
        temp.insert(MISSING);
        if((cur_per%per_len) == 0){
          cur_per = 1;
          cur_yr++;
        }
        else
          cur_per++;
      }
    } 
    
  }
  
  else {
    if(start_yr == year[0])
      num_missing = (int)period[0] - start_per;
    else{
      if(period[0] <= start_per)
        num_missing = ((int)year[0] - start_yr - 1)*per_len + ((int)period[0] - start_per + per_len);
      else
        num_missing = ((int)year[0] - start_yr)*per_len + ((int)period[0] - start_per);
    }
    printf("num_missing=%d\n",num_missing);
    for(j = 0; j < num_missing; j++)
      temp.insert(MISSING);
    
    i = 0;
    while( ( ( year[i] < end_yr ) ||
           ( year[i] == end_yr && period[i] <= end_per ) ) &&
           ( i < List.get_size() ) ) {
      temp.insert(Array[i]);
      printf("i=%d  cur_date: %d/%d  end_date: %d/%d\n",i,(int)period[i],(int)year[i],end_per,end_yr);
      i++;
    }
    if(i == List.get_size()){
      cur_yr = (int)year[i-1];
      cur_per = (int)period[i-1];
      if((cur_per%per_len) == 0){
        cur_per = 1;
        cur_yr++;
      }
      else
        cur_per++;   
      while( (cur_yr < end_yr) ||
             ( (cur_yr == end_yr) && (cur_per <= end_per)) ) {
        temp.insert(MISSING);
        printf("here i=%d  cur_date: %d/%d  end_date: %d/%d\n",i,cur_per,cur_yr,end_per,end_yr); 
        if((cur_per%per_len) == 0){
          cur_per = 1;
          cur_yr++;
        }
        else
          cur_per++;
      }
    }
  }
  delete [] Array;
  delete [] year;
  delete [] period;
  
  size = temp.get_size();
  return temp.returnArray();
  
}

inline int pdsi::is_int(char *string,int length) {
  int err=1; 
  for(int i=0; i<length; i++) 
    if(!isdigit(string[i])) err=0;
    return err; 
}
inline int pdsi::is_flt(char *string,int length) {
  int err=1;
  for(int i=0; i<length; i++)
    if(!isdigit(string[i]) && string[i]!='.') err=0;
    return err;
}

//-----------------------------------------------------------------------------
//**********   START OF FUNCTION DEFINITIONS FOR CLASS:  llist        *********
//-----------------------------------------------------------------------------
// The constructor for this function creates the sentinel for the
// double-linked list and points head to it.  The sentinel is a node that
// contains pointer to the beginning and to the end of the list.  Here both 
// of these pointers point back to the sentinel itself.
//-----------------------------------------------------------------------------
llist::llist() {
  head = new node;       // Creates the sentinel with head pointing to it
  head->next = head;     // Sets the sentinels front pointer to itself
  head->previous = head; // Sets the sentinels rear pointer to itself
  size = 0;
}
//-----------------------------------------------------------------------------
// The destructor for this function deallocates all of the memory
// for the entire list.  In order to do this it must move through
// the list and delete each of these nodes.  Then it deletes the
// sentinel.
//-----------------------------------------------------------------------------
llist::~llist() {
  node *mover;                // The temporary pointer to perform the work
  
  mover = head->next;         // mover starts at the node after the sentinel
  while (mover != head) {     // This loop occurs until mover has come complete
    // circle and is pointing at the sentinel
    mover = mover->next;      // mover becomes the next node
    delete mover->previous;   // The previous node is then deleted
  }
  delete mover;               // Finally the sentinel is deleted
}
//-----------------------------------------------------------------------------
// The insert function places a new node with the integer value x 
// between the sentinel and the first element in the list.  This
// effectively makes this new node the first element in the list.
//-----------------------------------------------------------------------------
void llist::insert(number x) {
  node *inserter;      // A new pointer to the node to be added
  inserter = new node; // This creates the node and points inserter to it
  inserter->key = x;   // The value in inserter is set to x
  
  inserter->next = head->next;        // The next field of the new node is
  // set to the node after the sentinel
  inserter->previous = head;          // The previous field is set to the 
  // sentinel
  inserter->next->previous = inserter;// The previous field of the node after
  // inserter is set to inserter
  inserter->previous->next = inserter;// The sentinels next field is set to 
  // inserter
  size++;                             // update size
}
int llist::get_size() {
  return size;
}

number* llist::returnArray() {
  
  node* cur;
  int i;
  number* A = new number[size];
  if(A != NULL){
    cur = head->previous;
    i = 0;
    while(cur != head){
      A[i] = cur->key;
      i++;
      cur=cur->previous;
    }
  }
  return A;
  
}
//-----------------------------------------------------------------------------
// The head_remove function removes the first node on the list and
// returns the value stored in that node
//-----------------------------------------------------------------------------
number llist::head_remove() {
  if(is_empty()) {
    return MISSING;
  }
  
  node *remover;
  number x;
  
  remover=head->next;
  x=remover->key;
  // First the previous field of the next node is set to head
  remover->next->previous = head;
  // Then the next field of head is set to the node after remover
  head->next = remover->next;
  size--;          //update size;
  delete remover;  // Finally the node can be deleted and function ends
  return x;  // The key is returned
}
//-----------------------------------------------------------------------------
// The tail_remove function removes the last nod on the list and
// returns the value stored in that node
//-----------------------------------------------------------------------------
number llist::tail_remove() {
  if(is_empty()) {
    return MISSING;
  }
  node *remover;
  number x;
  
  remover=head->previous;
  x=remover->key;
  //First the next field of the previous node is set to head
  remover->previous->next = head;
  //Then the previous field of head is set to the node before remover
  head->previous = remover->previous;
  size--;            //update size
  delete remover;    // Finally the node can be deleted
  return x; // The key is returned
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
node *llist::set_node(node *set, number x) {
  int error=1;
  node *comparer;
  
  if (set==NULL)
    set=head->next;
  comparer = head->next;
  while(comparer != head) {
    if(comparer == set) {
      error=0;
      break;
    }
    comparer = comparer->next;
  }
  
  if(error==1) {
    return NULL;
  }
  else {
    if(set->key != MISSING){
      set->key = x;
      return set->next;
    }
    else {
      //if the node is MISSING, then don't replace
      //that key.  instead, replace the first non-MISSING
      //node you come to.
      return set_node(set->next,x);
    }
  }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
int llist::is_empty() {
  if(head->next==head)
    return 1;
  else
    return 0;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void copy(llist &L1,const llist &L2) {
  while (!L1.is_empty()) 
    L1.tail_remove();
  
  node *comparer;
  comparer = L2.head->previous;
  while (comparer!=L2.head) {
    L1.insert(comparer->key);
    comparer = comparer->previous;
  }
}
number llist::sumlist(){
  
  number sum = 0;
  
  node* cur;
  cur = head->previous;
  while(cur != head){
    sum += cur->key;
    cur = cur->previous;
  }
  
  return sum;
}

void llist::sumlist(number &prev_sum, int sign){
  
  //printf("in sumlist(number &prev_sum, int sign)\n");
  number sum = 0;
  
  node* cur;
  cur = head->previous;
  while(cur != head){
    sum += cur->key;
    cur = cur->previous;
  }
  
  if(sign*sum > sign * prev_sum)
    prev_sum = sum;
  
  
  return;
}
number llist::maxlist(){
  number max = 0;
  
  node * cur;
  cur = head->previous;
  while(cur != head){
    if(cur->key > max)
      max = cur->key;
    cur = cur->previous;
  }
  
  return max;
}
number llist::minlist(){
  number min = 0;
  
  node * cur;
  cur = head->previous;
  while(cur != head){
    if(cur->key < min)
      min = cur->key;
    cur = cur->previous;
  }
  
  return min;
}
number llist::kthLargest(int k) {
  if(k < 1 || k > size)
    return MISSING;
  else if(k == 1)
    return minlist();
  else if(k == size) 
    return maxlist();
  
  else{
    //place the list in an array for
    //easier selection
    number *A;
    int i; 
    number return_value;
    node* cur = head->previous;
    A = new number[size];
    if(A != NULL){
      for(i = 0; i < size; i++){
        if(cur->key != MISSING)
          A[i] = cur->key;
        cur = cur->previous;
      }
      select(A,0,size-1,k);
      
      return_value = A[k-1];
      delete []A;
    }
    else {
      long_select(k);
    }
    return return_value;
  }
}
number llist::quartile(int q) {
  //q0 is the minimum
  if(q == 0)
    return minlist();
  //q4 is the maximum
  else if(q == 4)
    return maxlist();
  
  //q2 is the median
  else if(q == 2) {
    //if the size of the list is even, there is no exact median
    //so take the average of the two closest numbers
    if(size%2 == 0){
      double t1 = kthLargest(size/2);
      double t2 = kthLargest(size/2 + 1);
      return (t1+t2)/2;
    }
    else
      return kthLargest(1 + (size-1)/2);
  }
  
  //q1 is the first quartile, q3 is the third quartile
  else if(q == 1 || q == 3){
    //if (size+1) is not divisble by 4, there is no exact quartiles
    //so take the weighted average of the two closest numbers
    int k;
    if((k = ((size-1)%4)) != 0){
      int bottom = (int)floor(q*(size-1)/4);
      double t1 = (4-k) * kthLargest(bottom+1);
      double t2 = (k) * kthLargest(bottom+2);
      return (t1+t2)/4;
    }
    else
      return kthLargest(1 + q*(size-1)/4);
  }
  else
    return MISSING;
}
//safe percentile is a safer version of percentile that 
//takes MISSING values into account
number llist::safe_percentile(double percentage) {
  llist temp;
  node* cur = head->next; 
  while(cur != head){ 
    if(cur->key != MISSING) 
      temp.insert(cur->key);
    cur = cur->next; 
  } 
  return temp.percentile(percentage);
}
number llist::percentile(double percentage) {
  int k;
  
  //the argument may not be in correct demical 
  //representation of a percentage, that is,
  //it may be a whole number like 25 instead of .25
  if(percentage > 1)
    percentage = percentage / 100;
  k = (int)(percentage * size);
  return kthLargest(k);
}
number llist::long_select(int k) {
  //haven't gotten around to doing this function yet.
  //it's pretty low priority for me.
  printf("Low Memory.\n");
  return MISSING;
  
}
//-----------------------------------------------------------------------------
//**********   CLOSE OF FUNCTION DEFINITIONS FOR CLASS:  llist        *********
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//These three functions, partition(), select(), and exch() are used to select
//the kth largest number in an array.
//Partition partitions a subarray around the a key such that all entries above
//  the key are greater than the key, and all entries below are less than it.
//  The rightmost element in the subarray is used as the key.
//Select arranges the array in such a way that the kth largest item in the 
//  array is in the kth spot, that is at index # (k-1).
//Exch simply switches the values of the two arguments.
//To possibly speed up the process, these three functions could be combined, 
//  meaning there would be far fewer function calls.
//-----------------------------------------------------------------------------
int partition(number a[], int left, int right) {
  number val = a[right];
  int i = left - 1;
  for(int j = left; j < right; j++){
    if(a[j] <= val){
      i++;
      exch(a[i],a[j]);
    }
  }
  exch(a[i+1],a[right]);
  return i+1;
}

void select(number a[], int l, int r, int k) {
  int i;
  if (r <= l) 
    return;
  
  i = partition(a, l, r);
  if (i > k-1) 
    select(a, l, i-1, k);
  else if (i < k-1) 
    select(a, i+1, r, k);
  else 
    return;
}
void exch(number &x, number &y) {
  number temp;
  temp = x; 
  x = y;
  y = temp;
}
//-----------------------------------------------------------------------------
