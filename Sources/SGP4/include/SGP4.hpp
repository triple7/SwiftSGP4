#ifndef _SGP4h_
#define _SGP4h_

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include<stdbool.h>

/*
These are the function declarations for the SGP4 satellite orbit model. The SGP4 library is commonly used for predicting the positions of satellites based on orbital elements. The functions provide various utilities for initializing the satellite's orbit, propagating the orbit over time, converting coordinates, and performing calculations related to orbital elements and time.
The declarations include functions for initializing the satellite, propagating the orbit, converting between coordinate systems (e.g., Earth-centered Earth-fixed (ECEF) and geodetic), calculating satellite visibility (azimuth and elevation), and various mathematical operations.
*/

#define pi 3.14159265358979323846
#define SGP4Version  "SGP4 Version 2023-06-16"

typedef enum
{
  wgs72old,
  wgs72,
  wgs84
} gravconsttype;

typedef struct elsetrec
{
  char satnum[6];
  int epochyr, epochtynumrev;
  int error;
  char operationmode;
  char init, method;

  /* Near Earth */
  int isimp;
  double aycof  , con41  , cc1    , cc4      , cc5    , d2      , d3   , d4    ,
         delmo  , eta    , argpdot, omgcof   , sinmao , t       , t2cof, t3cof ,
         t4cof  , t5cof  , x1mth2 , x7thm1   , mdot   , nodedot, xlcof , xmcof ,
         nodecf;

  /* Deep Space */
  int    irez;
  double d2201  , d2211  , d3210  , d3222    , d4410  , d4422   , d5220 , d5232 ,
         d5421  , d5433  , dedt   , del1     , del2   , del3    , didt  , dmdt  ,
         dnodt  , domdt  , e3     , ee2      , peo    , pgho    , pho   , pinco ,
         plo    , se2    , se3    , sgh2     , sgh3   , sgh4    , sh2   , sh3   ,
         si2    , si3    , sl2    , sl3      , sl4    , gsto    , xfact , xgh2  ,
         xgh3   , xgh4   , xh2    , xh3      , xi2    , xi3     , xl2   , xl3   ,
         xl4    , xlamo  , zmol   , zmos     , atime  , xli     , xni;

  double a, altp, alta, epochdays, jdsatepoch, jdsatepochF, nddot, ndot,
     bstar, rcse, inclo, nodeo, ecco, argpo, mo, no_kozai;
  char  classification, intldesg[11];
  int   ephtype;
  long  elnum    , revnum;
  double no_unkozai;
  double am     , em     , im     , Om       , om     , mm      , nm;
  double tumin, mus, radiusearthkm, xke, j2, j3, j4, j3oj2;

  long dia_mm; // RSO dia in mm
  double period_sec; // Period in seconds
  unsigned char active; // "Active S/C" flag (0=n, 1=y) 
  unsigned char not_orbital; // "Orbiting S/C" flag (0=n, 1=y)  
  double rcs_m2; // "RCS (m^2)" storage  
} elsetrec;

bool sgp4init(gravconsttype whichconst,//*    whichconst  - which set of constants to use  72, 84
              char opsmode, //*    opsmode     - mode of operation afspc or improved 'a', 'i'
              const char satn[9], //*    satn        - satellite number
              const double epoch, //*    epoch       - epoch time in days from jan 0, 1950. 0 hr
              const double xbstar, //*    bstar       - sgp4 type drag coefficient              kg/m2er
              const double xndot, // NO IDEA
              const double xnddot,  // NO IDEA
              const double xecco, //*    ecco        - eccentricity
              const double xargpo, //*    argpo       - argument of perigee (output if ds)
              const double xinclo, //*    inclo       - inclination
              const double xmo, //*    mo          - mean anomaly (output if ds)
              const double xno_kozai, //*    no          - mean motion
              const double xnodeo, //*    nodeo       - right ascension of ascending node
              elsetrec* satrec);

/*
 satrec.no_kozai = satrec.no_kozai / xpdotp; // rad/min
 satrec.ndot = satrec.ndot / (xpdotp*1440.0);  // ? * minperday
 satrec.nddot = satrec.nddot / (xpdotp*1440.0 * 1440);
 satrec.inclo = satrec.inclo  * deg2rad;
 satrec.nodeo = satrec.nodeo  * deg2rad;
 satrec.argpo = satrec.argpo  * deg2rad;
 satrec.mo = satrec.mo     * deg2rad;


satrec.ndot = satrec.ndot / (xpdotp*1440.0);  // ? * minperday
satrec.nddot = satrec.nddot / (xpdotp*1440.0 * 1440);



_ = sgp4init(wgs72,
             opsMode,
             &genSatNum,
             jdEpoch,
             target.BSTAR,
             target.MEAN_MOTION_DOT/xpdotInv, // looks correct
             target.MEAN_MOTION_DDOT/xpdotInv2, // looks correct
             target.ECCENTRICITY*deg2rad,
             target.ARG_OF_PERICENTER*deg2rad,
             target.INCLINATION*deg2rad,
             target.MEAN_ANOMALY*deg2rad,
             target.MEAN_MOTION/xpdotp,
             target.RA_OF_ASC_NODE*deg2rad,
             &satrec)


*/





bool sgp4(elsetrec* satrec, double tsince,double r[3], double v[3]);

double  gstime(double jdut1        );

void getgravconst(gravconsttype whichconst,double* tumin,double* mus,double* radiusearthkm,double* xke,double* j2,double* j3,double* j4,double* j3oj2);

// older sgp4io methods
void twoline2rv(char      longstr1[130], char longstr2[130],char      typerun, char typeinput, char opsmode,gravconsttype       whichconst,double* startmfe, double* stopmfe, double* deltamin,elsetrec* satrec);

double asinh(double xval);

double sgn(double x);

double  mag(double x[3]);

void    cross(double vec1[3], double vec2[3], double outvec[3]);

double  dot(double x[3], double y[3]);

double  angle(double vec1[3], double vec2[3]);

void    newtonnu(double ecc, double nu, double* e0, double* m);

void    rv2coe(double r[3], double v[3], double mu, double* p, double* a, double* ecc, double* incl, double* omega, double* argp,          double* nu, double* m, double* arglat, double* truelon, double* lonper);

void jday(int year, int mon, int day, int hr, int minute, double sec,double* jd, double* jdFrac);

void days2mdhms(int year, double days, int* mon, int* day, int* hr, int* minute, double* sec);

void    invjday(double jd, int timezone, bool daylightsaving, int* year, int* mon, int* day,
          int* hour, int* minute, double* sec);

double floatmod(double a, double b);

bool summertime(int year, int month, int day, int hour, int tzHours);

void teme2ecef(double rteme[3], double jdut1, double recef[3]);

void teme2ecefOptimised(double rteme[3], double jdut1, double gmstCos, double gmstSin, double recef[3]);

void polarm(double jdut1, double pm[3][3]);

void ijk2ll(double r[3], double latlongh[3]);

void site(double latgd, double lon, double alt, double rs[3]);

void rv2azel(double ro[3], double latgd, double lon, double alt, double jdut1, double razel[3]);

void rot3(double invec[3], double xval, double outvec[3]);

void rot2(double invec[3], double xval, double outvec[3]);

double getJulianFromUnix(double unixSecs);

unsigned long getUnixFromJulian(double julian);

#ifdef __cplusplus
}
#endif

#endif
