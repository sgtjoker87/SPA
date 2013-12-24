/*

*/



//===========================================================================
#ifndef FMU_HXX
#define FMU_HXX

#ifndef PI
#define  PI       3.14159265358979323846    /* Circle circumference/Diameter */
#endif

#define  PIHALF   1.57079632679489661923    /* pi/2                          */
#define  PITHRD   1.04719755119659774615    /* pi/3                          */
#define  PI4TH    0.785398163397448309616   /* pi/4                          */

#ifndef TWOPI
#define  TWOPI    6.28318530717958647692    /* 2*pi                          */
#endif

#define  THREPI   9.42477796076937971539    /* 3*pi                          */

#define  PISQ     9.86960440108935861883    /* pi*pi                         */
#define  SQRTPI   1.772453850905516027298   /* sqrt(pi)                      */
#define  DEGREE  57.2957795130823208767     /* # of degrees/radian           */
#define  RADIAN   0.0174532925199432957692  /* # of radians/degree           */
#define  SQRT2    1.4142135623730950488     /* sqrt(2)                       */
#define  SQRT3    1.7320508075688772935     /* sqrt(3)                       */
#define  BASEE    2.71828182845904523536    /* e                             */
#define  BIG      1.e25                     /* infinity                      */
#define  SMALL    1.e-25                    /* infinitesimal                 */
#define  MAXSIGDIG 1.e-14                   /* 14 significant digits         */

struct fmu_TzContourLoop
{
    int     nNodes;
    int     iLoopDir;
    int     *aiNodes;
    double  *adCoordX;
    double  *adCoordY;
    double  *adG;
};

#endif
