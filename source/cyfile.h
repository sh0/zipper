/* module:  cyfile.h    echo image header file */

/* @(#)cyfile.h 1.30 */

/* globals */

/* Internal types, These modules all assume the following types:
 *
 *  char            1 byte signed integer, -128...127
 *  unsigned char   1 byte unsigned integer, 0...255
 *  short           2 byte signed integer, -32,768...32,767
 *  unsigned short  2 byte unsigned integer, 0...65,535
 *  long            4 byte signed integer, -2,147,483,648...2,147,483,647
 *  unsigned long   4 byte unsigned integer, 0...4,294,967,295
 *  real            a real variable natural to the machine
 *  int             at least as long as short
 *  unsigned int    at least as long as unsigned short
 *
 *  All other types are to be enclosed in #ifdefs.
 */

/* file constants, unpacked */

#define MAXR        (0x00007fff<<gs->rshift)
#define MAXRGS(gs)  (0x00007fff<<(gs)->rshift)
#define MINR        0
#define CYVOID      (0xffff8000<<gs->rshift)
#define CYVOIDGS(gs)    (0xffff8000<<(gs)->rshift)

/* various constants of pi */

#define RPI     3.141592654     /* floating point pi */
#define PI      3141593         /* pi in urads */
#define TWOPI   6283185         /* 2pi in urads */

#ifndef NULL
#define NULL    0               /* null address */
#endif

/* math tools */

#define MAX(a,b)    ((a)>(b)?(a):(b))       /* return greater of a and b */
#define MIN(a,b)    ((a)<(b)?(a):(b))       /* return lesser of a and b */
#define ABS(i)      ((i)<0?-(i):(i))        /* integer absolute value */
#define DELTA(a,b)  (ABS((a)-(b)))          /* int absolute difference */
#define SCALE(n,s)  ((((n)*(s))+50)/100)    /* int scale n by s percent */
#define WRAPP(n,m)  (if((n)>=(m))(n)-=(m))  /* modulo positive wrap */
#define WRAPN(n,m)  (if((n)<0)(n)+=(m))     /* modulo positive wrap */

/* unit conversions */

#define UMTOI(um)   ((real)(um)*3.937e-5)   /* microns (um) to float inch */
#define ITOUM(um)   ((int)((um)*2.54e4))    /* inches to int microns */
#define URTOD(ur)   ((real)(ur)*5.7296e-5)  /* urads to float degrees */
#define DTOUR(deg)  ((int)((deg)*1.74533e4) /* degrees to int urads */
#define DTOR(deg)   ((deg)*1.7453292e-2)    /* degrees to float radians */
#define RTOD(rad)   ((rad)*57.295779)       /* radians to float degrees */
#define URTOR(ur)   ((real)(ur)*1.e-6)      /* radians to urads */
#define RTOUR(ur)   (int)((ur)*1.e6)        /* radians to urads */

/* this structure defines 'grid file format'.  the file consists of
 * a parameter table followed immediatly by the data table.  the offset
 * to the start of the data table is the second parameter and is therefore
 * fifth thru eighth bytes of the file (msb first).
 *
 * the parameters nlg and nlt are important for accessing the data.  nlg
 * is the number of longitude entries in the table.  nlt is the number of
 * latitudes in the table.  nlt * nlg * 2 gives the number of bytes in the
 * table.
 *
 * the table is a set of radius values in a cylindrical coordinate space.
 * each radius value is stored in a 2 byte integer which when shifted
 * left by RSHIFT bits yields a radius in microns (4 byte long integer).
 * the radius values are stored in longitudnal groups of nlt values.  there
 * are nlg of these groups, one for each longitude of the cylinder.
 *
 * the functions GETR() and PUTR() defined below are usually all that is
 * required to fetch and store values in the table when it is in memory.
 * the parameters ltincr and lgincr define the distance between adjacent
 * latitudes (microns) and adjacent longitudes (microradians) respectively.
 *
 * There are two formats for this header, one portable, one not so
 * portable.  The older non-portable type is binary and has the value
 * 122 decimal ('z') in the fifth byte.  The portable header has a 'r'
 * in the fifth byte.  The portable header is in ascii and has the form
 * [name=value],... where name is a defined ascii symbol and value is a
 * string value for the symbol.  Format is variable and assignments are
 * separated by white space or commas.
 *
 * See header.c for details.
 */

#define NAMELEN     40
#define CREATE_MODE 0644        /* create image files with this mode */

typedef struct {

    /* internal private variables */
    short* base;            /* base of data buffer */
    long offset;                /* file offset to start of data, bytes */

    /* file parameters */
    char name[NAMELEN];         /* subject name */
    long time;                  /* original creation time */
    short camera;               /* camera id number */
    short setup;                /* camera setup code */
    char saved;                 /* file has been saved since modified */
    char valid;                 /* file buffer is valid */

    /* data parameters */
    short nlt;                  /* number of latitude intervals */
    short nlg;                  /* number of longitude intervals */
    short rshift;               /* shift to compress/expand radius data */
    short lgshift;              /* shift to extract longitude from addr */
    long flags;                 /* misc file state flags, see below */
    long ltincr;                /* distance between latitudes, um */
    long lgincr;                /* distance between longitudes, urad */
    long ltsize;                /* nlat * ltincr, um */
    long lgsize;                /* nlg * lgincr, urad (always 2pi in urads) */

    /* user parameters */
    char filled;                /* fill flag, useless */
    short smoothed;             /* smooth pass counter */
    short ltmin, ltmax;         /* latitude window limits, inclusive */
    short lgmin, lgmax;         /* longitude window limits, inclusive */
    long rmin, rmax;            /* radius range, from last run of rminmax */
    double scale;           /* current scale */
    double rprop;           /* current radius proportion */
} GSPEC;

/* macros for standardizing the use of the grid data. gs is a pointer to the
 * applicable GSSPEC table.  index is the offset of a data item in the
 * data. lt and lg are latitude and longitude indicies. r is the radius
 * in microns (um) of a data point. z is a position along the cylindrical
 * axis in microns. a is an angular coordinate around the cylinder in
 * microradians (urad).
 *
 * INDEX generates an index value from latitude and logitude indicies.
 * ADDR returns the absolute address of a data item.
 * PUTR and GETR are used to store and retrieve data from the image.
 */

#define INDEX(gs, lt, lg)   ((lg) * (gs)->nlt + (lt))
#define ADDR(gs, lt, lg)    ((gs)->base + INDEX(gs, lt, lg))

#define GETR(gs, lt, lg)          getr(gs,lt,lg)
#define PUTR(gs, lt, lg, r)       putr(gs,lt,lg,r)

/* flag bits for gs->flags */

#define FLAG_RESERVED   0x000000ff  /* older files have ones here, ignore */
#define FLAG_CARTESIAN  0x00000100  /* data is cartesian (vs. cyl) */
#define FLAG_OLDHEADER  0x00000200  /* please write file with old header */
#define FLAG_BILATERAL  0x00000400  /* bilateral image, ie: nus hands */
#define FLAG_COLOR      0x00000800  /* image has associated color file */
#define FLAG_THETARIGHT 0x00001000  /* theta is right hand rule */
#define FLAG_INSIDE_OUT 0x00002000  /* inside surface is outside */

/* non-int public functions */

extern GSPEC* cyread();
extern GSPEC* gsallo();
extern int cywrite();
extern void cyfree();
extern long getr();
extern void putr();
extern int geget();
extern int gsget();
extern int gdget();
extern int gdput();
extern int gsput();
extern int gdallo();
extern long getheader();
extern int getvalue();
extern int makegsheader();
extern int writegsheader();
