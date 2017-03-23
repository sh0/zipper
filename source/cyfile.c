/* module: grid.c grid format database file tools */

static char SccsId[] = "@(#)cyfile.c	1.45";

/* This module is used to hide the format of the Cyberware type image file
 * for historical reasons this file is in a compressed binary format and
 * may have a header which is inaccessable (easily) by some architectures,
 * especially iAPX xxx86 processors.  Also, the files from the digitizers
 * come in the old non-portable binary header style and in a newer portable
 * ASCII header style.  Use of this module will isolate you from all this
 * ugliness.
 *
 * Use these functions as follows:
 *
 *	#include "cyfile.h"
 *	GSPEC *cyread(int fd);
 *	int cywrite(GSPEC *gs, int fd);
 *	int cyfree(GSPEC *gs);
 *	long getr(GSPEC *gs, int latitude, int longitude);
 *	void putr(GSPEC *gs, int latitude, int longitude, long radius);
 *	int cyfree(GSPEC *gs);
 *
 *	Cyread() allocates a new set of buffers each time it is called. If
 *	this is not what you intend, be sure to call cyfree(gs) first.
 *	Opening and closing the file descriptors is up to the caller.
 *	GETR() and PUTR() are inline macros for getr() and putr(), they
 *	usually execute about twice as fast as the function versions.
 *
 * Use the header variables as follows:
 *
 *	char gs->name		Image name string (40 characters max).
 *	long gs->ltincr		Increment between latitudes, y, microns.
 *	long gs->lgincr		Increment between longitudes, ur or um.
 *	short gs->ltmin		Data window, min latitude, inclusive.
 *	short gs->ltmax		Data window, max latitude, inclusive.
 *	short gs->lgmin		Data window, min longitude, inclusive.
 *	short gs->lgmax		Data window, max longitude, inclusive.
 *	short gs->nlt		Total number of latitudes.
 *	short gs->nlg		Total number of longitudes.
 *	long gs->flags		Bit flags, as below.
 *
 *	FLAG_OLDHEADER		Force writing of old style header.
 *	FLAG_CARTESIAN		Indicates cartesian data, lgincr is
 *				in microns (x), radius becomes z. Think
 *				of getr() as  z = getr(gs, yi, xi).
 *	FLAG_BILATERAL		Cartesian and bilateral, ie: nus hands.
 *
 *	Use other header variable at own risk.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HIGHC
#	include <fcntl.h>
#	include <types.h>
#	include <stat.h>
#	include <io.h>
#	include <stdlib.h>
#endif

#include "strings.h"
#include "malloc.h"

#include "cyfile.h"

/* external declarations */

extern int errno;
extern char *ctime();
extern long lseek();
extern void perror();

GSPEC *gsallo();
long getheader();

/*************************** Public Functions ********************************/



/* Cyread optionally allocates buffer space for an image and its header; and
 * optionally reads an image file into these buffers.  If buffers are not
 * yet allocated then call with gs set to NULL.  The return value will be
 * a pointer to the header structure, gs.  If a file is to be read then
 * open a file and pass the descriptor in fd.  If fd is -1 then no files are
 * read and the buffers, if any, have undefined contents.
 */

GSPEC *cyread(gs, fd)	/* read data file from fd */

GSPEC *gs;
int fd;
{
	/* if gs is NULL allocate gs structure, if fd is not -1 read file */

	if (gs == NULL) {
		if ((gs = gsallo()) == NULL) {	/* allocate header memory */
			return(NULL);
		}
	} else {
		if (gs->base != NULL) {
			free((char *)gs->base);
			gs->base = NULL;
		}
	}
	if (fd != -1) {
		if (gsget(gs, fd) == -1) {	/* read header */
			return(NULL);
		}
		if (gs->base == NULL) {		/* not yet allocated ? */
			if (gdallo(gs) == -1) {	/* allocate data memory */
				return(NULL);
			}
		}
		if (gdget(gs, fd) == -1) {	/* read data */
			return(NULL);
		}
	}
	return(gs);
}



/* Cywrite writes the header and image data defined by the header to the
 * file with open descriptor fd.  The header and buffer contents are
 * not altered in any way.  Use cyfree to release the buffers if necessary.
 */

int cywrite(gs, fd)	/* write data file to fd */

GSPEC *gs;
int fd;
{
	if (gsput(gs, fd) == -1) {	/* write header */
		return(-1);
	}
	if (gdput(gs, fd) == -1) {
		return(-1);
	}
	return(0);
}



/* Cyfree will release any memory resources associated with the header gs
 * and its image buffer.
 */

void cyfree(gs)		/* free private resources */

GSPEC *gs;
{
	if (gs != NULL) {
		if (gs->base != NULL) {
			free((char *)gs->base);
		}
		free((char *)gs);
	}
}



/* Getr and putr are to be used to access the image data.  Please note
 * the histroically backwards order of the arguments lt and lg.
 * Getr may return the value CYVOID, which is a very large negative number.
 * The return value should always be tested for CYVOID unless you can be
 * sure that all have been filled.  Putr will accept CYVOID as the radius
 * argument to store a void value.
 */

long getr(gs, lt, lg)		/* a function version of GETR() used globally */

register GSPEC *gs;
register int lt, lg;
{
#ifdef HIGHC
	long offset;
	unsigned short range_packed;
	long range_um;

	offset = lg * gs->nlt + lt;
	range_packed = *(unsigned short *)(gs->base + offset);
	range_packed = range_packed << 8 | range_packed >> 8;
	range_um = (long)(short)range_packed << gs->rshift;
	return range_um;
#else
	return(GETR(gs, lt, lg));
#endif
}



void putr(gs, lt, lg, r)	/* a function version of PUTR() used globally */

register GSPEC *gs;
register int lt, lg, r;
{
#ifdef HIGHC
	long offset;
	unsigned short range_packed;
	offset = lg * gs->nlt + lt;

	range_packed = (unsigned short)(r >> gs->rshift);
	range_packed = range_packed << 8 | range_packed >> 8;
	*(gs->base + offset) = (short)range_packed;
#else
	PUTR(gs, lt, lg, r);
#endif
}



/*************************** Private Functions *******************************/

/* The loops around read()s facilitate the use of pipes, which may not always
 * read an entire header or data array at one time.  HP Integral i/o is
 * a little slow, so a line of dots is written across stderr to keep the user
 * awake.
 */


int gsget(gs, fd)		/* read GSPEC structure from file fd */

GSPEC *gs;
int fd;
{
	unsigned count = sizeof(GSPEC);		/* number of bytes in header */
	char *addr = (char *)gs;			/* start of header */
	int n;
	short *base_save = gs->base;		/* save address of start of data */

	if (lseek(fd, (long)0, 0) == -1L) { /* seek to beginning of file */
		perror(STR026);
		return(-1);
	}

	/* assume the header is of the older binary type */
	while (count > 0) {
		if ((n = read(fd, addr, count)) == -1) { /* n has number bytes read */
			perror(STR027);
			gs->base = base_save;
			return(-1);
		}
		count -= n; /* decrement count by number bytes read */
		addr += n;  /* update ptr to header structure */
	}

	/* determine header type */
	if (gs->offset != 122 && gs->offset != 114 && gs->offset != 128) {
		gs->flags |= FLAG_OLDHEADER;
		if (*((char *)gs + 4) == 'r') {
			/* reread header as portable type */
			if ((gs->offset = getheader(fd)) == -1) {
				puts(STR107);	/* some format problem */
				gs->base = base_save;
				return(-1);
			}
			if (makegsheader(gs) == -1) {
				puts(STR107);	/* some format problem */
				gs->base = base_save;
				return(-1);
			}
		} else {
			puts(STR106);		/* undefined header type */
			gs->base = base_save;
			return(-1);
		}
	}
	gs->base = base_save;
	gs->saved = 0;
	gs->valid = 0;
	return(0);
}



int gsput(gs, fd)		/* write GSPEC structure to file fd */

GSPEC *gs;
int fd;
{
	unsigned count = sizeof(GSPEC);

#ifdef HIGHC
	if (lseek(fd, (long)0, 0) == -1) {
		perror(STR026);
		return(-1);
	}
#else
	if (ftruncate(fd, (long)0) == -1 || lseek(fd, (long)0, 0) == -1) {
		perror(STR026);
		return(-1);
	}
#endif

	if (gs->flags & FLAG_OLDHEADER) {
		gs->offset = count;
		if (write(fd, (char *)gs, count) != count) {
			perror(STR028);
			return(-1);
		}
	} else {
		if (writegsheader(gs, fd) == -1) {
			return(-1);
		}
	}
	return(0);
}



int gdget(gs, fd)		/* read data from grid file fd */

GSPEC *gs;
int fd;
{
	unsigned long count = (long)sizeof(short) * (long)gs->nlt * (long)gs->nlg;
	unsigned int n;
	unsigned int readsize;
	char *addr;
	unsigned size = count;

	/* if unallocated, allocate image memory */
	if (gs->base == NULL) {
		if (gdallo(gs) == -1) {
			return(-1);
		}
	}

	if (lseek(fd, gs->offset, 0) == -1L) {
		perror(STR026);
		return(-1);
	}
	addr = (char *)gs->base;
	while (count > 0) {
		readsize = (unsigned int) MIN(size, count);
		if ((n = read(fd, addr, readsize)) == -1L) {
			perror(STR027);
			return(-1);
		}
		count -= (unsigned long)n;
		addr += n;
	}
	return(0);
}



int gdput(gs, fd)		/* write data to grid file */

GSPEC *gs;
int fd;
{
	unsigned long count;
	unsigned int n;
	unsigned int writesize;
	char *addr = (char *)gs->base;
	unsigned size = 64*1024;		/* a large block size */

	count = (long)sizeof(short) * (long)gs->nlt * (long)gs->nlg;

	if (gs->flags & FLAG_OLDHEADER) {
		if (lseek(fd, (long)gs->offset, 0) == (long)-1) {
			perror(STR026);
			return(-1);
		}
	}
	while (count > 0) {
		writesize = (unsigned int) MIN(size, count);
#ifdef HIGHC
		if ((n = _write(fd, addr, writesize)) == -1L) {  /* ??? really?? */
#else
		if ((n = write(fd, addr, writesize)) == -1L) {
#endif
			perror(STR028);
			return(-1);
		}
		count -= (unsigned long)n;
		addr += n;
	}
	return(0);
}



int gdallo(gs)				/* allocate a data buffer */

GSPEC *gs;
{
	unsigned long size;

	size = (unsigned long)gs->nlt *
					(unsigned long)gs->nlg * (unsigned long)sizeof(short);
	gs->base = (short *)malloc(size);

	if (gs->base == NULL) {
		puts(STR082);
		return(-1);
	} else {
		return(0);
	}
}



GSPEC *gsallo()				/* allocate a GSPEC structure */
{
	GSPEC *gs;

	gs = (GSPEC *)malloc((unsigned)sizeof(GSPEC));
	if (gs == NULL) {
		puts(STR082);
		return(NULL);
	}
	gs->base = NULL;
	return(gs);
}



#define MAXHEADER 4096		/* ??? might hang on very short files */
#define HEADEREND "DATA=\n"

static char *header = 0;

long getheader(fd)	/* get header and seek to data */

int fd;
{
	int count;
	char *end;
	char *h;
	char *endstr = HEADEREND;
	char *temp_header;
	char *addr;
	int n;

	temp_header = malloc(MAXHEADER);

	if (lseek(fd, (long)0, 0) == -1) {
		perror(STR108);
		return(-1);
	}
	addr = temp_header;
	for (count = 0; count < MAXHEADER; count += n) {
		if ((n = read(fd, addr, (unsigned)MAXHEADER)) == -1) {
			perror(STR109);
			return(-1);
		}
		addr += n;
	}

	/* end of header is eof or endstr string */
	end = temp_header + count;
	for (h = temp_header; h < end; ++h) {
		if (*h == endstr[0]) {
			if (strncmp(endstr, h, strlen(endstr)) == 0) {
				end = h + strlen(endstr);
				break;
			}
		}
	}
	count = end - temp_header;
	if (header != 0) {
		free(header);
	}
	header = malloc((unsigned)(count+1));
	strncpy(header, temp_header, count);
	header[count] = 0;			/* null terminate */
	free(temp_header);
	return(count);
}



int getvalue(name, dest, length)  /* copy value of name to dest or return -1 */

char *name;
char *dest;
int length;
{
	char *h = header;
	int n;
	char *p;

	if (header == 0) {								/* no header, oops! */
		puts("getvalue: no header");
		exit(-1);		/* fatal coding error */
	}
	n = strlen(name);
	while ((h = strchr(h, '\n')) != 0) {	/* move to next newline */
		h += 1;										/* skip over newline */
		if (strncmp(h, name, n) == 0) {	/* compare names */
			h += strlen(name);	/* skip over matched name */
			if (*h == '=') {	/* verify assignment char */
				h += 1;
				/* no value terminator ? */
				if ((p = strchr(h, '\n')) == 0) {
					puts(STR110);
					return(-1);
				}
				*p = 0;		/* temporary termination */
				strncpy(dest, h, length);
				*p = '\n';	/* restore terminator */
				return(0);
			}
		}
	}
	return(-1);				/* no match */
}



#define STRINGLEN	24

int makegsheader(gs)	/* fill GSPEC structure from portable header */

GSPEC *gs;
{
	char string[STRINGLEN+1];
	long i;

	string[STRINGLEN] = 0;

	/* defaults */
	gs->flags = 0;

	/* mandatory items */
	if (getvalue("NLT", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "NLT");
		return(-1);
	}
	gs->nlt = atoi(string);
	if (getvalue("NLG", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "NLG");
		return(-1);
	}
	gs->nlg = atoi(string);
	if (getvalue("LGSHIFT", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "LGSHIFT");
		return(-1);
	}
	gs->lgshift = atoi(string);
	if (getvalue("LTINCR", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "LTINCR");
		return(-1);
	}
	gs->ltincr = atol(string);
	if (getvalue("LGINCR", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "LGINCR");
		return(-1);
	}
	gs->lgincr = atol(string);
	if (getvalue("RSHIFT", string, STRINGLEN) == -1) {
		printf("%s: %s\n", STR111, "RSHIFT");
		return(-1);
	}
	gs->rshift = atoi(string);

	/* optional items */
	if (getvalue("NAME", gs->name, NAMELEN) == -1) {
		for (i = NAMELEN-1; i >= 0; --i) gs->name[i] = 0;
	}
	if (getvalue("LTMIN", string, STRINGLEN) == -1) {
		gs->ltmin = 0;
	} else {
		gs->ltmin = atoi(string);
	}
	if (getvalue("LTMAX", string, STRINGLEN) == -1) {
		gs->ltmax = gs->nlt - 1;
	} else {
		gs->ltmax = atoi(string);
	}
	if (getvalue("LGMIN", string, STRINGLEN) == -1) {
		gs->lgmin = 0;
	} else {
		gs->lgmin = atoi(string);
	}
	if (getvalue("LGMAX", string, STRINGLEN) == -1) {
		gs->lgmin = gs->nlg - 1;
	} else {
		gs->lgmax = atoi(string);
	}
	if (getvalue("RMIN", string, STRINGLEN) == -1) {
		gs->rmin = 0;
	} else {
		gs->rmin = atol(string);
	}
	if (getvalue("RMAX", string, STRINGLEN) == -1) {
		gs->rmax = 0;
	} else {
		gs->rmax = atol(string);
	}
	if (getvalue("SCALE", string, STRINGLEN) == -1) {
		gs->scale = 100.0;
	} else {
		gs->scale = atof(string);
	}
	if (getvalue("RPROP", string, STRINGLEN) == -1) {
		gs->rprop = 100.0;
	} else {
		gs->rprop = atof(string);
	}
	if (getvalue("FILLED", string, STRINGLEN) == -1) {
		gs->filled = 0;
	} else {
		gs->filled = 1;
	}
	if (getvalue("SMOOTHED", string, STRINGLEN) == -1) {
		gs->smoothed = 0;
	} else {
		gs->smoothed = 1;
	}
	if (getvalue("SPACE", string, STRINGLEN) == -1) {
		gs->flags = 0;
	} else {
		if (strcmp(string, "CARTESIAN") == 0) {
			gs->flags |= FLAG_CARTESIAN;
		} else if (strcmp(string, "CYLINDRICAL") == 0) {
			gs->flags &= ~FLAG_CARTESIAN;
		} else if (strcmp(string, "BILATERAL") == 0) {
			gs->flags |= FLAG_CARTESIAN;
			gs->flags |= FLAG_BILATERAL;
		} else {
			printf("%s: SPACE\n", STR112);
			return(-1);
		}
	}
	if (getvalue("INSIDE_OUT", string, STRINGLEN) != -1) {
		gs->flags |= FLAG_INSIDE_OUT;
	}
	if (getvalue("COLOR", string, STRINGLEN) != -1) {
		gs->flags |= FLAG_COLOR;
	}
	if (getvalue("THETA_RIGHTHAND", string, STRINGLEN) != -1) {
		gs->flags |= FLAG_THETARIGHT;
	}

	/* forced value items */
	gs->time = 0;
	gs->camera = 0;
	gs->setup = 0;
	gs->saved = 0;
	gs->valid = 0;
	gs->ltsize = gs->nlt * gs->ltincr;
	gs->lgsize = gs->nlg * gs->lgincr;
	return(0);
}



int writegsheader(gs, fd)	/* write a portable header from GSPEC struct */

GSPEC *gs;
int fd;
{
	/* Uses \012 instead of \n for the sake of DOS which outputs a
	 * CR-LF sequence with \n !!! */

	/* Open a stream on a duplicate file descriptor.  The duplicate is
	 * necessary to allow the stream to be closed with out closing the
	 * callers file descriptor */

	FILE *fp = fdopen(dup(fd), "w");

	/* error status only checked on first line and last line */
	if (fprintf(fp, "Cyberware Digitizer Data\012") < 0) {
		perror(STR028);
		fclose(fp);
		return(-1);
	}
	fprintf(fp, "NAME=%.40s\012", gs->name);
	fprintf(fp, "DATE=%s", ctime(&gs->time));
	if (gs->flags & FLAG_CARTESIAN) {
		if (gs->flags & FLAG_BILATERAL) {
			fprintf(fp, "SPACE=BILATERAL\012");
		} else {
			fprintf(fp, "SPACE=CARTESIAN\012");
		}
	} else {
		fprintf(fp, "SPACE=CYLINDRICAL\012");
	}
	if (gs->flags & FLAG_COLOR) {
		fprintf(fp, "COLOR=SGI\012");
	}
	if (gs->flags & FLAG_INSIDE_OUT) {
		fprintf(fp, "INSIDE_OUT=TRUE\012");
	}
	if (gs->flags & FLAG_THETARIGHT) {
		fprintf(fp, "THETA_RIGHTHAND=TRUE\012");
	}
	fprintf(fp, "NLG=%-.1d\012", gs->nlg);
	fprintf(fp, "LGINCR=%-.1ld\012", gs->lgincr);
	fprintf(fp, "LGMIN=%-.1d\012", gs->lgmin);
	fprintf(fp, "LGMAX=%-.1d\012", gs->lgmax);
	fprintf(fp, "NLT=%-.1d\012", gs->nlt);
	fprintf(fp, "LTINCR=%-.1ld\012", gs->ltincr);
	fprintf(fp, "LTMIN=%-.1d\012", gs->ltmin);
	fprintf(fp, "LTMAX=%-.1d\012", gs->ltmax);
	fprintf(fp, "RMIN=%-.1ld\012", gs->rmin);
	fprintf(fp, "RMAX=%-.1ld\012", gs->rmax);
	fprintf(fp, "RSHIFT=%-.1d\012", gs->rshift);
	fprintf(fp, "LGSHIFT=%-.1d\012", gs->lgshift);
	fprintf(fp, "SCALE=%-.2f\012", (float)gs->scale);
	fprintf(fp, "RPROP=%-.2f\012", (float)gs->rprop);
	if (gs->filled) {
		fprintf(fp, "FILLED=1\012");
	}
	if (gs->smoothed) {
		fprintf(fp, "SMOOTHED=1\012");
	}
	fprintf(fp, "DATA=\012");
	if (fflush(fp) == EOF) {
		perror(STR028);
		fclose(fp);
		return(-1);
	}
	fclose(fp);
	return(0);
}
