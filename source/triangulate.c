/*

Use a greedy algorithm to split a polygon into triangles.

Greg Turk

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.

Permission to use, copy, modify and distribute this software and its
documentation for any purpose is hereby granted without fee, provided
that the above copyright notice and this permission notice appear in
all copies of this software and that you do not sell the software.

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "triangulate.h"

/* screen stuff */

typedef struct Edge {
    int p1, p2;
    float len;
    float a, b, c;
    int final;
} Edge;

typedef struct Point {
    Vector pos;           /* location in plane */
    Vector pos3d;         /* old location in 3-space */
    int boundary;
    struct Edge** edges;
    int nedges;
    int max_edges;
    int index;
} Point;

/* the set of points to be triangulated */

static Point** points = NULL;
static int npoints = 0;
static int max_points = 20;



/* list of all possible edges */

static Edge** edges = NULL;
static int nedges = 0;
static int max_edges = 400;

/* list of edges that are in the final triangulation */

static Edge** final = NULL;
static int nfinal;

/* list of triangles */

typedef struct Triangle {
    int p1, p2, p3;
    Edge* e1, *e2, *e3;
} Triangle;

static Triangle* tris = NULL;
static int ntris;
static int tri_goal;

static int boundary_count = 0;  /* number of points that form the polygon */

static int rescale_flag = 1;    /* whether to re-scale polygon to fit screen */
static int shuffle_flag = 0;    /* whether to shuffle the edge list */
static int parallel_flag = 1;   /* print parallel edge warning? */

static int x_screen = 500;
static int y_screen = 500;

static Matrix trans_mat, trans_mat_inv;


/******************************************************************************
Initialize the polygon splitter.

Entry:
  a,b,c,d - coefficients of plane equation for plane on which to triangulate

Exit
  returns 1 if we were given a bad plane equation, 0 otherwise
******************************************************************************/

init_splitter(a, b, c, d)
float a, b, c, d;
{
    int i;
    int result;

    /* allocate the edge list or free up old edges */
    if (edges == NULL) {
        edges = (Edge**) malloc(sizeof(Edge*) * max_edges);
    } else {
        for (i = 0; i < nedges; i++)
            free(edges[i]);
    }

    /* allocate the point list or free up old points */
    if (points == NULL) {
        points = (Point**) malloc(sizeof(Point*) * max_points);
    } else {
        for (i = 0; i < npoints; i++) {
            free(points[i]->edges);
            free(points[i]);
        }
    }

    npoints = 0;
    nedges = 0;
    ntris = 0;
    boundary_count = 0;

    if (final != NULL) {
        free(final);
        final = NULL;
    }

    if (tris != NULL) {
        free(tris);
        tris = NULL;
    }

    /* create transformation matrix for mapping vertices onto the */
    /* plane on which the triangulation will take place */

    result = face_to_xy_plane(a, b, c, d, trans_mat, trans_mat_inv);

    return (result);
}


/******************************************************************************
Add a point that is a part of the polygonal boundary.

Entry:
  xx,yy,zz - 3-space position of point
  index    - data to save for calling program
******************************************************************************/

int add_boundary_point(xx, yy, zz, index)
float xx, yy, zz;
int index;
{
    Vector vec, tvec;
    Point* pt;

    if (npoints >= max_points) {
        max_points *= 2;
        points = (Point**) realloc(points, sizeof(Point*) * max_points);
    }

    /* transform (xx,yy,zz) onto plane */
    vec[X] = xx;
    vec[Y] = yy;
    vec[Z] = zz;
    vapply(trans_mat, vec, tvec);

    pt = (Point*) malloc(sizeof(Point));
    pt->pos[X] = tvec[X];
    pt->pos[Y] = tvec[Y];
    pt->pos[Z] = 0;
    pt->pos3d[X] = xx;
    pt->pos3d[Y] = yy;
    pt->pos3d[Z] = zz;
    pt->nedges = 0;
    pt->max_edges = 6;
    pt->edges = (Edge**) malloc(sizeof(Edge*) * pt->max_edges);
    pt->boundary = 1;
    pt->index = index;

    points[npoints] = pt;
    npoints++;

    boundary_count++;
}


/******************************************************************************
Add a point to the list of points interior to polygon.

Entry:
  xx,yy,zz - 3-space position of point
  index    - data to save for calling program

Exit:
  returns 0 if point is inside polygon, 1 if it is outside (bad)
******************************************************************************/

add_point(xx, yy, zz, index)
float xx, yy, zz;
int index;
{
    Vector vec, tvec;
    Point* pt;

    if (npoints >= max_points) {
        max_points *= 2;
        points = (Point**) realloc(points, sizeof(Point*) * max_points);
    }

    /* transform (xx,yy,zz) onto plane */
    vec[X] = xx;
    vec[Y] = yy;
    vec[Z] = zz;
    vapply(trans_mat, vec, tvec);

    /* check to see that point is inside polygon */

    if (point_in_split_poly(tvec[X], tvec[Y]) == 0) {
        return (1);
    }

    /* add the point */

    pt = (Point*) malloc(sizeof(Point));
    pt->pos[X] = tvec[X];
    pt->pos[Y] = tvec[Y];
    pt->pos[Z] = 0;
    pt->pos3d[X] = xx;
    pt->pos3d[Y] = yy;
    pt->pos3d[Z] = zz;
    pt->nedges = 0;
    pt->max_edges = 6;
    pt->edges = (Edge**) malloc(sizeof(Edge*) * pt->max_edges);
    pt->boundary = 0;
    pt->index = index;

    points[npoints] = pt;
    npoints++;

    return (0);  /* point is okay */
}


/******************************************************************************
Return the number of points in the polygon being split.
******************************************************************************/

int split_npoints()
{
    return (npoints);
}


/******************************************************************************
See if point is close to a given value.
******************************************************************************/

int close_orig(p, x, y, z)
Point* p;
float x, y, z;
{
    float len;
    float dx, dy, dz;

    dx = p->pos3d[X] - x;
    dy = p->pos3d[Y] - y;
    dz = p->pos3d[Z] - z;
    len = sqrt(dx * dx + dy * dy + dz * dz);
    if (len < 0.00001)
        return (1);
    else
        return (0);
}


/******************************************************************************
Set the value of the parallel edge flag.
******************************************************************************/

set_parallel_flag(val)
int val;
{
    parallel_flag = val;
}


/******************************************************************************
Set the value of the shuffle flag.
******************************************************************************/

set_shuffle_flag(val)
int val;
{
    shuffle_flag = val;
}


/******************************************************************************
Set the value of the rescale flag.

Entry:
  val - value of rescale flag
  x,y - size of screen
******************************************************************************/

set_rescale_flag(val, x, y)
int val;
int x, y;
{
    rescale_flag = val;
    x_screen = x;
    y_screen = y;
}


/******************************************************************************
Add an edge to the list of edges.
******************************************************************************/

add_edge(i, j)
int i, j;
{
    float dx, dy, dz;
    Edge* edge;

    if (nedges >= max_edges) {
        max_edges *= 2;
        edges = (Edge**) realloc(edges, sizeof(Edge*) * max_edges);
    }

    /* make the new edge */

    edge = (Edge*) malloc(sizeof(Edge));
    edge->p1 = i;
    edge->p2 = j;
    edge->final = 0;

    /* compute length of edge based on original 3-space position */
    dx = points[i]->pos3d[X] - points[j]->pos3d[X];
    dy = points[i]->pos3d[Y] - points[j]->pos3d[Y];
    dz = points[i]->pos3d[Z] - points[j]->pos3d[Z];
    edge->len = sqrt(dx * dx + dy * dy + dz * dz);

    compute_line(points[i]->pos[X], points[i]->pos[Y],
                 points[j]->pos[X], points[j]->pos[Y],
                 &edge->a, &edge->b, &edge->c);

    /* add to list of edges */
    edges[nedges] = edge;
    nedges++;
}


/******************************************************************************
Add an edge to the final edge list.
******************************************************************************/

add_final_edge(e)
Edge* e;
{
    Point* p1, *p2;

    /* add the edge to the final edge list */

    final[nfinal] = e;
    nfinal++;
    e->final = 1;  /* mark this edge as a final edge */

    /* have the edge's points refer to this edge */

    p1 = points[e->p1];
    p2 = points[e->p2];

    /* make sure there is enough room for the new edges */

    if (p1->nedges >= p1->max_edges) {
        p1->max_edges += 4;
        p1->edges = (Edge**) realloc(p1->edges, sizeof(Edge*) * p1->max_edges);
    }

    if (p2->nedges >= p2->max_edges) {
        p2->max_edges += 4;
        p2->edges = (Edge**) realloc(p2->edges, sizeof(Edge*) * p2->max_edges);
    }

    /* add the edges to the points' lists of edges */

    p1->edges[p1->nedges++] = e;
    p2->edges[p2->nedges++] = e;

    /* maybe draw the edge */

#if 0
    if (drawing_flag)
        draw_edge(e, color2);
#endif
}


/******************************************************************************
Copy a bunch of bytes.

Entry:
  dst - destination pointer to bytes
  src - source of bytes
  num - number of bytes to copy
******************************************************************************/

byte_copy(dst, src, num)
char* dst, *src;
int num;
{
    int i;

    for (i = 0; i < num; i++) {
        *dst = *src;
        dst++;
        src++;
    }
}


/******************************************************************************
Randomly shuffle a bunch of elements in a list.

Entry:
  list - array of elements to shuffle
  num  - number of elements in list
  size - size of an element

Exit:
  list is shuffled
******************************************************************************/

shuffle(list, num, size)
char* list;
int num;
int size;
{
    int i, j;
    char* temp;
    extern double drand48();

    temp = (char*) malloc(size);

    for (i = num - 1; i > 0; i--) {

        j = drand48() * (i + 1);

        /* swap i'th and j'th elements */
        byte_copy(temp, list + i * size, size);
        byte_copy(list + i * size, list + j * size, size);
        byte_copy(list + j * size, temp, size);
    }

    free(temp);
}


/******************************************************************************
Compare the length of two edges (for sorting by length).
******************************************************************************/

int edge_compare(e1, e2)
Edge** e1, ** e2;
{
    if ((*e1)->len < (*e2)->len)
        return (-1);
    else if ((*e1)->len > (*e2)->len)
        return (1);
    else
        return (0);
}


/******************************************************************************
See if this edge is inside the polygonal boundary.

Entry:
  e - edge to test

Exit:
  returns 1 if edge is inside, 0 if outside
******************************************************************************/

int inside_boundary(e)
Edge* e;
{
    float x, y;

    /* determine midpoint of edge */

    x = (points[e->p1]->pos[X] + points[e->p2]->pos[X]) * 0.5;
    y = (points[e->p1]->pos[Y] + points[e->p2]->pos[Y]) * 0.5;

    /* see if this midpoint is inside the polygon or not */

    if (point_in_split_poly(x, y))
        return (1);
    else
        return (0);
}


/******************************************************************************
See if any points are nearly on this edge.

Entry:
  e - edge to check

Exit:
  returns 1 if there is point nearly on the edge, 0 if not
******************************************************************************/

nearly_on_edge(e)
Edge* e;
{
    int i;
    float x, y;
    float val;
    float x1, y1;
    float x2, y2;
    float dx, dy;
    float t;

    /* look at each point in the set */

    for (i = 0; i < npoints; i++) {

        /* don't test the points that form the edge */

        if (i == e->p1 || i == e->p2)
            continue;

        x = points[i]->pos[X];
        y = points[i]->pos[Y];

        val = x * e->a + y * e->b + e->c;

        if (fabs(val) > 0.001)
            continue;

        x1 = points[e->p1]->pos[X];
        y1 = points[e->p1]->pos[Y];
        x2 = points[e->p2]->pos[X];
        y2 = points[e->p2]->pos[Y];
        dx = x2 - x1;
        dy = y2 - y1;

        if (fabs(dx) > fabs(dy)) {
            t = (x - x1) / dx;
            if (t > 0 && t < 1) {
#if 0
                if (parallel_flag)
                    fprintf(stderr, "point nearly on edge: val = %g, t = %g\n", val, t);
#endif
                return (1);
            }
        } else {
            t = (y - y1) / dy;
            if (t > 0 && t < 1) {
#if 0
                if (parallel_flag)
                    fprintf(stderr, "point nearly on edge: val = %g, t = %g\n", val, t);
#endif
                return (1);
            }
        }
    }

    return (0);
}


/******************************************************************************
See if this edge intersects any of the edges already in the final collection.

Entry:
  e - edge to check

Exit:
  return 0 if there are no intersections, 1 if there is an intersection
******************************************************************************/

int any_intersection(e)
Edge* e;
{
    int i;
    Edge* ee;
    float a, b, c;
    float aa, bb, cc;
    float x1, y1;
    float x2, y2;
    float xx1, yy1;
    float xx2, yy2;
    float value1, value2;
    float small = 0.00001;

    a = e->a;
    b = e->b;
    c = e->c;

    x1 = points[e->p1]->pos[X];
    y1 = points[e->p1]->pos[Y];
    x2 = points[e->p2]->pos[X];
    y2 = points[e->p2]->pos[Y];

    /* check the one edge against all others */

    for (i = 0; i < nfinal; i++) {

        ee = final[i];

        /* don't compare edges if they come from same point */
        if (e->p1 == ee->p1 || e->p1 == ee->p2 ||
            e->p2 == ee->p1 || e->p2 == ee->p2)
            continue;

        /* Plug the endpoints of each edge into the equation for the other. */
        /* If two endpoints have opposite sign, the edge straddles the line */
        /* of the other edge.  If this happens both ways, then the */
        /* edges intersect. */

        xx1 = points[ee->p1]->pos[X];
        yy1 = points[ee->p1]->pos[Y];
        xx2 = points[ee->p2]->pos[X];
        yy2 = points[ee->p2]->pos[Y];

        value1 = (xx1 * a + yy1 * b + c) * (xx2 * a + yy2 * b + c);

        if (value1 > 0)
            continue;

        aa = ee->a;
        bb = ee->b;
        cc = ee->c;

        value2 = (x1 * aa + y1 * bb + cc) * (x2 * aa + y2 * bb + cc);

        if (value2 > 0)
            continue;

        /* check the really ugly case that the lines of the edges */
        /* are nearly coincident */

        if (fabs(value1) < small && fabs(value2) < small) {

            float pa, pb;
            float d1, d2, d3, d4;
            float t;

            /* (pa,pb) is a unit vector parallel to both edges */

            pa = b;
            pb = -a;

            /* we can find out the order of the four endpoints along */
            /* their common line by examining the dot product of them */
            /* with the vector (pa,pb) */

            d1 = x1 * pa + y1 * pb;
            d2 = x2 * pa + y2 * pb;

            if (d1 > d2) {
                t = d1;
                d1 = d2;
                d2 = t;
            }

            d3 = xx1 * pa + yy1 * pb;
            d4 = xx2 * pa + yy2 * pb;

            if (d3 > d4) {
                t = d3;
                d3 = d4;
                d4 = t;
            }

            /*
            fprintf(stderr, "yucky intersect, values: %g %g %g %g\n", d1, d2, d3, d4);
            */

            if (d1 < d3 && d3 < d2) {
                return (1);
            }

            if (d1 < d4 && d4 < d2) {
                return (1);
            }

            if (d3 < d1 && d1 < d4) {
                return (1);
            }

            if (d3 < d2 && d2 < d4) {
                return (1);
            }

            continue;
        }

        /* if we get here, both edges straddle the other's line */
        /* so signal an intersection */

        return (1);
    }

    /* if we get here, the edge doesn't intersect any other edge */

    return (0);
}


/******************************************************************************
Use a greedy algorithm to connect the points into a set of triangles.

Exit:
  returns 1 if polygon self-intersects or something else went wrong, 0 if not
******************************************************************************/

int greedy_connect()
{
    int i, j;
    int p1, p2;
    int final_goal;
    Edge* e;
    int whoops_flag;

    /* maybe rescale point positions */

    if (rescale_flag)
        rescale_points();

    /* create all edges */

    for (i = 0; i < npoints; i++)
        for (j = i + 1; j < npoints; j++)
            add_edge(i, j);

    /* allocate space for the final collection of edges */

    final_goal = 3 * (npoints - 2) - boundary_count + 3;
    final = (Edge**) malloc(sizeof(Edge*) * final_goal);

#if 0
    if (final == NULL) {
        fprintf(stderr, "Memory allocation failed with nedges = %d\n", nedges);
    }
#endif

    /* add all the polygon's edges to the final edge list */

    nfinal = 0;

    for (i = 0; i < nedges; i++) {
        p1 = edges[i]->p1;
        p2 = edges[i]->p2;
        if (points[p1]->boundary && points[p2]->boundary &&
            (p2 == p1 + 1 || (p1 == 0 && p2 == boundary_count - 1))) {

            /* return from routine if the polygon self-intersects */

            if (any_intersection(edges[i])) {
                /*
                fprintf (stderr, "self-intersection at edge #%d\n", i);
                */
                return (1);
            } else
                add_final_edge(edges[i]);
        }
    }

    /* "adjust" the lengths to make very skinny triangles less likely */

    reorder_edges();

    /* sort the edges by increasing length */

    qsort(edges, nedges, sizeof(Edge*), edge_compare);

    /* maybe randomly shuffle the edges */

    if (shuffle_flag)
        shuffle(edges, nedges, sizeof(Edge*));

#if 0
    fprintf(stderr, "%d edges on boundary\n", nfinal);
#endif

#if 0
    fprintf(stderr, "npoints nedges: %d %d\n", npoints, nedges);
#endif

    /* keep adding to the collection of final edges until done */

    for (i = 0; i < nedges; i++) {

        e = edges[i];

        /* skip over edges already in final list */

        if (e->final)
            continue;

        /* see if this edge is outside the polygonal boundary */
        /* and if so, don't add it */

        if (!inside_boundary(e))
            continue;

        /* see if any points are very nearly on this edge */

        if (nearly_on_edge(e))
            continue;

        /* if this edge doesn't intersect any final edges then */
        /* add it to the list of final edges (greedy algorithm) */

        if (any_intersection(e) == 0)
            add_final_edge(e);

        /* quit if we have enough edges */

        if (nfinal == final_goal)
            break;
    }

#if 0
    fprintf(stderr, "nfinal final_goal: %d %d\n", nfinal, final_goal);
#endif

    whoops_flag = 0;

    if (nfinal != final_goal) {
        fprintf(stderr, "Whoops, nfinal = %d, final_goal = %d\n",
                nfinal, final_goal);
        collect_triangles(1);
        whoops_flag = 1;
    } else
        whoops_flag = collect_triangles(0);

    return (whoops_flag);
}


/******************************************************************************
Collect together edges to form triangles.
******************************************************************************/

int collect_triangles(whoops_flag)
int whoops_flag;
{
    int i, j, k;
    int p1, p2, p3;
    Edge* e1, *e2, *e3;
#define TRI_GOAL_SLOP 8

    /* allocate space for triangles */

    tri_goal = nfinal - npoints + 1;
    tris = (Triangle*) malloc(sizeof(Triangle) * (tri_goal + TRI_GOAL_SLOP));
    ntris = 0;

    /* check for triangle allocation */

    if (tris == 0) {
        fprintf(stderr, "collect_triangles: cannot allocate %d tris\n", tri_goal);
        return;
    }

    /* note: the edges e1, e2, e3 are numbered so they are across the */
    /* triangle from the correspondingly numbered points p1, p2, p3 */
    /* (this is an assumption of the "maybe_make_tri" routine) */

    /* examine each edge to see what triangles it is a part of */

    for (i = 0; i < nfinal; i++) {

        e1 = final[i];

        /* remember the two points at either end of this edge */

        p2 = e1->p1;
        p3 = e1->p2;

        /* look at all edges of first point of the edge "e1" */

        for (j = 0; j < points[p2]->nedges; j++) {

            e3 = points[p2]->edges[j];

            /* we'll only look at larger numbered edges to make sure */
            /* we only try to form each triangle once */

            if ((size_t) e3 <= (size_t) e1)
                continue;

            /* pick third point from opposite side of the edge "e3" */

            if (e3->p1 == p2)
                p1 = e3->p2;
            else if (e3->p2 == p2)
                p1 = e3->p1;
            else {
                fprintf(stderr, "collect_triangles: edge doesn't agree with point\n");
                continue;
            }

            /* look at all edges of second point of the edge "e1" */

            for (k = 0; k < points[p3]->nedges; k++) {

                e2 = points[p3]->edges[k];

                /* only try to form each triangle once */

                if ((size_t) e2 <= (size_t) e1)
                    continue;

                /* if point at opposite side of "e3" is that same third point, */
                /* we may have a valid triangle */

                if (e2->p1 == p1 || e2->p2 == p1) {
                    maybe_make_tri(p1, p2, p3, e1, e2, e3);
                    break;
                }
            }
        }
    }

#if 0
    fprintf(stderr, "ntris tri_goal: %d %d\n", ntris, tri_goal);
#endif

    /* consistancy check on number of triangles */

    if (ntris != tri_goal) {
        fprintf(stderr, "Whoops, ntris = %d, tri_goal = %d\n", ntris, tri_goal);
        whoops_flag = 1;
    }

    /* orient all triangles the same direction as the boundary polygon */

    orient_triangles();

    /* debug drawing */

    if (whoops_flag) {

#if 0
        clear_screen();

        fprintf(stderr,
                "npoints = %d, boundary_count = %d\n", npoints, boundary_count);

#if 1
        for (i = 0; i < nfinal; i++)
            draw_edge(final[i], color1);
#endif

#if 1
        for (i = 0; i < ntris; i++)
            draw_small_tri(tris[i].p1, tris[i].p2, tris[i].p3, color2);
#endif

        for (i = 0; i < npoints; i++) {
            fprintf(stderr, "index x y: %4d %f %f\n", points[i]->index,
                    points[i]->pos[X], points[i]->pos[Y]);
            x = points[i]->pos[X];
            y = points[i]->pos[Y];
            draw_fat_pixel((int) x, (int) y, color2);
        }

        fprintf(stderr, "next> ");
        gets(line);
#endif
    }

    return (whoops_flag);
}


/******************************************************************************
Make a triangle out of three points if the triangle does not enclose any
other points.

Entry:
  p1,p2,p3 - points of prospective triangle
  e1,e2,e3 - edges of triangle
******************************************************************************/

maybe_make_tri(p1, p2, p3, e1, e2, e3)
int p1, p2, p3;
Edge* e1, *e2, *e3;
{
    int i;
    float v;
    float v1, v2, v3;
    float* pos;

    /* calculate values of the points by plugging into the linear edge */
    /* equations, to be used to find which side of the edges other */
    /* points are on */

    v1 = points[p1]->pos[X] * e1->a + points[p1]->pos[Y] * e1->b + e1->c;
    v2 = points[p2]->pos[X] * e2->a + points[p2]->pos[Y] * e2->b + e2->c;
    v3 = points[p3]->pos[X] * e3->a + points[p3]->pos[Y] * e3->b + e3->c;

    /* see if any other points are inside the potential triangle */

    for (i = 0; i < npoints; i++) {

        /* don't look at points of triangle */
        if (i == p1 || i == p2 || i == p3)
            continue;

        pos = points[i]->pos;

        v = pos[X] * e1->a + pos[Y] * e1->b + e1->c;
        if (v * v1 < 0)
            continue;

        v = pos[X] * e2->a + pos[Y] * e2->b + e2->c;
        if (v * v2 < 0)
            continue;

        v = pos[X] * e3->a + pos[Y] * e3->b + e3->c;
        if (v * v3 < 0)
            continue;

        /* if we get here, the point is inside the triangle, so this */
        /* triangle is no good */

        return;
    }

    /* if we get here then the triangle is okay, so add it to the */
    /* list of triangles */

    if (ntris < tri_goal + TRI_GOAL_SLOP) {
        tris[ntris].p1 = p1;
        tris[ntris].p2 = p2;
        tris[ntris].p3 = p3;

        tris[ntris].e1 = e1;
        tris[ntris].e2 = e2;
        tris[ntris].e3 = e3;

        ntris++;
    } else {
        fprintf(stderr, "Whoa!  Way too many triangles in maybe_make_tri\n");
    }

    /* maybe draw triangle */

#if 0
    if (drawing_flag)
        draw_small_tri(p1, p2, p3, color1);
#endif

}


/******************************************************************************
Make all the created triangles oriented the same way (clockwise vs. counter-
clockwise) as the original boundary polygon.
******************************************************************************/

new_not_used_orient_triangles()
{
    int i;
    int p1, p2, p3;
    int b1, b2, b3;
    int orient;
    int t;
    int one, neg;

    /* search for a triangle with two boundary points */

    one = neg = 0;

    for (i = 0; i < ntris; i++) {

        p1 = tris[i].p1;
        p2 = tris[i].p2;
        p3 = tris[i].p3;
        b1 = points[p1]->boundary;
        b2 = points[p2]->boundary;
        b3 = points[p3]->boundary;

        if (b1 && b2) {
            if (p2 == (p1 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            } else if (p1 == (p2 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            }
        }

        if (b2 && b3) {
            if (p3 == (p2 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            } else if (p2 == (p3 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            }
        }

        if (b3 && b1) {
            if (p1 == (p3 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            } else if (p3 == (p1 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                if (orient == 1)
                    one++;
                else
                    neg++;
                continue;
            }
        }
    }

    /* make sure we've found a border triangle */

    if (one + neg == 0) {
        fprintf(stderr, "orient_triangles: couldn't find triangle on boundary\n");
        fprintf(stderr, "ntris = %d\n", ntris);
    }

    /* pick the majority orientation */

    if (one >= neg)
        orient = 1;
    else
        orient = -1;

    /* go through all triangles and flip them to match the orientation */
    /* of the border triangle */

    for (i = 0; i < ntris; i++) {
        t = triangle_direction(&tris[i]);
        if (t != orient)
            flip_triangle(&tris[i]);
    }
}


/******************************************************************************
Make all the created triangles oriented the same way (clockwise vs. counter-
clockwise) as the original boundary polygon.
******************************************************************************/

orient_triangles()
{
    int i;
    int found;
    int p1, p2, p3;
    int b1, b2, b3;
    int orient;
    int t;

    /* search for a triangle with two boundary points */

    found = 0;

    for (i = 0; i < ntris; i++) {

        p1 = tris[i].p1;
        p2 = tris[i].p2;
        p3 = tris[i].p3;
        b1 = points[p1]->boundary;
        b2 = points[p2]->boundary;
        b3 = points[p3]->boundary;

        if (b1 && b2) {
            if (p2 == (p1 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            } else if (p1 == (p2 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            }
        }

        if (b2 && b3) {
            if (p3 == (p2 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            } else if (p2 == (p3 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            }
        }

        if (b3 && b1) {
            if (p1 == (p3 + 1) % boundary_count) {     /* correct direction */
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            } else if (p3 == (p1 + 1) % boundary_count) { /* backwards */
                flip_triangle(&tris[i]);
                orient = triangle_direction(&tris[i]);
                found = 1;
                break;
            }
        }
    }

    /* make sure we've found a border triangle */

    if (!found) {
        fprintf(stderr, "orient_triangles: couldn't find triangle on boundary\n");
        fprintf(stderr, "ntris = %d\n", ntris);
        /*
        exit (-1);
        */
    }

    /* go through all triangles and flip them to match the orientation */
    /* of the border triangle */

    for (i = 0; i < ntris; i++) {
        t = triangle_direction(&tris[i]);
        if (t != orient)
            flip_triangle(&tris[i]);
    }
}


/******************************************************************************
Return which direction a triangle is facing (clockwise or counter-clockwise).

Entry:
  tri - triangle to find orientation of

Exit:
  returns 1 if one way, -1 for other (I'm too lazy to figure which is which)
******************************************************************************/

int triangle_direction(tri)
Triangle* tri;
{
    Vector v1, v2, v3;
    Vector cross;

    /* get the triangle vertices */

    vcopy(points[tri->p1]->pos, v1);
    vcopy(points[tri->p2]->pos, v2);
    vcopy(points[tri->p3]->pos, v3);

    /* find directions of two edges and take their crossproduct */

    vsub(v2, v1, v2);
    vsub(v3, v1, v3);
    vcross(v2, v3, cross);

    /* return the orientation */

    if (cross[Z] < 0)
        return (-1);
    else if (cross[Z] > 0)
        return (1);
    else {
        fprintf(stderr, "triangle_direction: degenerate triangle\n");
        return (1);  /* return something arbitrary */
    }
}


/******************************************************************************
Flip the order of the vertices and edges in a triangle.
******************************************************************************/

flip_triangle(tri)
Triangle* tri;
{
    int temp;
    Edge* etemp;

    /* swap p2 and p3 */

    temp = tri->p2;
    tri->p2 = tri->p3;
    tri->p3 = temp;

    /* swap e2 and e3 */

    etemp = tri->e2;
    tri->e2 = tri->e3;
    tri->e3 = etemp;
}


/******************************************************************************
Return the number of triangles that were formed.
******************************************************************************/

int get_ntris()
{
    return (ntris);
}


/******************************************************************************
Return the indices of a given triangle.

Entry:
  num - number of the desired triangle

Exit:
  p1,p2,p3 - indices of the triangle's vertices
******************************************************************************/

int get_triangle(num, p1, p2, p3)
int num;
int* p1, *p2, *p3;
{
    *p1 = points[tris[num].p1]->index;
    *p2 = points[tris[num].p2]->index;
    *p3 = points[tris[num].p3]->index;
}


/******************************************************************************
Print polygon info.
******************************************************************************/

print_poly()
{
    int i;

    fprintf(stderr, "boundary_count: %d\n", boundary_count);

    for (i = 0; i < boundary_count; i++) {
        fprintf(stderr, " x y: %f %f\n", points[i]->pos[X], points[i]->pos[Y]);
    }
}


/******************************************************************************
Determine which triangle a point is in.

Entry:
  x,y - position of point

Exit:
  b1,b2,b3 - barycentric coordinates of point in the triangle
  returns index of triangle the point is in, or -1 if error
******************************************************************************/

int point_in_which_triangle(x, y, b1, b2, b3)
float x, y;
float* b1, *b2, *b3;
{
    int i;
    int result;
    Triangle* tri;
    Edge* e1, *e2, *e3;
    int p1, p2, p3;
    float v1, v2, v3;
    float v;
    int found;
    int index;
    Vector n, v1p, v12, v13;
    Vector sc, tc;
    float r, s, t;
    float* pos1, *pos2, *pos3;
    float len;

    /* check to see that the point is inside the polygon */

    result = point_in_split_poly(x, y);

    if (!result) {
        fprintf(stderr, "point_in_which_triangle: point not in polygon\n");
        fprintf(stderr, "x y: %f %f\n", x, y);
        for (i = 0; i < npoints; i++) {
            fprintf(stderr, "poly x y: %f %f\n",
                    points[i]->pos[X], points[i]->pos[Y]);
        }
        return (-1);
    }

    /* look through each triangle in set */

    found = 0;

    for (i = 0; i < ntris; i++) {

        tri = &tris[i];

        e1 = tri->e1;
        e2 = tri->e2;
        e3 = tri->e3;

        p1 = tri->p1;
        p2 = tri->p2;
        p3 = tri->p3;

        v1 = points[p1]->pos[X] * e1->a + points[p1]->pos[Y] * e1->b + e1->c;
        v2 = points[p2]->pos[X] * e2->a + points[p2]->pos[Y] * e2->b + e2->c;
        v3 = points[p3]->pos[X] * e3->a + points[p3]->pos[Y] * e3->b + e3->c;

        /* see if the point is inside this triangle */

        v = x * e1->a + y * e1->b + e1->c;
        if (v * v1 < 0)
            continue;

        v = x * e2->a + y * e2->b + e2->c;
        if (v * v2 < 0)
            continue;

        v = x * e3->a + y * e3->b + e3->c;
        if (v * v3 < 0)
            continue;

        /* if we get here, the point is inside the triangle */

        if (found) {
            fprintf(stderr, "point_in_which_triangle: point in two tris\n");
            return (-1);
        } else {
            found = 1;
            index = i;
        }
    }

    /* we should have found a triangle that the point is in */

    if (!found) {
        fprintf(stderr, "point_in_which_triangle: point not in any tri\n");
        return (-1);
    }

#if 0
    draw_small_tri(tris[index].p1, tris[index].p2, tris[index].p3, color1);
#endif

    /* determine barycentric cooridinates of point in the triangle */

    tri = &tris[index];

    p1 = tri->p1;
    p2 = tri->p2;
    p3 = tri->p3;

    pos1 = points[p1]->pos;
    pos2 = points[p2]->pos;
    pos3 = points[p3]->pos;

    v12[X] = pos2[X] - pos1[X];
    v12[Y] = pos2[Y] - pos1[Y];
    v12[Z] = 0;

    v13[X] = pos3[X] - pos1[X];
    v13[Y] = pos3[Y] - pos1[Y];
    v13[Z] = 0;

    v1p[X] = x - pos1[X];
    v1p[Y] = y - pos1[Y];
    v1p[Z] = 0;

    vcross(v12, v13, n);
    vcross(v1p, v13, sc);
    vcross(v12, v1p, tc);
    len = vlen(n);

    s = fabs(sc[Z]) / len;
    t = fabs(tc[Z]) / len;

    r = 1 - (s + t);

    *b1 = r;
    *b2 = s;
    *b3 = t;

#if 0
    fprintf(stderr, "r s t: %f %f %f\n", r, s, t);

    xx = pos1[X] * r + pos2[X] * s + pos3[X] * t;
    yy = pos1[Y] * r + pos2[Y] * s + pos3[Y] * t;
    zz = pos1[Z] * r + pos2[Z] * s + pos3[Z] * t;

    writepixel(xx, yy, color1);
#endif

    return (index);
}


/******************************************************************************
Compute a line passing through two points.

Entry:
  x1,y1 - first point to put line through
  x2,y2 - second point

Exit:
  aa,bb,cc - equation of line
******************************************************************************/

compute_line(x1, y1, x2, y2, aa, bb, cc)
float x1, y1;
float x2, y2;
float* aa, *bb, *cc;
{
    float a, b, c;
    float len;

    a = y2 - y1;
    b = x1 - x2;
    c = y1 * x2 - x1 * y2;

    len = sqrt(a * a + b * b);

    if (len > 0) {
        a /= len;
        b /= len;
        c /= len;
    }

    *aa = a;
    *bb = b;
    *cc = c;
}


/******************************************************************************
Re-scale the points so they fit in the window.
******************************************************************************/

rescale_points()
{
    int i;
    float x, y;
    float xmin, xmax;
    float ymin, ymax;
    float dx, dy;
    float cx, cy;
    float scale;

    xmin = ymin = 1e20;
    xmax = ymax = -1e20;

    for (i = 0; i < npoints; i++) {
        x = points[i]->pos[X];
        y = points[i]->pos[Y];
        if (x < xmin) xmin = x;
        if (y < ymin) ymin = y;
        if (x > xmax) xmax = x;
        if (y > ymax) ymax = y;
    }

    dx = fabs(xmax - xmin);
    dy = fabs(ymax - ymin);

    if (x_screen / dx < y_screen / dy)
        scale = 0.9 * x_screen / dx;
    else
        scale = 0.9 * y_screen / dy;

    cx = (xmin + xmax) * 0.5;
    cy = (ymin + ymax) * 0.5;

    for (i = 0; i < npoints; i++) {
        points[i]->pos[X] = (points[i]->pos[X] - cx) * scale + x_screen * 0.5;
        points[i]->pos[Y] = (points[i]->pos[Y] - cy) * scale + y_screen * 0.5;
    }
}


/******************************************************************************
Check to see that none of the old polygons in an area to be re-tiled were
folded over one another.

Entry:
  x,y - projected position of central point

Exit:
  returns 1 if there was folding over, 0 if none of the polygons overlapped
******************************************************************************/

int fold_in_poly_check(x, y)
float x, y;
{
    int i;
    Vector n;
    Vector v1, v2;
    Vector center;
    float cross1, cross2;

    /* We're going to determine folds by examining the cross-products of */
    /* the vectors radiating from the central point.  If the sign of this */
    /* ever changes, then there is a fold. */

    center[X] = x;
    center[Y] = y;
    center[Z] = 0;

    /* determine sign of initial cross-product */

    vsub(points[0]->pos, center, v1);
    vsub(points[1]->pos, center, v2);
    vcross(v1, v2, n);
    cross1 = n[Z];

    if (fabs(cross1) < 0.0001) {
        fprintf(stderr, "small area in fold_in_poly_check: %f\n", cross1);
        return (1);
    }

    /* examine other cross-products */

    for (i = 1; i < npoints; i++) {

        /* get next cross-product */

        vsub(points[i]->pos, center, v1);
        vsub(points[(i + 1) % npoints]->pos, center, v2);
        vcross(v1, v2, n);
        cross2 = n[Z];

        /* see if the sign changed */

        if (cross1 * cross2 < 0)
            return (1);
    }

    return (0);
}


/******************************************************************************
See if point is in the given polygon.
******************************************************************************/

int point_in_split_poly(x, y)
float x, y;
{
    int result;

    result = point_in_poly(x, y, boundary_count, points);

    return (result);
}


/******************************************************************************
Point in polygon test from Ken McElvain.

Entry:
  x,y     - point to test against polygon
  cnt     - number of vertices in polygon
  polypts - vertices of polygon

Exit:
  returns zero if point is outside polygon, non-zero for points in polygon
******************************************************************************/

int point_in_poly(x, y, cnt, polypts)
float x, y;
int cnt;
Point** polypts;
{
    int i;
    int oldquad, newquad;
    float* thispt, *lastpt;
    float a, b;     /* these can be int's if x and y are */
    int wind;       /* current winding number */

    wind = 0;
    lastpt = polypts[cnt - 1]->pos;
    oldquad = whichquad(lastpt, x, y);  /* get starting angle */

    /* examine each point of the polygon and follow the winding number */
    /* as we move around the given point to test */

    for (i = 0; i < cnt; i++) {

        thispt = polypts[i]->pos;
        newquad = whichquad(thispt, x, y);

        /* adjust the winding number */

        if (oldquad != newquad) {

            /* Use mod 4 comparsions to see if we have advanced or */
            /* backed up one quadrant.  */

            if (((oldquad + 1) & 3) == newquad)
                wind++;
            else if (((newquad + 1) & 3) == oldquad)
                wind--;
            else {

                /* Upper left to lower right, or upper right to lower left. */
                /* Determine direction of winding by intersection with x = 0.  */

                a = lastpt[Y] - thispt[Y];
                a *= (x - lastpt[X]);
                b = lastpt[X] - thispt[X];
                a += lastpt[Y] * b;
                b *= y;

                if (a > b) wind += 2;
                else wind -= 2;
            }
        }
        lastpt = thispt;
        oldquad = newquad;
    }

    /* non-zero return value means the point is in the polygon */

    return (wind);
}


/******************************************************************************
Figure out which quadrent pt is in with respect to (x,y).
******************************************************************************/

int whichquad(pt, x, y)
Vector pt;
float x, y;
{
    int quadrant;

    if (pt[X] < x) {
        if (pt[Y] < y)
            quadrant = 2;
        else
            quadrant = 1;
    } else {
        if (pt[Y] < y)
            quadrant = 3;
        else
            quadrant = 0;
    }

    return (quadrant);
}

/******************************************************************************
Give transformation from arbitrary plane into xy-plane, and also give the
inverse transformation.

Entry:
  a,b,c,d - coefficients of plane equation, (a,b,c) assumed unit vector

Exit:
  mat  - transformation from arbitrary plane to xy-plane
  imat - inverse transformation
  returns 1 if we've got bad plane equation, 0 if everything is okay
******************************************************************************/

int face_to_xy_plane(a, b, c, d, mat, imat)
float a, b, c, d;
Matrix mat, imat;
{
    int i;
    float dx, dy, dz;
    Vector v1, v2, v3;
    Matrix t;

    /* translational component of matrix */

    dx = a * d;
    dy = b * d;
    dz = c * d;

    /* rotational component can be any orthonormal matrix */
    /* where this is the third column */

    v3[X] = a;
    v3[Y] = b;
    v3[Z] = c;

    /* make a second vector that's not a multiple of v3 */

    vcopy(v3, v2);

    if (v2[X] != 0)
        v2[Y] += 1;
    else if (v2[Y] != 0)
        v2[Z] += 1;
    else if (v2[Z] != 0)
        v2[X] += 1;
    else {
        fprintf(stderr, "degenerate plane equation in face_to_xy_plane\n");
        fprintf(stderr, "a b c d: %g %g %g %g\n", a, b, c, d);
        return (1);
    }

    /* do a few cross products to make the vectors orthonormal */

    vcross(v2, v3, v1);
    vnorm(v1);
    vcross(v3, v1, v2);
    vnorm(v2);

    /* make a matrix out of the vectors */

    mat_ident(mat);
    for (i = 0; i < 3; i++) {
        mat[i][0] = v1[i];
        mat[i][1] = v2[i];
        mat[i][2] = v3[i];
    }

    /* now create the matrix that does the trick */

    mat_translate(t, dx, dy, dz);
    mat_mult(mat, mat, t);
    mat_copy(imat, mat);
    mat_invert(imat);

    return (0);
}


/******************************************************************************
Re-order the edges by down-weighting edges that are nearly parallel to the
boundary.
******************************************************************************/

reorder_edges()
{
    int i, j;
    float max_len;
    Edge* e;
    Edge* edge;
    Point* p1, *p2;
    float dot, dot_max;
    float near_one = 0.99;

    /* find the maximum edge length */
    max_len = -1e20;
    for (i = 0; i < nedges; i++)
        if (edges[i]->len > max_len)
            max_len = edges[i]->len;

#if 0
    printf("max_len = %f\n", max_len);
#endif

    max_len *= 2;

    /* examine each edge to see if it is nearly parallel to the */
    /* boundary of the polygon */

    for (i = 0; i < nedges; i++) {

        edge = edges[i];

        /* don't look at boundary edges */
        if (edge->final)
            continue;

        p1 = points[edge->p1];
        p2 = points[edge->p2];

        /* each of these points should have exactly two final edges */
        /* in their lists, and these are the boundary edges */

        if (p1->nedges != 2 || p2->nedges != 2) {
            fprintf(stderr, "reorder_edges: bad number of edges: %d %d\n",
                    p1->nedges, p2->nedges);
            exit(-1);
        }

        /* compare the boundary edges with the edge in question */
        /* and if they are close to parallel, increase the "length" */
        /* of the edge to make it unlikely to be used */

        dot_max = -1;

        for (j = 0; j < p1->nedges; j++) {
            e = p1->edges[j];
            dot = fabs(e->a * edge->a + e->b * edge->b);
#if 0
            printf("dot = %f\n", dot);
#endif
            if (dot > dot_max)
                dot_max = dot;
        }

        for (j = 0; j < p2->nedges; j++) {
            e = p2->edges[j];
            dot = fabs(e->a * edge->a + e->b * edge->b);
#if 0
            printf("dot = %f\n", dot);
#endif
            if (dot > dot_max)
                dot_max = dot;
        }

        /* this new "length" penalizes edges for being nearly parallel to */
        /* any of the boundary edges, and this penalty is greater for those */
        /* edges that are more nearly parallel to the boundary */

        if (dot_max > near_one) {
            edge->len = max_len * (1.0 + dot_max);
        }

#if 0
        printf("len = %f\n", edge->len);
#endif

    }
}

