#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "clipper.h"

#define info_fmt(x)   "[INFO] " x
#define warn_fmt(x)   "[WARN] " x
#define erro_fmt(x)   "[ERROR] " x
#define trc_fmt(x) "[INFO] %s:%d " x
#define clip_trc(fmt, args...) printf(trc_fmt(fmt), __func__, __LINE__, ##args)
#define clip_log(fmt, args...) printf(info_fmt(fmt), ##args)
#define clip_warn(fmt, args...) printf(warn_fmt(fmt), ##args)
#define clip_err(fmt, args...) printf(erro_fmt(fmt), ##args)
#define clip_info_h(fmt, args...)
#define clip_info_l(fmt, args...)

#define FIXED_POINT_SHIFT (8)
#define FIXED_POINT_ONE (1 << FIXED_POINT_SHIFT)
#define FIXED_POINT_HALF (1 << (FIXED_POINT_SHIFT-1))

#define ROUND_DOWN(x) ((x + FIXED_POINT_HALF) >> FIXED_POINT_SHIFT)

#define PUSH_BACK_THRES (8)
#define PUSH_BACK_MASK (7)
#define clip_assert(cond, msg) \
    do {\
        if (!(cond)) {\
            clip_log("%s",msg);\
            assert(cond);\
        }\
    } while(0)

#define poly_free(x) (free(x))
#define sign(x) (((x) < 0) ? -1 : 1)
float max(float a, float b)
{
    return ( a > b ) ? a : b;
}

float min(float a, float b)
{
    return ( a < b ) ? a : b;
}

float mean3way(float a, float b, float c)
{
    return (a + b + c) * 0.333333f;    
}

#define define_clamp(type) \
inline type clamp_##type(type a, type low, type high) \
{ \
    return ( a < low ) ? low : ( ( a > high ) ? high : a ); \
}

#define clampT(type) clamp_##type
define_clamp(int)
define_clamp(float)

int clamp(int a, int low, int high)
{
    return ( a < low ) ? low : ( ( a > high ) ? high : a );
}

const char *poly_type_str[] = {
    "Undetermined",
    "Convex",
    "Concave"
};

typedef union {
    Point *pts[2];
    struct {
        Point *prev;
        Point *curr;
    };
    struct {
        Point *from;
        Point *to;
    };
} Line_Ptr;


typedef struct Intersection Intersection;
typedef struct Edge_Ptr Edge_Ptr;

struct Edge_Ptr {
    Point *prev;
    Point *curr;
    Intersection **intersection_point;
    Edge_Ptr *next;
    int size;
};

typedef struct {
    Edge_Ptr edges[MAX_POINTS];
    int size; 
} Edge_Ptr_List;

struct Intersection{
    int type;
    Edge_Ptr *edges[2];
    Point pt;
};

typedef struct {
    Intersection pts[MAX_POINTS];
    int size;
} Intersection_Point_List;

typedef enum { VERT_UNDERTERMINED=-1, VERT_NORM=0, VERT_IN, VERT_OUT, VERT_SAME, VERT_CHECKED, VERT_COMPLETE } Vert_Type;


void edge_freeList(Edge_Ptr_List *list)
{
    for (int i = 0; i < list->size; ++i) {
        Edge_Ptr *edge = &list->edges[i];
        if (edge->size)
            free(edge->intersection_point);
    }
    free(list);
}

void edge_pushBackIntersection(Edge_Ptr *edge, Intersection *pt)
{
    if (edge->size == 0) {
        edge->intersection_point = (Intersection**)calloc(sizeof(Intersection*)*PUSH_BACK_THRES,1);
    } else if ((edge->size & PUSH_BACK_MASK) == 0) {
        memmove((void *)edge->intersection_point, (void *)edge->intersection_point, sizeof(Intersection*)*PUSH_BACK_THRES);
    }
    edge->intersection_point[edge->size++] = pt;
}

static int line_onSegment(const Point* p, const Point* q, const Point* r) 
{ 
    if (q->x <= max(p->x, r->x) && q->x >= min(p->x, r->x) && 
            q->y <= max(p->y, r->y) && q->y >= min(p->y, r->y)) 
        return 1; 
    return 0; 
} 

typedef enum {
    LINE_SAME = 0,
    LINE_LEFT = 1,
    LINE_RIGHT = 2,
    LINE_IN = 1,
    LINE_OUT = 2
} LINE_DIRECTION;

typedef enum {
    LINE_CLOCKWISE = 1,
    LINE_COUNTER = 2
} LINE_ORIENTATION;

/***
 *@brief To find orientation of ordered triplet (p, q, r). 
 *@param[in] p line start point
 *@param[in] q line end point
 *@param[out] r query point
 *@retval 0 : p, q and r are colinear, 1 : Clockwise, 2 : Counterclockwise
 *@retval 0 : p, q and r are colinear, 1 : left, 2 : right
**/
int line_orientation(const Point* p, const Point* q, const Point* r) 
{
    int val = (q->y - p->y) * (r->x - q->x) - 
              (q->x - p->x) * (r->y - q->y); 
  
    if (val == 0) return 0;  // colinear 
    return (val > 0)? 1: 2; // clock or counterclock wise 
}

// p - start point 
// q - end point
// r - query point 
// output 0 : on the line , 1 : 
int line_orientation_v2(const Point* p, const Point* q, const Point* r) 
{
    int val = (q->y - p->y) * (r->x - q->x) - 
              (q->x - p->x) * (r->y - q->y); 
  
    return (val <= 0)? LINE_LEFT: LINE_RIGHT; // clock or counterclock wise 
}

int line_orientation_v3(const Point* p, const Point* q, const Point* r) 
{
    int val = (q->y - p->y) * (r->x - q->x) - 
              (q->x - p->x) * (r->y - q->y); 
  
    return (val < 0)? LINE_LEFT: ((val == 0) ? LINE_SAME : LINE_RIGHT); // clock or counterclock wise 
}

// The function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
int line_doIntersect(const Point *p1, const Point *q1, const Point *p2, const Point *q2) 
{ 
    // Find the four orientations needed for general and 
    // special cases 
    int o1 = line_orientation(p1, q1, p2); 
    int o2 = line_orientation(p1, q1, q2); 
    int o3 = line_orientation(p2, q2, p1); 
    int o4 = line_orientation(p2, q2, q1); 
  
    // General case 
    if (o1 != o2 && o3 != o4) 
        return 1; 
  
    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && line_onSegment(p1, p2, q1)) return 1; 
  
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && line_onSegment(p1, q2, q1)) return 1; 
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && line_onSegment(p2, p1, q2)) return 1; 
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && line_onSegment(p2, q1, q2)) return 1; 
  
    return 0; // Doesn't fall in any of the above cases 
} 
  
int line_doIntersectv2(const Point *p1, const Point *q1, const Point *p2, const Point *q2) 
{ 
    // Find the four orientations needed for general and 
    // special cases 
    #if 0
    if (!memcmp(p1, p2, sizeof(Point)) || !memcmp(q1, p2, sizeof(Point)) ||
        !memcmp(p1, q2, sizeof(Point)))
    {
        
        clip_trc("%.4f, %.4f -  %.4f, %.4f vs  %.4f, %.4f - %.4f, %.4f \n",
            p1->x, p1->y, q1->x, q1->y, p2->x, p2->y, q2->x, q2->y);
        clip_trc("%d %d\n",memcmp(p1, p2, sizeof(Point)), memcmp(p1, q2, sizeof(Point)));
        return 0; 
    }
    #endif
    int o1 = line_orientation(p1, q1, p2); 
    int o2 = line_orientation(p1, q1, q2); 
    int o3 = line_orientation(p2, q2, p1); 
    int o4 = line_orientation(p2, q2, q1); 
  
    // General case
    clip_info_l("o[%d %d %d %d]\n", o1, o2, o3, o4);
    if (o1 == 0 || o3 == 0)
        return 0;
    if (o1 != o2 && o3 != o4)
        return 1; 
  
    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && line_onSegment(p1, p2, q1)) return 1; 
  
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && line_onSegment(p1, q2, q1)) return 1; 
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && line_onSegment(p2, p1, q2)) return 1; 
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && line_onSegment(p2, q1, q2)) return 1; 
  
    return 0; // Doesn't fall in any of the above cases 
} 
// Returns true if the point p lies inside the polygon[] with n vertices 
/*int isInside(const Polygon* polygon, int n, const Point *p) 
{ 
    // There must be at least 3 vertices in polygon[] 
    if (n < 3)  return false; 
  
    // Create a point for line segment from p to infinite 
    Point extreme = {INF, p.y}; 
  
    // Count intersections of the above line with sides of polygon 
    int count = 0, i = 0; 
    do
    { 
        int next = (i+1)%n; 
  
        // Check if the line segment from 'p' to 'extreme' intersects 
        // with the line segment from 'polygon[i]' to 'polygon[next]' 
        if (doIntersect(polygon[i], polygon[next], p, extreme)) 
        { 
            // If the point 'p' is colinear with line segment 'i-next', 
            // then check if it lies on segment. If it lies, return true, 
            // otherwise false 
            if (orientation(polygon[i], p, polygon[next]) == 0) 
               return onSegment(polygon[i], p, polygon[next]); 
  
            count++; 
        } 
        i = next;
    } while (i != 0); 
  
    // Return true if count is odd, false otherwise 
    return count&1;  // Same as (count%2 == 1) 
} 
*/
/*
static Point line_computeIntersection(const Point* pt1, const Point* pt2, const Point edge[2])
{
    Point intersection_point = { };
    const Point *pt3 = edge;
    const Point *pt4 = edge + 1;

    float x1y2_y1x2 = pt1->x * pt2->y - pt1->y * pt2->x;
    float x3_x4 = pt3->x - pt4->x;
    float x1_x2 = pt1->x - pt2->x;
    float x3y4_y3x4 = pt3->x * pt4->y - pt3->y * pt4->x;
    float y3_y4 = pt3->y - pt4->y;
    float y1_y2 = pt1->y - pt2->y;

    float denom = x1_x2 * y3_y4 - y1_y2 * x3_x4;

    if (denom == 0.0f ) return intersection_point;

    intersection_point.x =
    (x1y2_y1x2 * x3_x4 - x1_x2 * x3y4_y3x4) /
    denom;
    intersection_point.y =
    (x1y2_y1x2 * y3_y4 - y1_y2 * x3y4_y3x4) /
    denom;

    return intersection_point;
}
*/
static Point line_computeIntersection_v2(const Point* pt1, const Point* pt2, const Line_Ptr *edge)
{
    Point intersection_point = { };
    const Point *pt3 = edge->pts[0];
    const Point *pt4 = edge->pts[1];

    float x1y2_y1x2 = pt1->x * pt2->y - pt1->y * pt2->x;
    float x3_x4 = pt3->x - pt4->x;
    float x1_x2 = pt1->x - pt2->x;
    float x3y4_y3x4 = pt3->x * pt4->y - pt3->y * pt4->x;
    float y3_y4 = pt3->y - pt4->y;
    float y1_y2 = pt1->y - pt2->y;

    float denom = x1_x2 * y3_y4 - y1_y2 * x3_x4;

    if (denom == 0.0f ) return intersection_point;

    intersection_point.x =
    (x1y2_y1x2 * x3_x4 - x1_x2 * x3y4_y3x4) /
    denom;
    intersection_point.y =
    (x1y2_y1x2 * y3_y4 - y1_y2 * x3y4_y3x4) /
    denom;

    return intersection_point;
}

static int inRange(float a, float b, float c)
{
    float rmax = max(b, c);
    float rmin = min(b, c);
    return !(a < rmin || a > rmax);
}

static Intersection line_computeIntersection_v3(Edge_Ptr *edge, Edge_Ptr *query_edge)
{
    Intersection intersection_point = {VERT_UNDERTERMINED, {NULL, NULL}, {-1, -1} };
    const Point *pt1 = edge->prev;
    const Point *pt2 = edge->curr;
    const Point *pt3 = query_edge->prev;
    const Point *pt4 = query_edge->curr;

    // filter out of box
    //if (pt4->x < pt1->x || pt4->y < pt1->y ||
     //   pt3->x > pt2->x || pt3->y > pt2->y)
     //   return intersection_point;

    int inter = line_doIntersectv2(pt1, pt2, pt3, pt4);
    if (!inter) return intersection_point;

    int pt3_orig = line_orientation_v3(pt1, pt2, pt3);
    int pt4_orig = line_orientation_v3(pt1, pt2, pt4);

    if (pt4_orig == LINE_SAME) {
        intersection_point.type = VERT_SAME;
        intersection_point.pt = *pt4;
        intersection_point.edges[0] = edge;
        intersection_point.edges[1] = query_edge;
        return intersection_point;
    }

    float x1y2_y1x2 = pt1->x * pt2->y - pt1->y * pt2->x;
    float x3_x4 = pt3->x - pt4->x;
    float x1_x2 = pt1->x - pt2->x;
    float x3y4_y3x4 = pt3->x * pt4->y - pt3->y * pt4->x;
    float y3_y4 = pt3->y - pt4->y;
    float y1_y2 = pt1->y - pt2->y;

    float denom = x1_x2 * y3_y4 - y1_y2 * x3_x4;

    if (denom == 0.0f ) return intersection_point;

    intersection_point.type = (pt4_orig == LINE_LEFT) ? VERT_IN : VERT_OUT;

    intersection_point.edges[0] = edge;
    intersection_point.edges[1] = query_edge;

    intersection_point.pt.x =
    (x1y2_y1x2 * x3_x4 - x1_x2 * x3y4_y3x4) /
    denom;
    intersection_point.pt.y =
    (x1y2_y1x2 * y3_y4 - y1_y2 * x3y4_y3x4) /
    denom;

    return intersection_point;
}

// To solve example6_2
static int line_isIntersectUnique(const Intersection *inter, const Intersection_Point_List *list)
{
    for (int i = 0; i < list->size; ++i) {
        if (inter->pt.x == list->pts[i].pt.x &&
            inter->pt.y == list->pts[i].pt.y) {
            return 0;
        }
    }
    return 1;
}

static Array2di arr_init2di(int w, int h)
{
    Array2di arr;
    arr.w = w;
    arr.h = h;
    arr.data = (int*)calloc(w*h*sizeof(int),1);
    arr.size = w * h;
    return arr;
}

static inline int arr_get(Array2di *arr, int x, int y)
{
    return arr->data[y*arr->w+x];
}

static inline void arr_fill(Array2di *arr, int x, int y, int val)
{
    arr->data[y*arr->w+x] = val;
}


static void arr_free(Array2di arr)
{
    free(arr.data);
}

static Line poly_getLine(const Polygon *polygon, int idx)
{
    assert(idx < polygon->size);
    Line edge = { 
                .pts = {
                polygon->pts[(idx + polygon->size - 1) % polygon->size],
                polygon->pts[idx]
                }
            };
    return edge;
}

static Line_Ptr poly_getLinePtr(Polygon *polygon, int idx)
{
    assert(idx < polygon->size);
    Line_Ptr edge = { 
                .pts = {
                &polygon->pts[(idx + polygon->size - 1) % polygon->size],
                &polygon->pts[idx]
                }
            };
    return edge;
}

static Edge_Ptr poly_getEdgePtr(Polygon *polygon, int idx)
{
    assert(idx < polygon->size);
    Edge_Ptr edge = { 
                .prev = &polygon->pts[(idx + polygon->size - 1) % polygon->size],
                .curr = &polygon->pts[idx],
                .intersection_point = NULL,
                .size = 0,
                .next = NULL
                };
    return edge;
}

static Edge_Ptr_List* poly_getEdgePtrList(Polygon *polygon)
{
    Edge_Ptr_List *edge_list = (Edge_Ptr_List*)malloc(sizeof(Edge_Ptr_List));
    edge_list->size = polygon->size;
    for (int i = 0; i < polygon->size; ++i) {
        edge_list->edges[i] = poly_getEdgePtr(polygon, i);
        int next_idx = (i + 1) % polygon->size;
        edge_list->edges[i].next = &edge_list->edges[next_idx];
    }
    return edge_list;
}

static inline void poly_empty(Polygon *polygon)
{
    polygon->size = 0;
}


static inline void poly_copy(const Polygon *src, Polygon *dst)
{
    dst->size = src->size;
    dst->type = src->type;
    memcpy(dst->pts, src->pts, sizeof(Point) * src->size);
}

static inline void poly_addPoint(Polygon* polygon, const Point *pt)
{
    assert(polygon->size < MAX_POINTS);
    polygon->pts[polygon->size++] = *pt;
    clip_info_l("poly size:%d [%.2f, %2f]\n", polygon->size, pt->x, pt->y);
}

static int poly_checkType(Polygon* polygon)
{
    const int size = polygon->size;
    for (int i = 0; i < size; ++i) {
        const Point *prev = &polygon->pts[(i + size - 1) % size];
        const Point *cur = &polygon->pts[i];
        const Point *next = &polygon->pts[(i + 1) % size];
        if (line_orientation(prev, next, cur) == LINE_RIGHT)
            return CONCAVE;
    }
    return CONVEX;
}

int compare_func_hori_acc(const void *a, const void *b)
{
    return (*(Intersection**)a)->pt.x > (*(Intersection**)b)->pt.x;
}

int compare_func_hori_dec(const void *a, const void *b)
{
    return (*(Intersection**)a)->pt.x < (*(Intersection**)b)->pt.x;
}

int compare_func_vert_acc(const void *a, const void *b)
{
    return (*(Intersection**)a)->pt.y > (*(Intersection**)b)->pt.y;
}

int compare_func_vert_dec(const void *a, const void *b)
{
    return (*(Intersection**)a)->pt.y < (*(Intersection**)b)->pt.y;
}

static float edge_calcSlope(const Edge_Ptr *edge)
{
    float x_diff = edge->curr->x - edge->prev->x;
    if (!x_diff)
        return 999999999.0f;
    return (edge->curr->y - edge->prev->y)/x_diff;
}
static void edge_sortIntersection(Edge_Ptr *edge)
{
    int (*compare_func)(const void*, const void*);
    float slope = edge_calcSlope(edge);
    int i = 0;
    if (abs(slope) > 1.0) {
        if (edge->curr->y >= edge->prev->y) {
            compare_func = compare_func_vert_acc;
            i = 0;
        } else {
            compare_func = compare_func_vert_dec;
            i = 1;
        }
    } else {
        if (edge->curr->x >= edge->prev->x) {
            compare_func = compare_func_hori_acc;
            i = 2;
        } else {
            compare_func = compare_func_hori_dec;
            i = 3;
        }
    }
    qsort(edge->intersection_point, edge->size, sizeof(Intersection**),compare_func);
}

float CLIP_getArea(const Polygon *polygon)
{
    const int number_vertex = polygon->size;
    const Point *pts = polygon->pts;
    float area = 0;
    int next = 0;
    for (int i = 0; i<number_vertex; ++i) {
        next = (i + 1) % number_vertex;
        area += (pts[i].x * (pts[next].y)) - ((pts[i].y) * pts[next].x);
    }
    return abs((float)area/2.f);
}

/*
static void vert_addVert(Vert_list *list, const Vert *vertex)
{
    assert(list->size < MAX_POINTS);
    list->vertex[list->size++] = *vertex;
}

static void vert_addPointPtr(Vert_Ptr_List *list, Point *pt, int type)
{
    assert(list->size < MAX_POINTS);
    list->pts[list->size] = pt;
    list->type[list->size++] = type;
}


static void vert_addVertPtr(Vert_Ptr_List *list, Vert *vert)
{
    assert(list->size < MAX_POINTS);
    list->pts[list->size] = &vert->pt;
    list->type[list->size++] = vert->type;
}

static Vert *vert_getLastPtr(Vert_List *list)
{
    return &list->vert[list->size-1];
}
*/

static void poly_incrementResultPoly(Clipper_result *result)
{
    if (result->size == 0) {
        result->polygon = (Polygon*)calloc(sizeof(Polygon) * PUSH_BACK_THRES,1);
    } else if ((result->size & PUSH_BACK_MASK) == 0) {
        memmove((void*)result->polygon, (void*)result->polygon, sizeof(Polygon) * (result->size + PUSH_BACK_THRES));
    }
    result->size++;
}

static void poly_addResultPoly(Clipper_result *result, const Polygon *polygon)
{
    if (result->size == 0) {
        result->polygon = (Polygon*)calloc(sizeof(Polygon) * PUSH_BACK_THRES,1);
    } else if ((result->size & PUSH_BACK_MASK) == 0) {
        memmove((void*)result->polygon, (void*)result->polygon, sizeof(Polygon) * (result->size + PUSH_BACK_THRES));
    }
    result->polygon[result->size++] = *polygon;
}

static int poly_checkIfPointInside(Polygon *poly, const Point *pt)
{
    const Line line = {.from=*pt, .to={99999.9f,pt->y}};
    Line_Ptr edge; 
    int intersection_count = 0;
    for (int i = 0; i < poly->size; ++i) {
        edge = poly_getLinePtr(poly, i);
        intersection_count += line_doIntersect(edge.from, edge.to, &line.from, &line.to);
    }
    if (intersection_count & 0x1)
        return 1;
    return 0;
}

static Clipper_result clip_genericPolygon(const Polygon *in_subj, const Polygon *in_clipper)
{
    Polygon _subj = *in_subj;
    Polygon _clipper = *in_clipper;
    Polygon *subj = &_subj;
    Polygon *clipper = &_clipper;

    Clipper_result result = { };
    Intersection_Point_List intersection_list = { { } };
    Intersection *intersection;
    Polygon *poly;
    // get all edges
    Edge_Ptr_List *subj_edges = poly_getEdgePtrList(subj);
    Edge_Ptr_List *clipper_edges = poly_getEdgePtrList(clipper);

    // get all intersection and orientation for subj polygon
    for (int c = 0; c < clipper->size; ++c) {
        Edge_Ptr *clip_edge = &clipper_edges->edges[c];
        for (int s = 0; s < subj->size; ++s) {
            Edge_Ptr *subj_edge = &subj_edges->edges[s];
            Intersection intersection_point = line_computeIntersection_v3(subj_edge, clip_edge);
            if (intersection_point.type != VERT_UNDERTERMINED) {// && line_isIntersectUnique(&intersection_point, &intersection_list)) {
                intersection = &intersection_list.pts[intersection_list.size++];
                *intersection = intersection_point;
                edge_pushBackIntersection(subj_edge, intersection);
                //for (int ii = 0; ii < subj_edge->size; ii++)
                //    printf("subj subj_edge#(%d) :%d(%p) -> inter %p\n", s, ii, subj_edge, subj_edge->intersection_point[ii]);
                edge_pushBackIntersection(clip_edge, intersection);
            } else continue;
            clip_info_l("C:%d,%p,[%2.f,%2.f][%2.f,%2.f] - S:%d,%p,[%2.f,%2.f][%2.f,%2.f] Inter:(%d)(%p)[%2.2f-%2.2f] L:%p,R:%p\n",
                c,clip_edge,clip_edge->prev->x,clip_edge->prev->y,
                clip_edge->curr->x,clip_edge->curr->y,

                s,subj_edge,subj_edge->prev->x,subj_edge->prev->y,
                subj_edge->curr->x,subj_edge->curr->y,
                intersection_point.type,intersection,intersection_point.pt.x,intersection_point.pt.y,
                intersection_point.edges[0], intersection_point.edges[1]);
        }
    }
    if (intersection_list.size == 0) {
        Point *pt = &clipper->pts[0];
        Edge_Ptr *edge = &subj_edges->edges[0];
        if (line_orientation_v3(edge->prev, edge->curr, pt) == LINE_LEFT) {
            // clipper polygon is inclusive.
            edge_freeList(subj_edges);
            edge_freeList(clipper_edges);
            poly_addResultPoly(&result, clipper);
        }
        return result;
    }

    for (int c = 0; c < clipper_edges->size; ++c) {
        Edge_Ptr *clip_edge = &clipper_edges->edges[c];
        if (clip_edge->size > 1)
            edge_sortIntersection(clip_edge);
        clip_info_l("clip_edge %d:%p-> %p[%.2f,%.2f] %p[%.2f,%.2f]\n", c, clip_edge,
            (clip_edge->size >= 1) ? *clip_edge->intersection_point :0,
            (clip_edge->size >= 1) ? clip_edge->intersection_point[0]->pt.x :0,
            (clip_edge->size >= 1) ? clip_edge->intersection_point[0]->pt.y :0,
            (clip_edge->size >= 2) ? *(clip_edge->intersection_point+1) :0,
            (clip_edge->size >= 2) ? clip_edge->intersection_point[1]->pt.x :0,
            (clip_edge->size >= 2) ? clip_edge->intersection_point[1]->pt.y :0);
    }
    for (int s = 0; s < subj_edges->size; ++s) {
        Edge_Ptr *subj_edge = &subj_edges->edges[s];
        if (subj_edge->size > 1)
            edge_sortIntersection(subj_edge);
        clip_info_l("subj_edge %d:%p-> %p[%.2f,%.2f] %p[%.2f,%.2f]\n", s, subj_edge,
          (subj_edge->size >= 1) ? *subj_edge->intersection_point :0,
          (subj_edge->size >= 1) ? subj_edge->intersection_point[0]->pt.x :0,
          (subj_edge->size >= 1) ? subj_edge->intersection_point[0]->pt.y :0,
          (subj_edge->size >= 2) ? *(subj_edge->intersection_point+1) :0,
          (subj_edge->size >= 2) ? subj_edge->intersection_point[1]->pt.x :0,
          (subj_edge->size >= 2) ? subj_edge->intersection_point[1]->pt.y :0);
    }

    // sort intersection points



    // assign intersection point to clipper polygon
    Edge_Ptr const *edge;
    //int *intersect_marker = (int*)calloc(intersection_list.size * sizeof(int));
    int *vert_in_index = (int*)calloc(intersection_list.size * sizeof(int),1);
    int vert_in_num = 0;

    for (int i = 0; i < intersection_list.size; ++i) {
        if (intersection_list.pts[i].type == VERT_IN)
            vert_in_index[vert_in_num++] = i;
        clip_info_l("inter:%d %p -> edge:%p %p\n", i, &intersection_list.pts[i],
            intersection_list.pts[i].edges[0],
            intersection_list.pts[i].edges[1]);
    }
    
    // if find vert_in in intersection_list // tranverse clipper list
    while (1) {
        int i;
        for (i = 0; i < intersection_list.size; ++i) {
            clip_info_l("%d inter:[%.2f,%.2f] type:%d\n", i,
                intersection_list.pts[i].pt.x,
                intersection_list.pts[i].pt.y,
                intersection_list.pts[i].type);
        }
        for (i = 0; i < vert_in_num; ++i) {
            if (intersection_list.pts[vert_in_index[i]].type == VERT_IN)
                break;
        }
        if (i == vert_in_num)
            break;
        intersection = &intersection_list.pts[vert_in_index[i]];
        // set origin intersection point // checked
        Intersection *orig = intersection;
        Intersection *curr = orig;
        // go to the subj edge
        poly_incrementResultPoly(&result);
        clip_info_l("polygon result size:%d\n", result.size);
        poly = &result.polygon[result.size-1];
        poly_addPoint(poly, &curr->pt);
        curr->type = VERT_CHECKED;
        edge = intersection->edges[1];
        int jump = 1;
        int same = 0;
        Point *curr_pt = NULL;
        do {
            // we can first get the pointer of pt and intersection if any
            clip_info_l("edge size %d.jump:%d\n", edge->size, jump);
            int i = -1;
            if (edge->size == 1) {
               if (!jump) {
                    curr = edge->intersection_point[0];
                    same = curr->type == VERT_SAME;
                    curr_pt = &curr->pt;
                    if (curr!=orig) poly_addPoint(poly, &curr->pt);
                    edge = (curr->edges[0] == edge) ?
                        curr->edges[1] :
                        curr->edges[0];
                    curr->type = VERT_CHECKED;
                    jump = 1;
                } else {
                    if (!same)
                        poly_addPoint(poly, edge->curr);
                    //curr_pt = edge->curr;
                    same = 0;
                    edge = edge->next;
                    curr = NULL;
                    jump = 0;
                }
        //   2 if edge size == 0
            } else if (edge->size == 0) { // ?
                if (!same)
                    poly_addPoint(poly, edge->curr);
                same=0;
                //curr_pt = edge->curr;
                edge = edge->next;
                curr = NULL;
                jump = 0;
        //   3 if edge size > 1
            } else if (edge->size > 1) {
                if (jump) {
                    for (i = 0; i < edge->size; ++i) {
                        if (edge->intersection_point[i] == curr)
                            break;
                    }
                    if (i == (edge->size-1)) {
                        if (!same) {
                            poly_addPoint(poly, edge->curr);
                        }
                        //curr_pt = edge->curr;
                        same=0;
                        curr = NULL;
                        edge = edge->next;
                        jump = 0;
                    } else if (i < edge->size-1) {
                        i++;
                        curr = edge->intersection_point[i];
                        if (curr == orig)
                            break;
                        //curr_pt = &curr->pt;
                        if (curr!=orig) poly_addPoint(poly, &curr->pt);
                        edge = (curr->edges[0] == edge) ?
                            curr->edges[1] :
                            curr->edges[0];
                        same = curr->type == VERT_SAME;
                        curr->type = VERT_CHECKED;
                        jump = 1;
                    } else {
                        assert(0);
                    }
                } else {
                    curr = edge->intersection_point[0];
                    if (curr!=orig) poly_addPoint(poly, &curr->pt);
                    edge = (curr->edges[0] == edge) ?
                            curr->edges[1]:
                            curr->edges[0];
                    same = curr->type == VERT_SAME;
                    curr->type = VERT_CHECKED;
                    jump = 1;
                }
            }
            clip_info_l("poly size:%d (%p,%p) edge size %d.jump:%d,i:%d\n", poly->size, curr, orig, edge->size, jump, i);
        }  while (curr != orig);
        // if curr point is origin intersection
        //   turn all marked intersection complete
        // for next vert in 
        //   go back to while loop 

    //for (int i = 0; i < result.size; ++i) {
    //    poly_removeConsecutiveDuplicate(result.polygon[i]);
    }
    free(vert_in_index);
    return result;
}

#define INCRE_NUM 4
#define INCRE_MASK 3

static Clipper_result clip_convexPolygon(const Polygon *in_subj, const Polygon *in_clipper)
{
    Polygon _subj = *in_subj;
    Polygon _clipper = *in_clipper;
    Polygon *subj = &_subj;
    Polygon *clipper = &_clipper;

    Clipper_result result;
    result.size = 1;
    result.polygon = (Polygon*)malloc(sizeof(Polygon));
    Polygon *output_polygon = result.polygon;

    poly_copy(subj, output_polygon);
    for (int e = 0; e < clipper->size; ++e) {
        Line_Ptr clipEdge = poly_getLinePtr(clipper, e);
        Polygon input_polygon;
        poly_copy(output_polygon, &input_polygon);
        poly_empty(output_polygon);

        for (int i = 0; i < input_polygon.size; ++i) {
            Point *current_point = &input_polygon.pts[i];
            Point *prev_point = &input_polygon.pts[(i + input_polygon.size - 1) % input_polygon.size];            
            Point intersection_point = line_computeIntersection_v2(prev_point, current_point, &clipEdge);

            int current_orient = line_orientation_v2(clipEdge.pts[0], clipEdge.pts[1], current_point);
            int prev_orient = line_orientation_v2(clipEdge.pts[0], clipEdge.pts[1], prev_point);
            /*clip_log("[subj:%d][prev](%.1f, %.1f)[%d] "
                "[cur](%.1f, %.1f)[%d] "
                "vs edge[%d](%.1f, %1.f)(%.1f, %.1f)\n",
                i,prev_point->x, prev_point->y, prev_orient,
                current_point->x, current_point->y, current_orient,
                e, clipEdge.pts[0].x, clipEdge.pts[0].y,
                clipEdge.pts[1].x, clipEdge.pts[1].y);
            */
            if (current_orient == LINE_LEFT) { //https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
                if (prev_orient == LINE_RIGHT) {
                    poly_addPoint(output_polygon, &intersection_point);
                }
                poly_addPoint(output_polygon, current_point);
            } else if (prev_orient == LINE_LEFT) {
                poly_addPoint(output_polygon, &intersection_point);
            }
        }

        #if 0
        for (int i = 0; i < output_polygon->size; ++i) {
            printf("[(%d)%.1f,%1.f],", i, 
                output_polygon->pts[i].x,
                output_polygon->pts[i].y);
        }
        printf("\n");
        #endif
    }
    return result;
}

Clipper_result CLIP_clipPolygon(const Polygon *subj, const Polygon *clipper)
{
    if (subj->type == clipper->type && subj->type == CONVEX ) {
        return clip_convexPolygon(subj, clipper);
    } else if (subj->type == CONCAVE || clipper->type == CONCAVE) {
        return clip_genericPolygon(subj, clipper);
    } else if (subj->type == UNDETERMINED || clipper->type == UNDETERMINED) {
        clip_warn("polygon is not determined yet!");
        return clip_genericPolygon(subj, clipper);
    } else {
        assert(0);
    }
    return clip_genericPolygon(subj, clipper);
}

int CLIP_drawLine(const Point line[2], Array2di* arr, int thickness, int val)
{
    int h = arr->h;
    int w = arr->w;

    int sx = clamp((int)line[0].x, 0, w-1);
    int sy = clamp((int)line[0].y, 0, h-1);
    int ex = clamp((int)line[1].x, 0, w-1);
    int ey = clamp((int)line[1].y, 0, h-1);

    const int dy = (ey - sy);
    const int dx = (ex - sx);

    int slope_fixed = 0;
    if (dy == 0)      { slope_fixed = 0; }
    else if (dx == 0) { slope_fixed = 0x7fffffff; }
    else              { slope_fixed = (dy << FIXED_POINT_SHIFT)
                       /dx ; }

    const int slope_fixed_abs = abs(slope_fixed);

    int x_fixed = sx << FIXED_POINT_SHIFT;
    int y_fixed = sy << FIXED_POINT_SHIFT;
    int dx_fixed = 1, dy_fixed = 1;
    int thickhalf = thickness >> 1;
    int signx = 1, signy = 1;
    if (dy < 0 && dx < 0) {
        signx = -1; signy = -1;
        dy_fixed = -1; dx_fixed = -1;
    }
    else if (dy < 0) { signy = -1; dy_fixed = -1; }
    else if (dx < 0) { signx = -1; dx_fixed = -1; }
    int x, y;
    //printf("xy: %d %d => %d %d dxy: %d %d\n", sx, sy, ex, ey, dx, dy);
    //printf("slope fixed => %d vs one %d \n", slope_fixed, FIXED_POINT_ONE);
    //printf("sign %d %d th(%d)\n", signx, signy,  thickhalf);
    if (slope_fixed_abs >= FIXED_POINT_ONE) {
        dx_fixed *= (FIXED_POINT_ONE << FIXED_POINT_SHIFT) / slope_fixed_abs;
        while (1) {
            x = ROUND_DOWN(x_fixed);
            int start = clamp(x - thickhalf, 0, w);
            int end = clamp(x - thickhalf + thickness, 0, w);
            for (int i = start; i < end; ++i) {
                arr->data[sy*w + i] = val;
            }
            if (sy == ey)
                break;
            x_fixed += dx_fixed;
            sy+=signy;
        }
    } else {
        dy_fixed *= ROUND_DOWN(FIXED_POINT_ONE * slope_fixed_abs);
        while (1) {
            y = ROUND_DOWN(y_fixed);
            int start = clamp(y - thickhalf, 0, h);
            int end = clamp(y - thickhalf + thickness, 0, h);
            for (int i = start; i < end; ++i) {
                arr->data[i*w + sx] = val;
            }
            if (sx == ex)
                break;
            y_fixed += dy_fixed;
            sx+=signx;
        }
    }
    return 0;
}

int CLIP_drawLinePtr(const Point *line[2], Array2di* arr, int thickness, int val)
{
    int h = arr->h;
    int w = arr->w;

    int sx = clamp((int)line[0]->x, 0, w-1);
    int sy = clamp((int)line[0]->y, 0, h-1);
    int ex = clamp((int)line[1]->x, 0, w-1);
    int ey = clamp((int)line[1]->y, 0, h-1);

    const int dy = (ey - sy);
    const int dx = (ex - sx);

    int slope_fixed = 0;
    if (dy == 0)      { slope_fixed = 0; }
    else if (dx == 0) { slope_fixed = 0x7fffffff; }
    else              { slope_fixed = (dy << FIXED_POINT_SHIFT)
                       /dx ; }

    const int slope_fixed_abs = abs(slope_fixed);

    int x_fixed = sx << FIXED_POINT_SHIFT;
    int y_fixed = sy << FIXED_POINT_SHIFT;
    int dx_fixed = 1, dy_fixed = 1;
    int thickhalf = thickness >> 1;
    int signx = 1, signy = 1;
    if (dy < 0 && dx < 0) {
        signx = -1; signy = -1;
        dy_fixed = -1; dx_fixed = -1;
    }
    else if (dy < 0) { signy = -1; dy_fixed = -1; }
    else if (dx < 0) { signx = -1; dx_fixed = -1; }
    int x, y;
    //printf("xy: %d %d => %d %d dxy: %d %d\n", sx, sy, ex, ey, dx, dy);
    //printf("slope fixed => %d vs one %d \n", slope_fixed, FIXED_POINT_ONE);
    //printf("sign %d %d th(%d)\n", signx, signy,  thickhalf);
    if (slope_fixed_abs >= FIXED_POINT_ONE) {
        dx_fixed *= (FIXED_POINT_ONE << FIXED_POINT_SHIFT) / slope_fixed_abs;
        while (1) {
            x = ROUND_DOWN(x_fixed);
            int start = clamp(x - thickhalf, 0, w);
            int end = clamp(x - thickhalf + thickness, 0, w);
            for (int i = start; i < end; ++i) {
                arr->data[sy*w + i] = val;
            }
            if (sy == ey)
                break;
            x_fixed += dx_fixed;
            sy+=signy;
        }
    } else {
        dy_fixed *= ROUND_DOWN(FIXED_POINT_ONE * slope_fixed_abs);
        while (1) {
            y = ROUND_DOWN(y_fixed);
            int start = clamp(y - thickhalf, 0, h);
            int end = clamp(y - thickhalf + thickness, 0, h);
            for (int i = start; i < end; ++i) {
                arr->data[i*w + sx] = val;
            }
            if (sx == ex)
                break;
            y_fixed += dy_fixed;
            sx+=signx;
        }
    }
    return 0;
}

void CLIP_printArr(const Array2di* arr)
{
    const int h = arr->h;
    const int w = arr->w;
    printf(" ");
    for (int j = 0; j < w; ++j)
        printf("%d", j % 10);
    printf("\n");
    for (int i = 0; i < h; ++i) {
        printf("%d", i % 10);
        for (int j = 0; j < w; ++j) {
            int val = arr->data[i*w + j];
            if (val == 1)
                printf("X");
            else if (val == 2)
                printf("O");
            else if (val == 'C')
                printf("C");
            else if (val == 'S')
                printf("S");
            else
                printf("-");
        }
        printf("\n");
    }
}


int CLIP_drawPoint(const Point *pt, Array2di *arr, int val)
{
    arr_fill(arr, pt->x, pt->y, val);
    return 0;
}

int CLIP_drawPoly(const Polygon *polygon, Array2di *arr, int val)
{
    for (int i = 0; i < polygon->size; ++i) {
        Line line = poly_getLine(polygon, i);
        CLIP_drawLine(line.pts, arr, 1, val);
    }
    return 0;
}


static void poly_floodFill(Array2di *arr, int x, int y, int val)
{
    if (arr_get(arr, x, y)) return;
    arr_fill(arr, x, y, val);
    poly_floodFill(arr, x-1, y, val);
    poly_floodFill(arr, x+1, y, val);
    poly_floodFill(arr, x, y-1, val);
    poly_floodFill(arr, x, y+1, val);
}

int CLIP_fillPoly(const Polygon *polygon, Array2di *arr, int val)
{
    const int size = polygon->size;
    Polygon poly = *polygon;
    Line_Ptr *line = (Line_Ptr*)malloc(sizeof(Line_Ptr)*size);
    Point fill_point = {}; 
    const Point *prev = NULL;
    const Point *cur = NULL;
    const Point *next = NULL;

    for (int i = 0; i < size; ++i) {
        line[i] = poly_getLinePtr(&poly, i);
        CLIP_drawLinePtr((const Point**)line[i].pts, arr, 1, val);
    }
    for (int i = 0; i < size; ++i) {
        prev = &poly.pts[(i + size - 1) % size];
        cur = &poly.pts[i];
        next = &poly.pts[(i + 1) % size];
        if (line_orientation(prev, next, cur) == LINE_RIGHT)
            continue;
        else {
            fill_point.x = mean3way(prev->x, cur->x, next->x);
            fill_point.y = mean3way(prev->y, cur->y, next->y);
            int signx = sign(fill_point.x-cur->x);
            int signy = sign(fill_point.y-cur->y);
            fill_point.x = cur->x + signx * 1.5f;
            fill_point.y = cur->y + signy * 1.5f;
            int isinside = poly_checkIfPointInside(&poly, &fill_point);
            clip_trc("3 way[%.1f,%.1f,%.1f]->%.1f,[%.1f,%.1f,%.1f]->%.1f inside?[%d]\n",
                prev->x, cur->x, next->x, fill_point.x,
                prev->y, cur->y, next->y, fill_point.y, isinside);
            if (isinside)
                poly_floodFill(arr, (int)fill_point.x, (int)fill_point.y, val);
            //CLIP_drawPoint(&fill_point, arr, 2);
        }
    }
    
    return 0;   
}

void example1(char **argv)
{
    Array2di arr = arr_init2di(50,35);
    Point pts[2];
    pts[0].x = atof(argv[2]);
    pts[0].y = atof(argv[3]);
    pts[1].x = atof(argv[4]);
    pts[1].y = atof(argv[5]);
    int th = atoi(argv[6]);
    CLIP_drawLine(pts, &arr, th, 1);
    CLIP_printArr(&arr);
    arr_free(arr);
}

void example2(char **argv)
{
    int num = atoi(argv[2]);
    clip_assert(num >= 3, "polygon points should be greater than 3\n");
    Polygon poly = { };
    int idx = 3;
    for (int i = idx; i < idx + (num) * 2; i+=2) {
        poly.pts[poly.size++] = (Point){atof(argv[i]), atof(argv[i+1])};
        printf("(%2.1f  %2.1f) ", poly.pts[poly.size-1].x, poly.pts[poly.size-1].y);
    }
    printf("\n");
    float area = CLIP_getArea(&poly);
    int type = poly_checkType(&poly);
    clip_log("num %d, area is %.4f\n", num, area);
    clip_log("type is %s\n", poly_type_str[type]);

    Array2di arr = arr_init2di(50,25);
    for (int i = 0; i < num; ++i) {
        Line line = poly_getLine(&poly, i);
        CLIP_drawLine(line.pts, &arr, 1, 1);
    }
    CLIP_printArr(&arr);
    arr_free(arr);
}

void example3(char **argv)
{
    int num = atoi(argv[2]);
    clip_assert(num >= 3, "polygon points should be greater than 3\n");
    Polygon poly = { };
    int idx = 3;
    for (int i = idx; i < idx + (num) * 2; i+=2) {
        poly.pts[poly.size++] = (Point){atof(argv[i]), atof(argv[i+1])};
        printf("(%2.1f  %2.1f) ", poly.pts[poly.size-1].x, poly.pts[poly.size-1].y);
    }
    printf("\n");
    float area = CLIP_getArea(&poly);
    int type = poly_checkType(&poly);
    clip_log("num %d, area is %.4f\n", num, area);
    clip_log("type is %s\n", poly_type_str[type]);

    Array2di arr = arr_init2di(50,25);
    CLIP_fillPoly(&poly, &arr, 1);
    CLIP_printArr(&arr);
    arr_free(arr);
}

void example4(char **argv)
{
    // check point inside polygon
    int num = atoi(argv[2]);
    clip_assert(num >= 3, "polygon points should be greater than 3\n");
    Polygon poly = { };
    int idx = 3;
    int i = idx;
    printf("Polygon:\n");
    for (i = idx; i < idx + (num) * 2; i+=2) {
        poly.pts[poly.size++] = (Point){atof(argv[i]), atof(argv[i+1])};
        printf("(%2.1f  %2.1f) ", poly.pts[poly.size-1].x, poly.pts[poly.size-1].y);
    }
    printf("\n");
    Point pt = {atof(argv[i]), atof(argv[i+1])};
    printf("Point:\n");
    printf("(%2.1f  %2.1f) \n", pt.x, pt.y);
    float area = CLIP_getArea(&poly);
    int type = poly_checkType(&poly);
    clip_log("num %d, area is %.4f\n", num, area);
    clip_log("type is %s\n", poly_type_str[type]);

    Array2di arr = arr_init2di(50,25);
    CLIP_drawPoly(&poly, &arr, 1);
    CLIP_drawPoint(&pt, &arr, 2);
    CLIP_printArr(&arr);
    int ifinside = poly_checkIfPointInside(&poly, &pt);
    clip_log("Result is %d\n", ifinside);
}

void example5(char **argv)
{
   Polygon subj = {.size=4,
        .pts = {
            {5, 0},
            {10, 0},
            {10, 20},
            {5, 10}
        }};
    Polygon clipper = {.size=4,
        .pts = {
            {0, 5},
            {20, 5},
            {20, 15},
            {0, 15}
        }};
    subj.type = poly_checkType(&subj);
    clipper.type = poly_checkType(&clipper);
    Array2di arr = arr_init2di(30,25);
    CLIP_drawPoly(&subj, &arr, 1);
    CLIP_drawPoly(&clipper, &arr, 1);
    CLIP_printArr(&arr);
    arr_free(arr);

    Clipper_result result = CLIP_clipPolygon(&subj, &clipper);
    clip_log("subj: %s vs clipper: %s\n",
        poly_type_str[subj.type],
        poly_type_str[clipper.type]);
    Polygon *poly = result.polygon;
    for (int i = 0; i < result.polygon[0].size; ++i) {
        clip_log("(%d: %.2f %.2f)\n", i, poly->pts[i].x,
            poly->pts[i].y);
    }
    clip_log("output area is %.2f\n",CLIP_getArea(poly));
    arr = arr_init2di(30,25);
    CLIP_drawPoly(poly, &arr, 1);
    CLIP_printArr(&arr);
    arr_free(arr);
}


void example6_lambda(const Polygon *in_subj, const Polygon *in_clipper)
{
    Polygon subj = *in_subj;
    Polygon clipper = *in_clipper;

    subj.type = poly_checkType(&subj);
    clipper.type = poly_checkType(&clipper);
    Array2di arr = arr_init2di(60,30);
    CLIP_drawPoly(&subj, &arr, 'S');
    CLIP_drawPoly(&clipper, &arr, 'C');
    CLIP_printArr(&arr);
    arr_free(arr);

    Clipper_result result = CLIP_clipPolygon(&subj, &clipper);
    clip_log("subj: %s vs clipper: %s\n",
        poly_type_str[subj.type],
        poly_type_str[clipper.type]);
    for (int j = 0; j < result.size; ++j) { 
        Polygon *poly = &result.polygon[j];
        for (int i = 0; i < poly->size; ++i) {
            clip_log("res:%d (%d: %.4f %.4f)\n", j, i, poly->pts[i].x,
                poly->pts[i].y);
        }
        clip_log("output area is %.2f\n",CLIP_getArea(poly));
    }
    arr = arr_init2di(60,30);
    for (int j = 0; j < result.size; ++j) { 
        Polygon *poly = &result.polygon[j];
        CLIP_drawPoly(poly, &arr, j+1);
    }
    CLIP_printArr(&arr);
    arr_free(arr);   
}

void example6(char **argv)
{
   clip_log("Example6_1===================\n");
   Polygon subj = {.size=6,
        .pts = {
            {10, 15},
            {30, 2},
            {40, 2},
            {20, 15},
            {40, 28},
            {30, 28}
        }};
    Polygon clipper = {.size=6,
        .pts = {
            {40, 15},
            {20, 28},
            {10, 28},
            {30, 15},
            {10, 2},
            {20, 2},
        }};
    example6_lambda(&subj, &clipper);
    clip_log("Example6_2===================\n");
    subj = (Polygon){.size=4,
        .pts = {
            {15, 5},
            {40, 5},
            {40, 15},
            {15, 15}
        }};
    clipper = (Polygon){.size=6,
        .pts = {
            {40, 15},
            {20, 28},
            {10, 28},
            {30, 15},
            {10, 2},
            {20, 2},
        }};
    example6_lambda(&subj, &clipper);
    clip_log("Example6_3===================\n");
    subj = (Polygon){.size=4,
        .pts = {
            {15, 5},
            {40, 5},
            {40, 15},
            {15, 15}
        }};
    clipper = (Polygon){.size=6,
        .pts = {
            {40, 15},
            {20, 28},
            {10, 28},
            {30, 15},
            {10, 2},
            {20, 2},
        }};
    example6_lambda(&subj, &clipper);
}

const char *helper_str = {
    "\t>>>>>>>>>>  helper str >>>>>>>>>>>>>>>>>\n"
    "\t-l\tdraw line from <sx> <sy> <ex> <ey>\n"
    "\t-a\tcalc and draw polygon outline pt >= 3\n"
    "\t-f\tcalc and fill polygon pt >= 3\n"
    "\t-ch\ttest if pt inside polygon <num> <polygon> <pt>\n"
    "\t-c\tclip convex polygons example\n"
    "\t-ca\tclip concave polygons example\n"
};

int main(int argc, char **argv)
{
    // draw line app
    clip_assert(argc>1, helper_str);

    if (strcmp(argv[1], "-l")==0) {
        clip_assert(argc==7, helper_str);
        example1(argv);
    } else if (strcmp(argv[1], "-a")==0) {
        clip_assert(((argc-2) & 1) == 1 && argc > 8, helper_str);
        example2(argv);
    } else if (strcmp(argv[1], "-f")==0) {
        clip_assert(((argc-2) & 1) == 1 && argc > 8, helper_str);
        example3(argv);
    } else if (strcmp(argv[1], "-ch")==0) {
        example4(argv);
    } else if (strcmp(argv[1], "-c")==0) {
        example5(argv);
    } else if (strcmp(argv[1], "-ca")==0) {
        example6(argv);
    } else {
        clip_log("%s",helper_str);
    }
    return 0;
}
