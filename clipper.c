#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "clipper.h"

#define info_fmt(x)   "[INFO] " x
#define warn_fmt(x)   "[WARN] " x
#define erro_fmt(x)   "[ERROR] " x
#define trc_fmt(x) "[INFO] %s:%d " x
#define clip_trc(x, args...) printf(trc_fmt(x), __func__, __LINE__, ##args)
#define clip_log(x, args...) printf(info_fmt(x), ##args)
#define clip_warn(x, args...) printf(warn_fmt(x), ##args)
#define clip_info_h
#define clip_info_l 

#define FIXED_POINT_SHIFT (8)
#define FIXED_POINT_ONE (1 << FIXED_POINT_SHIFT)
#define FIXED_POINT_HALF (1 << (FIXED_POINT_SHIFT-1))

#define ROUND_DOWN(x) ((x + FIXED_POINT_HALF) >> FIXED_POINT_SHIFT)

#define poly_free(x) (free(x))
#define PUSH_BACK_THRES (8)
#define PUSH_BACK_MASK (7)

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
} Line_Ptr;

typedef struct Edge_Ptr {
    Point *prev;
    Point *curr;
    Intersection **intersection_point;
    Edge_Ptr *next;
    int size;
} Edge_Ptr;

typedef struct {
    Edge_Ptr edges[MAX_POINTS];
    int size; 
} Edge_Ptr_List;

typedef struct {
    int type;
    Edge *edges[2];
    Point pt;
} Intersection;

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
}

void edge_pushBackIntersection(Edge_Ptr *edge, Intersection *pt)
{
    if (edge->size == 0)
        edge->intersection_point = (**Intersection)calloc(sizeof(*Intersection_point)*PUSH_BACK_THRES);
    else if ((edge->size & PUSH_BACK_MASK) == 0)
        memmove(edge->intersection_point, sizeof(*Intersection_point)*PUSH_BACK_THRES);
    edge->intersection_point[edge->size++] = pt;
}

float max(float a, float b)
{
    return ( a > b ) ? a : b;
}

float min(float a, float b)
{
    return ( a < b ) ? a : b;
}

int clamp(int a, int low, int high)
{
    return ( a < low ) ? low : ( ( a > high ) ? high : a );
}

static int line_onSegment(const Point* p, const Point* q, const Point* r) 
{ 
    if (q->x <= max(p->x, r->x) && q->x >= min(p->x, r->x) && 
            q->y <= max(p->y, r->y) && q->y >= min(p->y, r->y)) 
        return 1; 
    return 0; 
} 

typedef enum {
    LINE_SAME = 0;
    LINE_LEFT = 1,
    LINE_RIGHT = 2,
    LINE_IN = 1,
    LINE_OUT = 2
} LINE_DIRECTION;

typedef enum {
    LINE_CLOCKWISE = 1,
    LINE_COUNTER = 2,
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


static Intersection line_computeIntersection_v3(const Edge_Ptr *edge, const Edge_Ptr *query_edge)
{
    Intersection intersection_point = {UNDETERMINED, {NULL, NULL}, {-1, -1} };
    const Point *pt1 = edge->prev;
    const Point *pt2 = edge->curr;
    const Point *pt3 = query_edge->prev;
    const Point *pt4 = query_edge->curr;

    // filter out of box
    if (pt4->x < pt1->x || pt4->y < pt1->y ||
        pt3->x > pt2->x || pt3->y > pt2->y)
        return intersection_point;

    int pt3_orig = line_orientation_v3(pt1, pt2, pt3);
    int pt4_orig = line_orientation_v3(pt1, pt2, pt4);


    if (pt4_orig == LINE_SAME) {
        intersection_point.type = VERT_SAME;
        intersection_point.pt = *pt4;
        intersection_point.edges[0] = edge;
        intersection_point.edges[1] = query_edge;
        return intersection_point;
    }

    // filter no intersection
    if (pt3_orig == pt4_orig)
        return intersection_point;

    float x1y2_y1x2 = pt1->x * pt2->y - pt1->y * pt2->x;
    float x3_x4 = pt3->x - pt4->x;
    float x1_x2 = pt1->x - pt2->x;
    float x3y4_y3x4 = pt3->x * pt4->y - pt3->y * pt4->x;
    float y3_y4 = pt3->y - pt4->y;
    float y1_y2 = pt1->y - pt2->y;

    float denom = x1_x2 * y3_y4 - y1_y2 * x3_x4;

    if (denom == 0.0f ) return intersection_point;

    intersection_point.type = (pt4_orig == LINE_LEFT) ? VERT_IN : VERT_OUT;

    intersection_point.pt.x =
    (x1y2_y1x2 * x3_x4 - x1_x2 * x3y4_y3x4) /
    denom;
    intersection_point.pt.y =
    (x1y2_y1x2 * y3_y4 - y1_y2 * x3y4_y3x4) /
    denom;

    return intersection_point;
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
    Edge_Ptr_List *edge_list = malloc(sizeof(Edge_Ptr_List));
    edge_list->size = polygon->size;
    for (int i = 0; i < polygon->size; ++i) {
        edge_list->edges[i] = poly_getEdgePtr(polygon, i);
        int next_idx = (i + polygon->size - 1) % polygon->size;
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

int compare_func_hori_acc(void *a, void *b)
{
    return ((Intersection*)a)->x > ((Intersection*)b)->x;
}

int compare_func_hori_dec(void *a, void *b)
{
    return ((Intersection*)a)->x < ((Intersection*)b)->x;
}

int compare_func_vert_acc(void *a, void *b)
{
    return ((Intersection*)a)->y > ((Intersection*)b)->y;
}

int compare_func_vert_dec(void *a, void *b)
{
    return ((Intersection*)a)->y < ((Intersection*)b)->y;
}

static void edge_sortIntersection(Edge_Ptr *edge)
{
    int (*compare_func)(const void*, const void*);

    if (edge->prev.x == edge->curr.x) {
        if (edge->prev.y < edge->curr.y) {
            compare_func = compare_func_vert_acc;
        } else {
            compare_func = compare_func_vert_dec;
        }
    } else {
        if (edge->prev.x < edge->curr.x) {
            compare_func = compare_func_hori_acc;
        } else {
            compare_func = compare_func_hori_dec;
        }
    }
    qsort(edge->intersection_point, edge->size, sizeof(Intersection),compare_func);
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

static void poly_incrementResultPoly(Clipper_result *result)
{
    if (result->size == 0) {
        result->polygon = (Polygon*)calloc(sizeof(Polygon) * INCRE_NUM);
    } else if (result->size & INCRE_MASK == 0) {
        memmove(result->polygon, sizeof(Polygon) * (result->size + INCRE_MASK));
    }
    result->size++;
}

static void poly_addResultPoly(Clipper_result *result, const Polygon *polygon)
{
    if (result->size == 0) {
        result->polygon = (Polygon*)calloc(sizeof(Polygon) * INCRE_NUM);
    } else if (result->size & INCRE_MASK == 0) {
        memmove(result->polygon, sizeof(Polygon) * (result->size + INCRE_MASK));
    }
    result->polygon[result->size++] = *polygon;
}


static Clipper_result clip_genericPolygon(const Polygon *subj, const Polygon *clipper)
{
    Clipper_result result = { };
    Vert_Ptr_List subj_list, clipper_list;
    Intersection_point_list intersection_list = { };
    Intersection *intersection;
    Point *pt;
    Polygon *poly;
    // get all edges
    Edge_Ptr_List subj_edges = poly_getEdgePtrList(subj);
    Edge_Ptr_List clipper_edges = poly_getEdgePtrList(clipper);

    // get all intersection and orientation for subj polygon
    for (int c = 0; c < clipper->size; ++e) {
        Edge_Ptr *clip_edge = &clipper->edges[c];
        for (int s = 0; s < subj->size; ++s) {
            Edge_Ptr *subj_edge = &subj_edges->edges[s];
            Intersection intersection_point = line_computeIntersection_v3(subjEdge, clipEdge);
            if (intersection_point.type != UNDETERMINED) {
                intersection = &intersection_list.pts[intersection_list.size++];
                *intersection = intersection_point;
                edge_pushBackIntersection(subjEdge, pt);
                edge_pushBackIntersection(clipEdge, pt);
            }
        }
    }
    if (intersection_list.size == 0) {
        Point *pt = &clipper->pts[0];
        Edge_Ptr *edge = &subj_edges->edges[0];
        if (line_orientation_v3(edge->prev, edge->curr, pt) == LINE_LEFT) {
            // clipper polygon is inclusive.
            edge_freeList(&subj_edges);
            edge_freeList(&poly_edges);
            poly_addResultPoly(&result, clipper);
        }
        return result;
    }

    // sort intersection points
    for (int c = 0; c < clipper_edges->size; ++e) {
        Edge_Ptr clip_edge = clipper_edges->edges[c];
        if (clip_edge.size >= 2)
            edge_sortIntersection(clip_edge);
    }
    for (int s = 0; s < subj_edges->size; ++e) {
        Edge_Ptr subj_edge = subj_edges->edges[c];
        if (subj_edge.size >= 2)
            edge_sortIntersection(subj_edge);
    }

    // assign intersection point to clipper polygon
    int edge_index = 0;
    Edge_Ptr *edge;
    int *intersect_marker = (int*)calloc(intersection_list.size * sizeof(int));
    int *vert_in_index = (int*)calloc(intersection_list.size * sizeof(int));
    int vert_in_num = 0;

    for (int i = 0; i < intersection_list.size; ++i) {
        if (intersection_list.pts[i].type == VERT_IN)
            vert_in_index[vert_in_num++] = i;
    }

    // get first intersection point vert in
    // set origin intersection point // checked
    // go to the subj edge
    // while curr point is not origin intersection
    //   1. if edge size == 1
    //      if not jump
    //      a. add intersection point // checked
    //         curr pt = intersection point
    //         jump to opposite edge
    //         jump = 1
    //      if jump:
    //         add edge->curr point 
    //         curr pt = edge->curr
    //         edge = edge->ext
    //         jump = 0
    //   2 if edge size == 0
    //      a. add edge->curr point 
    //         curr pt == edge->curr
    //         edge = edge->next
    //         jump = 0
    //   3 if edge size > 1
    //     if jump
    //     for i++;;
    //         if edge.intersection[i].pt == curr
    //            if i == edge->size-1
    //                add edge->curr point
    //                curr pt = edge->curr
    //                edge = edge->next
    //                jump = 0
    //            else
    //                add intersection point // checked
    //                curr pt = intersection point
    //                jump to opposite edge
    //                jump = 1
    //     if not jump       
    //         add intersection point // checked
    //         curr pt = intersection point
    //         jump to opposite edge
    //         jump = 1
    // if curr point is origin intersection
    //   turn all marked intersection complete
    // for next vert in 
    //   go back to while loop 
    intersection = &intersection_list.pts[vert_in_index[0]];
    intersect_marker[vert_in_index[0]] = 1; // marked
    poly_incrementResultPoly(&result); // incre
    poly = result.polygon[result->size-1]; // get polygon for result
    poly_addPoint(poly, &intersection->pt); // add intersection pt;
    edge = intersection->edge[1]; // get clipper edge;
    while (1) {
        if (edge->size == 0) {
            if (start) {
                // copy to polygon
                poly = result.polygon[result->size-1];
                poly_addPoint(poly, edge->curr);
            } else 
                continue;
        } else if (edge->size == 1) {
            poly_addPoint(poly, edge->curr);
        } else {
            for (int i = 0; i < edge->size; ++i) {
                intersection = edge->intersection_point[i];
                if (!start++) {
                    if (intersection->type == VERT_IN) {
                        // start polygon
                        poly_incrementResultPoly(&result);
                        poly = result.polygon[result->size-1];
                        poly_addPoint(poly, &intersection->pt);
                        intersection->type = VERT_CHECKED;
                    }
                    continue;
                } else if (start) {
                    poly = result.polygon[result->size-1];
                    if (intersect->type == VERT_OUT) {
                        poly_addPoint(poly, &intersection->pt);
                        edge = edge == intersection->edges[0] ? intersection->edges[1] : intersection->edges[0];  
                        edge = edge->next;
                        intersection->type = VERT_CHECKED;
                        break;
                    } else if (intersection->type == VERT_SAME) {
                        poly_addPoint(poly, &intersection->pt);
                        intersection->type = VERT_CHECKED;
                    } else if (intersection->type == VERT_CHECKED) {
                        intersection->type = VERT_COMPLETE;
                        start = 0;
                    } else
                       continue;
                }
            }
            poly = result.polygon[result->size-1];
            poly_addPoint(poly, edge->curr);
        }
        edge = edge.next;
    return result;
}

#define INCRE_NUM 4
#define INCRE_MASK 3

static Clipper_result clip_convexConvexPolygon(const Polygon *subj, const Polygon *clipper)
{
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
    if (subj->type == CONVEX) {
        return clip_convexPolygon(subj, clipper);
    } else if (subj->type == CONCAVE) {
        return clip_concavePolygon(subj, clipper);
    } else if (subj->type == UNDETERMINED) {
        clip_warn("polygon is not determined yet!");
        return clip_convexPolygon(subj, clipper);
    } else {
        assert(0);
    }
    return clip_convexPolygon(subj, clipper);
}

int CLIP_drawLine(const Point line[2], Array2di* arr, int thickness)
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
                arr->data[sy*w + i] = 1;
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
                arr->data[i*w + sx] = 1;
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
            if (arr->data[i*w + j] != 0)
                printf("X");
            else
                printf("-");
        }
        printf("\n");
    }
}


int CLIP_drawPoly(const Polygon *polygon, Array2di *arr)
{
    for (int i = 0; i < polygon->size; ++i) {
        Line line = poly_getEdge(polygon, i);
        CLIP_drawLine(line.pts, arr, 1);
    }
    return 0;
}
int main(int argc, char **argv)
{
    // draw line app
    if (strcmp(argv[1], "-l")==0) {
        Array2di arr = arr_init2di(50,35);
        Point pts[2];
        pts[0].x = atof(argv[2]);
        pts[0].y = atof(argv[3]);
        pts[1].x = atof(argv[4]);
        pts[1].y = atof(argv[5]);
        int th = atoi(argv[6]);
        CLIP_drawLine(pts, &arr, th);
        CLIP_printArr(&arr);
        arr_free(arr);
    } else if (strcmp(argv[1], "-a")==0) {
        int num = atoi(argv[2]);
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

        Array2di arr = arr_init2di(30,20);
        for (int i = 0; i < num; ++i) {
            Line line = poly_getEdge(&poly, i);
            CLIP_drawLine(line.pts, &arr, 1);
        }
        CLIP_printArr(&arr);
        arr_free(arr);
    } else if (strcmp(argv[1], "-c")==0) {
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
        CLIP_drawPoly(&subj, &arr);
        CLIP_drawPoly(&clipper, &arr);
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
        CLIP_drawPoly(poly, &arr);
        CLIP_printArr(&arr);
        arr_free(arr);

    }
    return 0;
}
