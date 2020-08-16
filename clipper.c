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

const char *poly_type_str[] = {
    "Undetermined",
    "Convex",
    "Concave"
};

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
    LINE_LEFT = 1,
    LINE_RIGHT = 2,
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
    Point intersecting_point = { };
    const Point *pt3 = edge;
    const Point *pt4 = edge + 1;

    float x1y2_y1x2 = pt1->x * pt2->y - pt1->y * pt2->x;
    float x3_x4 = pt3->x - pt4->x;
    float x1_x2 = pt1->x - pt2->x;
    float x3y4_y3x4 = pt3->x * pt4->y - pt3->y * pt4->x;
    float y3_y4 = pt3->y - pt4->y;
    float y1_y2 = pt1->y - pt2->y;

    float denom = x1_x2 * y3_y4 - y1_y2 * x3_x4;

    if (denom == 0.0f ) return intersecting_point;

    intersecting_point.x =
    (x1y2_y1x2 * x3_x4 - x1_x2 * x3y4_y3x4) /
    denom;
    intersecting_point.y =
    (x1y2_y1x2 * y3_y4 - y1_y2 * x3y4_y3x4) /
    denom;

    return intersecting_point;
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

static Line poly_getEdge(const Polygon *polygon, int idx)
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

static inline void poly_empty(Polygon *polygon)
{
    polygon->size = 0;
    polygon->type = 0;
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
    polygon->pts[polygon->size-1] = *pt;
    polygon->size++;
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

float CLIP_getArea(const Polygon *polygon)
{
    const int number_vertex = polygon->size;
    const Point *pts = polygon->pts;
    float area = 0;
    int next = 0;
    for (int i = 0; i<number_vertex; ++i) {
        next = (i + 1) % number_vertex;
        area += (pts[i].x * (pts[next].y)) - ((pts[i].y) * pts[next].x);
        printf("%d, area:%.5f\n", i, area);
    }
    printf("\n");
    return abs((float)area/2.f);
}

static Clipper_result clip_concavePolygon(const Polygon *subj, const Polygon *clipper)
{
    Clipper_result result = { };
    return result;
}

static Clipper_result clip_convexPolygon(const Polygon *subj, const Polygon *clipper)
{
    Clipper_result result;
    result.size = 1;
    result.polygon = (Polygon*)malloc(sizeof(Polygon));
    Polygon *output_polygon = result.polygon;

    poly_copy(subj, output_polygon);
    for (int e = 0; e < clipper->size; ++e) {
        Line clipEdge = poly_getEdge(clipper, e);
        Polygon input_polygon;
        poly_copy(&input_polygon, output_polygon);
        poly_empty(output_polygon);

        for (int i = 0; i < input_polygon.size; ++i) {
            Point *current_point = &input_polygon.pts[i];
            Point *prev_point = &input_polygon.pts[(i + input_polygon.size - 1) % input_polygon.size];            
            Point intersection_point = line_computeIntersection(prev_point, current_point, clipEdge.pts);

            if (line_orientation(&clipEdge.pts[0], &clipEdge.pts[1], current_point) == LINE_LEFT) { //https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
                if (line_orientation(&clipEdge.pts[0], &clipEdge.pts[1], current_point) == LINE_RIGHT) {
                    poly_addPoint(output_polygon, &intersection_point);
                }
                poly_addPoint(output_polygon, current_point);
            } else if (line_orientation(&clipEdge.pts[0], &clipEdge.pts[1], prev_point) == LINE_LEFT) {
                poly_addPoint(output_polygon, &intersection_point);
            }
        }
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

    }
    return 0;
}
