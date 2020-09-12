#ifndef CLIPPER_H_
#define CLIPPER_H_

#define MAX_POINTS 20

typedef struct {
    int w;
    int h;
    int size;
    int *data;
} Array2di;

typedef struct {
    float x;
    float y;
} Point;

typedef enum {UNDETERMINED, CONVEX, CONCAVE} Polygon_Type;

typedef struct {
    Point pts[MAX_POINTS];
    int size;
    int type;
} Polygon;

typedef union {
    Point pts[2];
    struct {
        Point from;
        Point to;
    };
} Line;

typedef struct {
    Point pts[3];
} Triangle;

typedef struct {
    Triangle tris[MAX_POINTS];                                                 
    int size;
} Vertex;

typedef struct {
    Vertex vertex;
    Polygon polygon;
    Point edges[MAX_POINTS][2];
    int convex;
} Clipper;

typedef struct {
    Polygon *polygon;
    int size;
} Clipper_result;

int CLIP_checkIsConvex(const Polygon *polygon);
float CLIP_getArea(const Polygon *polygon);
int CLIP_drawPoly(const Polygon *polygon, Array2di* arr, int val);
int CLIP_fillPoly(const Polygon *polygon, Array2di* arr, int val);
int CLIP_drawLine(const Point line[2], Array2di* arr, int thickness, int val);
void CLIP_printArr(const Array2di* arr);
Clipper_result CLIP_clipPolygon(const Polygon *subj, const Polygon *clipper);

#endif