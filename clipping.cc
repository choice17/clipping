#include <algorithm>
#include <cstdio>
#include <cstring>

#include "clipping.hh"

using namespace std;

#define ROUND_DOWN(x) ((x + FIXED_POINT_HALF) >> FIXED_POINT_SHIFT)

inline clamp(int a, int low, int high)
{
    return ( a < low ) ? low : ( ( a > high ) ? high : a );
}
namespace clip {

template<class T>
const Polygon<T>& Clipper<T>::get_polygon(void)
{
    return this->polygon;
}


template<class T>
const Vertex<T>& Clipper<T>::get_vertex(void)
{
    return this->vertex;
}


template<class T>
float Clipper<T>::get_area(void)
{
    int number_vertex = this->polygon.pts.size();
    std::vector<Point<T>>& pts = this->polygon.pts;
    T area = 0;
    int next = 0;
    for (int i = 0; i<number_vertex; ++i) {
        next = (i + 1) % number_vertex;
        area += (pts[i].x * (pts[next].y)) - ((pts[i].y) * pts[next].x);
    }
    return abs((float)area/2.f);
}

int draw_line(const Point<int> pts[2], vector2di& arr, int thickness)
{
    int h = arr.size();
    int w = arr[0].size();

    int sx = clamp(pts[0].x, 0, w-1);
    int sy = clamp(pts[0].y, 0, h-1);
    int ex = clamp(pts[1].x, 0, w-1);
    int ey = clamp(pts[1].y, 0, h-1);

    const int dy = (ey - sy + 1);
    const int dx = (ex - sx + 1);

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
        while (sy != ey) {
            x_fixed += dx_fixed;
            x = ROUND_DOWN(x_fixed);
            int start = clamp(x - thickhalf, 0, w);
            int end = clamp(x - thickhalf + thickness, 0, w);
            for (int i = start; i < end; ++i) {
                arr[sy][i] = 1;
            }
            sy+=signy;
        }
    } else {
        dy_fixed *= ROUND_DOWN(FIXED_POINT_ONE * slope_fixed_abs);
        while (sx != ex) {
            y_fixed += dy_fixed;
            y = ROUND_DOWN(y_fixed);
            int start = clamp(y - thickhalf, 0, h);
            int end = clamp(y - thickhalf + thickness, 0, h);
            for (int i = start; i < end; ++i) {
                arr[i][sx] = 1;
            }
            sx+=signx;
        }
    }
    return 0;
}

void print_arr(const vector2di& arr)
{
    int h = arr.size();
    int w = arr[0].size();
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            if (arr[i][j] != 0)
                printf("X");
            else
                printf("-");
        }
        printf("\n");
    }
}

} //clip

int main(int argc, char **argv)
{
    // draw line app
    if (strcmp(argv[1], "-l")==0) {
        clip::vector2di arr(35, vector<int>(50,0));
        clip::Point<int> pts[2];
        pts[0].x = atoi(argv[2]);
        pts[0].y = atoi(argv[3]);
        pts[1].x = atoi(argv[4]);
        pts[1].y = atoi(argv[5]);
        int th = atoi(argv[6]);
        clip::draw_line(pts, arr, th);
        clip::print_arr(arr);
    } else if (strcmp(argv[1], "-a")==0) {
        int num = atoi(argv[2]);
        clip::Polygon<int> poly;
        int idx = 3;
        for (int i = idx; i < (idx*2+num); i+=2) {
            clip::Point<int> pt{atoi(argv[i]), atoi(argv[i+1])};
            poly.pts.push_back(pt); 
        }
        printf("%d %d %d %d %d %d\n", poly.pts[0].x,poly.pts[0].y,poly.pts[1].x,poly.pts[1].y, poly.pts[2].x,poly.pts[2].y);
        clip::Clipper<int> clips(poly);
        float area = clips.get_area();
        printf("area is %.4f\n", area);
    }
    return 0;
}
