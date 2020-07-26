#pragma once

const int FIXED_POINT_SHIFT = 8;
const int FIXED_POINT_ONE = (1 << FIXED_POINT_SHIFT);
const int FIXED_POINT_HALF = (1 << (FIXED_POINT_SHIFT-1));

#include <vector>

namespace clip {

using vector2du = std::vector<std::vector<unsigned char>>;
using vector2dc = std::vector<std::vector<signed char>>;
using vector2di = std::vector<std::vector<int>>;
using vector2df = std::vector<std::vector<float>>;

template<class T>
struct Point {
    T x;
    T y;
};

template<class T>
struct Polygon {
    std::vector<Point<T>> pts;
};

template<class T>
struct Triangle {
    Point<T> pts[3];
};

template<class T>
struct Vertex {
    std::vector<Triangle<T>> tris;
};

template<class T>
struct Clipper_Var {
    Polygon<T> polygon;
    Vertex<T> vertex;
};

template<class T>
class Clipper : public Clipper_Var<T> {
public:
    Clipper(void)
    {};
    Clipper(Polygon<T> poly)
    {
        this->polygon = poly;
    }
    const Polygon<T>& get_polygon(void);
    const Vertex<T>& get_vertex(void);
    float get_area(void);
    Clipper clipping(const Polygon<T>& poly);
    float clippingArea(const Polygon<T>& poly);
};

int draw_line(const Point<int> pts[2], vector2di& arr, int thickness);
void print_arr(const vector2di& arr);

template class Clipper<int>;
template class Clipper<float>;
} // clip
