/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#ifndef STL_UTILS_H
#define STL_UTILS_H

#include <cmath>
#include "matar.h"

using namespace mtr;


// a vector type with 3 components
struct vec_t{
    double x;
    double y;
    double z;
    
    // default constructor
    KOKKOS_INLINE_FUNCTION
    vec_t (){};
    
    // overloaded constructor
    KOKKOS_INLINE_FUNCTION
    vec_t(const double x_in, const double y_in, const double z_in){
        x = x_in;
        y = y_in;
        z = z_in;
    };

    // overload opperators
    KOKKOS_INLINE_FUNCTION
    vec_t operator+(const vec_t& b) const { return {x + b.x, y + b.y, z + b.z}; }

    KOKKOS_INLINE_FUNCTION
    vec_t operator-(const vec_t& b) const { return {x - b.x, y - b.y, z - b.z}; }

    KOKKOS_INLINE_FUNCTION
    vec_t operator*(double s) const { return {x * s, y * s, z * s}; }
    
}; // end vec_t


// a triangle data type
struct triangle_t {
    
    vec_t normal; // surface normal
    
    vec_t p[3];   // three nodes with x,y,z coords
    
    // default constructor
    KOKKOS_INLINE_FUNCTION
    triangle_t(){};

    // overloaded constructor to accept array having 3 vectors
    KOKKOS_INLINE_FUNCTION
    triangle_t (const vec_t p_in[3])
    {
        // store the coords
        p[0]=p_in[0]; 
        p[1]=p_in[1]; 
        p[2]=p_in[2];

        // calculate the normal to this surface

        //A = p1 - p0;
        //B = p2 - p0;
        vec_t A;
        A.x = p[1].x - p[0].x;
        A.y = p[1].y - p[0].y;
        A.z = p[1].z - p[0].z;
        
        vec_t B;
        B.x = p[2].x - p[0].x;
        B.y = p[2].y - p[0].y;
        B.z = p[2].z - p[0].z;
        
        // normal = cross product of A and B
        normal.x = A.y * B.z - A.z * B.y;
        normal.y = A.z * B.x - A.x * B.z;
        normal.z = A.x * B.y - A.y * B.x;
        
        const double mag = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
        
        // save the unit normal
        normal.x /= mag;
        normal.y /= mag;
        normal.z /= mag;
    };
    
}; // end triangle_t


/////////////////////////////////////////////////////////////////////////////
///
/// \fn cross product
///
/// \brief returns the cross product of two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
vec_t cross(const vec_t &a, const vec_t &b) {
    return {a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x};
} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn dot product
///
/// \brief the dot product of two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double dot(const vec_t &a, const vec_t &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
} // end function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn magnitude
///
/// \brief the magnitude of a vector
///
/// \param vec_t the vector having coordinates in x, y, z space
///
///////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double magnitude(const vec_t &a){
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
} // end function


////////////////////////////////////////////////////////////////////////////
///
/// \fn distance
///
/// \brief returns the magnitude of a the difference tween two vectors
///
/// \param vec_t the first vector having coordinates in x, y, z space
/// \param vec_t the second vector having coordinates in x, y, z space
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double distance(const vec_t &a, const vec_t &b){
    return sqrt((a.x-b.x)*(a.x-b.x) + 
                (a.y-b.y)*(a.y-b.y) + 
                (a.z-b.z)*(a.z-b.z));
} // end function


////////////////////////////////////////////////////////////////////////////
///
/// \fn normalize
///
/// \brief returns a unit vector
///
/// \param vec_t the vector having coordinates in x, y, z space
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
vec_t normalize(const vec_t& v) {
    double len = magnitude(v);
    if (len == 0.0) return {0,0,0};
    return v * (1.0 / len);
}


////////////////////////////////////////////////////////////////////////////
///
/// \fn clamp
///
/// \brief returns a value clamped between low and high values
///
/// \param double the value passed into the function to be clamped
/// \param double the low value
/// \param double the high value
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double clamp(double v, double lo, double hi) {
    return fmax(lo, fmin(v, hi));
} // end function


////////////////////////////////////////////////////////////////////////////
///
/// \fn Closest point on segment ab to point p
///
/// \brief returns the closest point on line segment between point a and b 
///        to point p
///
/// \param double the point
/// \param double the starting point on the line sement
/// \param double the ending point on the line sement
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
vec_t closest_point_on_segment(const vec_t& p, 
                               const vec_t& a, 
                               const vec_t& b) {
    vec_t ab = b - a;
    double t = dot(p - a, ab) / dot(ab, ab);
    t = clamp(t, 0.0, 1.0);
    return a + ab * t;
} // end function


////////////////////////////////////////////////////////////////////////////
///
/// \fn Closest point on triangle abc to point p (Ericson method)
///
/// \brief Returns the closest point on a trianglular facet (defined by 
///        points a, b, and c) to point p 
///
/// \param double the point
/// \param double the first point defining the triangle
/// \param double the second point defining the triangle
/// \param double the third point defining the triangle
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
vec_t closest_point_on_triangle(const vec_t& p, 
                                const vec_t& a, 
                                const vec_t& b, 
                                const vec_t& c) {
    vec_t ab = b - a;
    vec_t ac = c - a;
    vec_t ap = p - a;

    double d1 = dot(ab, ap);
    double d2 = dot(ac, ap);

    if (d1 <= 0.0 && d2 <= 0.0) return a;

    vec_t bp = p - b;
    double d3 = dot(ab, bp);
    double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        return closest_point_on_segment(p, a, b);

    vec_t cp = p - c;
    double d5 = dot(ab, cp);
    double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        return closest_point_on_segment(p, a, c);

    double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        return closest_point_on_segment(p, b, c);

    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;

    return a + ab * v + ac * w;
} // end function


////////////////////////////////////////////////////////////////////////////
///
/// \fn Signed distance from point p to triangle abc
///
/// \brief Returns the signed distance value to a trianglular facet (defined 
///        by points a, b, and c) from point p 
///
/// \param double the point
/// \param double the first point defining the triangle
/// \param double the second point defining the triangle
/// \param double the third point defining the triangle
///
////////////////////////////////////////////////////////////////////////////// 
KOKKOS_INLINE_FUNCTION
double signed_distance_to_triangle(const vec_t& p,
                                   const vec_t& a,
                                   const vec_t& b,
                                   const vec_t& c) {

    vec_t closest = closest_point_on_triangle(p, a, b, c);

    vec_t ab = b - a;
    vec_t ac = c - a;
    vec_t normal = normalize(cross(ab, ac));

    vec_t diff = p - closest;

    double dist = magnitude(diff);

    // sign from orientation
    double sign = (dot(diff, normal) >= 0.0) ? 1.0 : -1.0;

    return sign * dist;
} // end function




#endif