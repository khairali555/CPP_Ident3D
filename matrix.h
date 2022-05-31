#ifndef MATRIX_H
#define MATRIX_H

#include <math.h>
#include <string.h>
//MFC #include <afxtempl.h>
#include "PI.h"

class R2Vector;
class R2Point;
class R3Vector;

// Fight againts Microsoft
#include "MFC.h"        

class R2Matrix;
class R3Matrix;

class ByteR3Vector;

class R3Box;

class McDrawable;
class McGC;

const double R2_EPSILON = 0.0000001;
const double R3_EPSILON = 0.0000001;

const int NormalLength = 120;
const double RNormalLength = double(NormalLength);
const double InvRNormalLength = 1. / RNormalLength;

class R2Vector
{
public:
    double x, y;

public:
    R2Vector() : x( 0. ), y( 0. ) {}
    R2Vector( double X, double Y ) : x( X ), y( Y ) {}
    R2Vector( const R2Vector& v ) : x( v.x ), y( v.y ) {}
    ~R2Vector() {}

    const R2Vector& operator=( const R2Vector& v )
    {
        x = v.x; y = v.y; return *this;
    }

    R2Vector operator+( const R2Vector& v ) const
    {
        return R2Vector( x + v.x, y + v.y );
    }

    R2Vector operator-( const R2Vector& v ) const
    {
        return R2Vector( x - v.x, y - v.y );
    }

    const R2Vector& operator+=( const R2Vector& v )
    {
        x += v.x; y += v.y; return *this;
    }

    const R2Vector& operator-=( const R2Vector& v )
    {
        x -= v.x; y -= v.y; return *this;
    }

    R2Vector operator*( double r ) const
    {
         return R2Vector( r * x,  r * y );
    }

    const R2Vector& operator*=( double r )
    {
        x *= r;  y *= r; return *this;
    }

    double operator*( const R2Vector& v ) const
    {
        return  x * v.x + y * v.y;
    }

    double length2() const { return x*x + y*y; }
    double length() const { return sqrt(x*x + y*y); }
    double Length() const { return length(); }

    R2Vector& operator*=(const R2Matrix&);

    R2Vector Normal() const
    {
        return R2Vector(-y, x);
    }
    R2Vector normal() const { return Normal(); }

    // Make length = 1.
    R2Vector& Normalize()
    {
        double l = Length();
        if (l > R2_EPSILON)
        {
            x /= l;
            y /= l;
        }
        return *this;
    }
    R2Vector& normalize() { return Normalize(); }

    bool Parallel(const R2Vector& v) const
    {
        return (fabs(x*v.y - y*v.x) <= R2_EPSILON);
    }
    bool parallel(const R2Vector& v) const { return Parallel(v); }

    double Angle(const R2Vector& v) const;
    double angle(const R2Vector& v) const { return Angle(v); }

    // Comparings
    bool operator==(const R2Vector& v) const {
        //... return (x == v.x && y == v.y);
        return (
            fabs(x - v.x) <= R2_EPSILON &&
            fabs(y - v.y) <= R2_EPSILON
        );
    }
    bool closeTo(const R2Vector& v, double eps = R2_EPSILON) const {
        return (
            fabs(x - v.x) <= eps &&
            fabs(y - v.y) <= eps
        );
    }
    bool operator!=(const R2Vector& v) const { return !operator==(v); }
    bool operator>=(const R2Vector& v) const {
        //... return (x > v.x || (x == v.x && y >= v.y));
        return (x > v.x || (x >= v.x && y >= v.y));
    }
    bool operator>(const R2Vector& v) const {
        //... return (x > v.x || (x == v.x && y > v.y));
        return (x > v.x || (x >= v.x && y > v.y));
    }
    bool operator<(const R2Vector& v) const { return !operator>=(v); }
    bool operator<=(const R2Vector& v) const { return !operator>(v); }

    // Area of oriented parallelogram (determinant)
    double signed_area(const R2Vector& v) const {
        return (x * v.y - y * v.x);
    }

    static double signed_area(
        const R2Vector& a, const R2Vector& b
    ) {
        return a.signed_area(b);
    }
};

inline R2Vector operator*(double c, const R2Vector& v) {
    return R2Vector(c*v.x, c*v.y);
}

class R2Point {
public:
    double x;
    double y;

    R2Point():                         // Default constructor
        x(0.),
        y(0.)
    {}

    R2Point(const R2Point& p):        // Copy-constructor
        x(p.x),
        y(p.y)
    {}

    R2Point(double xx, double yy):
        x(xx),
        y(yy)
    {}

    R2Point& operator=(const R2Point& p) {    // Copy-operator
        x = p.x; y = p.y;
        return *this;
    }

    ~R2Point() {}                              // Destructor

    R2Point operator+(const R2Point& p) const {
        return R2Point(x+p.x, y+p.y);
    }

    R2Point operator+(const R2Vector& v) const {
        return R2Point(x+v.x, y+v.y);
    }

    R2Point& operator+=(const R2Point& p) {
        x += p.x;
        y += p.y;
        return *this;
    }

    R2Point& operator+=(const R2Vector& v) {
        x += v.x;
        y += v.y;
        return *this;
    }

    R2Vector operator-(const R2Point& p) const {
        return R2Vector(x-p.x, y-p.y);
    }

    R2Point operator-(const R2Vector& v) const {
        return R2Point(x-v.x, y-v.y);
    }

    R2Point& operator-=(const R2Vector& v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }

    R2Point& operator-=(const R2Point& p) {
        x -= p.x;
        y -= p.y;
        return *this;
    }

    R2Point operator*(double c) const {
        return R2Point(x*c, y*c);
    }

    R2Point& operator*=(double c) {
        x *= c;
        y *= c;
        return *this;
    }

    // Comparings
    bool operator==(const R2Point& p) const {
        //... return (x == p.x && y == p.y);
        return (
            fabs(x - p.x) <= R2_EPSILON &&
            fabs(y - p.y) <= R2_EPSILON
        );
    }
    bool closeTo(const R2Point& v, double eps = R2_EPSILON) const {
        return (
            fabs(x - v.x) <= eps &&
            fabs(y - v.y) <= eps
        );
    }
    bool operator!=(const R2Point& p) const { return !operator==(p); }
    bool operator>=(const R2Point& p) const {
        //... return (x > p.x || (x == p.x && y >= p.y));
        return (x > p.x || (x >= p.x && y >= p.y));
    }
    bool operator>(const R2Point& p) const {
        //... return (x > p.x || (x == p.x && y > p.y));
        return (x > p.x || (x >= p.x && y > p.y));
    }
    bool operator<(const R2Point& p) const { return !operator>=(p); }
    bool operator<=(const R2Point& p) const { return !operator>(p); }

    // Area of oriented triangle
    static double signed_area(
        const R2Point& a, const R2Point& b, const R2Point& c
    ) {
        return 0.5 * R2Vector::signed_area(b-a, c-a);
    }

    static double area(
        const R2Point& a, const R2Point& b, const R2Point& c
    ) {
        return fabs(signed_area(a, b, c));
    }
    static double Area(
        const R2Point& a, const R2Point& b, const R2Point& c
    ) {
        return area(a, b, c);
    }

    bool between(const R2Point& a, const R2Point& b) const {
        R2Vector v(b - a);
        R2Vector m(*this - a);
        return (
            fabs(v.normal() * m) <= R2_EPSILON &&   // point on line(a, b)
            m * v >= 0. &&  (*this - b) * v <= 0.   // between (a, b)
        );
    }

    static bool on_line(
        const R2Point& a, const R2Point& b, const R2Point& c
    ) {
        return (area(a, b, c) <= R2_EPSILON);
    }

    // Angle from this point between points a and b (counterclockwise)
    double angle(const R2Point& a, const R2Point& b) const {
        return (a - *this).angle(b - *this);
    }
    double Angle(const R2Point& a, const R2Point& b) const {
        return angle(a, b);
    }

    // Angle with vertex A from AB to AC counterclockwise
    static double angle(
        const R2Point& A, const R2Point& B, const R2Point& C
    ) {
        return A.angle(B, C);
    }

    double distance(const R2Point& p) const {
        return (p - *this).length();
    }
    double Distance(const R2Point& p) const { return distance(p); }

    static double distance(const R2Point& a, const R2Point& b) {
        return a.distance(b);
    }
    static double Distance(const R2Point& a, const R2Point& b) {
        return distance(a, b);
    }

    // Square of distance
    double distanceSquare(const R2Point& p) const {
        double dx = p.x - x;
        double dy = p.y - y;
        return (dx*dx + dy*dy);
    }
    static double distanceSquare(
        const R2Point& p0, const R2Point& p1
    ) {
        return p0.distanceSquare(p1);
    }

    // Max. of abs. values of coord. differencies
    double distanceAbs(const R2Point& p) const {
        double dx = fabs(p.x - x);
        double dy = fabs(p.y - y);
        if (dx >= dy)
            return dx;
        else
            return dy;
    }
    static double distanceAbs(
        const R2Point& p0, const R2Point& p1
    ) {
        return p0.distanceAbs(p1);
    }

    double distanceToLine(
        const R2Point& p, const R2Vector& v
    ) const;

    double distanceToLine(
        const R2Point& p0, const R2Point& p1
    ) const {
        return distanceToLine(p0, p1 - p0);
    }

    double distanceToLineSegment(
        const R2Point& p0, const R2Point& p1
    ) const;
};

inline R2Point operator*(double c, const R2Point& p) {
    return R2Point(c*p.x, c*p.y);
}

// The R2Point with strict comparings
class R2Point0: public R2Point {
public:
    R2Point0(const R2Point& p): R2Point(p) {}
    R2Point0& operator=(const R2Point& p) {
        x = p.x; y = p.y;
    }
    //... operator R2Point() { return R2Point(x, y); }

    // Comparings
    bool operator==(const R2Point& p) const {
        //... return (
        //...     fabs(x - p.x) <= R2_EPSILON &&
        //...     fabs(y - p.y) <= R2_EPSILON
        //... );
        return (x == p.x && y == p.y);
    }
    bool operator!=(const R2Point& p) const { return !operator==(p); }
};

// Un-oriented line segment in R2-plane
class R2LineSegment {
public:
    static double precision;
    R2Point ends[2];            // Ordered: ends[0] <= ends[1]

    static void setPrecision(double prec);

    R2LineSegment(const R2Point& p0, const R2Point& p1) {
        if (p0 <= p1) {
            ends[0] = p0; ends[1] = p1;
        } else {
            ends[0] = p1; ends[1] = p0;
        }
    }

    R2LineSegment(const R2Point v[2]) {
        if (v[0] <= v[1]) {
            ends[0] = v[0]; ends[1] = v[1];
        } else {
            ends[0] = v[1]; ends[1] = v[0];
        }
    }

    bool operator==(const R2LineSegment& s) const {
        return (
            //... ends[0] == s.ends[0] &&
            //... ends[1] == s.ends[1]
            ends[0].distanceAbs(s.ends[0]) <= precision &&
            ends[1].distanceAbs(s.ends[1]) <= precision
        );
    }
    bool operator!=(const R2LineSegment& s) const {
        return (!operator==(s));
    }
    bool operator<=(const R2LineSegment& s) const {
        return (
            operator==(s)
            ||
            (
                ends[0] < s.ends[0] ||
                (ends[0] <= s.ends[0] && ends[1] <= s.ends[1])
            )
        );
    }
    bool operator<(const R2LineSegment& s) const {
        return (
            !operator==(s) &&
            (
                ends[0] < s.ends[0] ||
                (ends[0] <= s.ends[0] && ends[1] < s.ends[1])
            )
        );
    }
    bool operator>(const R2LineSegment& s) const {
        return (!operator<=(s));
    }
    bool operator>=(const R2LineSegment& s) const {
        return (!operator<(s));
    }
};

// Oriented line segment in R2-plane
class R2OrientedLineSegment {
public:
    static double precision;
    R2Point ends[2];            // Unordered

    static void setPrecision(double prec);

    R2OrientedLineSegment(const R2Point& p0, const R2Point& p1) {
        ends[0] = p0; ends[1] = p1;
    }

    R2OrientedLineSegment(const R2Point v[2]) {
        ends[0] = v[0]; ends[1] = v[1];
    }

    bool operator==(const R2OrientedLineSegment& s) const {
        return (
            //... ends[0] == s.ends[0] &&
            //... ends[1] == s.ends[1]
            ends[0].distanceAbs(s.ends[0]) <= precision &&
            ends[1].distanceAbs(s.ends[1]) <= precision
        );
    }
    bool operator!=(const R2OrientedLineSegment& s) const {
        return (!operator==(s));
    }
    bool operator<=(const R2OrientedLineSegment& s) const {
        return (
            operator==(s)
            ||
            (
                ends[0] < s.ends[0] ||
                (ends[0] <= s.ends[0] && ends[1] <= s.ends[1])
            )
        );
    }
    bool operator<(const R2OrientedLineSegment& s) const {
        return (
            !operator==(s) &&
            (
                ends[0] < s.ends[0] ||
                (ends[0] <= s.ends[0] && ends[1] < s.ends[1])
            )
        );
    }
    bool operator>(const R2OrientedLineSegment& s) const {
        return (!operator<=(s));
    }
    bool operator>=(const R2OrientedLineSegment& s) const {
        return (!operator<(s));
    }
};

class R3Vector
{
public:
    double x, y, z;

public:
    R3Vector() : x( 0. ), y( 0. ), z( 0. ) {}
    R3Vector( double X, double Y, double Z ) : x( X ), y( Y ), z( Z ) {}
    R3Vector( const R3Vector& v ) : x( v.x ), y( v.y ), z( v.z ) {}
    R3Vector( const ByteR3Vector& v );
    ~R3Vector() {}

    R3Vector& operator=( const R3Vector& v )
    {
        x = v.x; y = v.y; z = v.z; return *this;
    }

    bool operator==(const R3Vector& v) const
    {
        /*...
        return (
            x == v.x && y == v.y && z == v.z
        );
        ...*/
        return (
            fabs(x - v.x) <= R3_EPSILON &&
            fabs(y - v.y) <= R3_EPSILON &&
            fabs(z - v.z) <= R3_EPSILON
        );
    }
    bool closeTo(const R3Vector& v, double eps = R3_EPSILON) const {
        return (
            fabs(x - v.x) <= eps &&
            fabs(y - v.y) <= eps &&
            fabs(z - v.z) <= eps
        );
    }
    bool operator!=(const R3Vector& v) const { return !operator==(v); }
    bool operator>=(const R3Vector& v) const {
        return (
            x > v.x ||
            (x >= v.x && y > v.y) ||
            (x >= v.x && y >= v.y && z >= v.z)
        );
    }
    bool operator>(const R3Vector& v) const {
        return (
            x > v.x ||
            (x >= v.x && y > v.y) ||
            (x >= v.x && y >= v.y && z > v.z)
        );
    }
    bool operator<(const R3Vector& v) const { return !operator>=(v); }
    bool operator<=(const R3Vector& v) const { return !operator>(v); }

    R3Vector& operator=( const ByteR3Vector& v );

    R3Vector operator+( const R3Vector& v ) const
    {
        return R3Vector( x + v.x, y + v.y, z + v.z );
    }

    R3Vector operator-( const R3Vector& v ) const
    {
        return R3Vector( x - v.x, y - v.y, z - v.z );
    }

    const R3Vector& operator+=( const R3Vector& v )
    {
        x += v.x; y += v.y; z += v.z; return *this;
    }

    const R3Vector& operator-=( const R3Vector& v )
    {
        x -= v.x; y -= v.y; z -= v.z; return *this;
    }

    R3Vector operator*( double r ) const
    {
        return R3Vector( r * x,  r * y, r * z );
    }

    const R3Vector& operator*=( double r )
    {
        x *= r;  y *= r; z *= r; return *this;
    }

    // Scalar product
    double operator*( const R3Vector& v ) const
    {
        return  x * v.x + y * v.y + z * v.z;
    }

    double ScalarProduct(const R3Vector& v) const
    {
        return operator*(v);
    }

    // Vector product
    R3Vector operator&(const R3Vector& v) const
    {
        return VectorProduct(v);
    }
    R3Vector VectorProduct(const R3Vector& v) const
    {
        return R3Vector(
            y * v.z - z * v.y,
            -(x * v.z) + z * v.x,
            x * v.y - y * v.x
        );
    }
    R3Vector vectorProduct(const R3Vector& v) const
    {
        return VectorProduct(v);
    }

    double Length() const { return sqrt( x*x + y*y + z*z ); }
    double length() const { return Length(); }
    double lengthSquare() const { return (x*x + y*y + z*z); }
    double LengthSquare() const { return lengthSquare(); }

    R3Vector& operator*=(const R3Matrix&);

    // Make length = 1.
    R3Vector& Normalize()
    {
        double l = Length();
        if (l <= R3_EPSILON)
            return *this;
        l = 1./l;
        x *= l;
        y *= l;
        z *= l;
        return *this;
    }
    R3Vector& normalize() { return Normalize(); }

    // Make length = RNormalLength (120.)
    R3Vector& Normalize64()
    {
        double l = Length();
        if (l <= R3_EPSILON)
            return *this;
        l = RNormalLength / l;
        x *= l;
        y *= l;
        z *= l;
        return *this;
    }

    bool Parallel(const R3Vector& v) const
    {
        R3Vector p = VectorProduct(v);
        return (
            fabs(p.x) <= R3_EPSILON &&
            fabs(p.y) <= R3_EPSILON &&
            fabs(p.z) <= R3_EPSILON
        );
    }

    R3Vector Normal() const     // Return arbitrary normal to vector
    {
        if (
            fabs(x) <= R3_EPSILON &&
            fabs(y) <= R3_EPSILON &&
            fabs(z) <= R3_EPSILON
        )
            return R3Vector(1., 0., 0.);
        if (fabs(x) >= fabs(y))
        {
            if (fabs(x) >= fabs(z))
                return R3Vector(-y, x, 0.);     // x maximal
            else
                return R3Vector(z, 0., -x);     // z maximal
        }
        else
        {
            if (fabs(y) >= fabs(z))
                return R3Vector(y, -x, 0.);     // y maximal
            else
                return R3Vector(0., z, -y);     // z maximal
        }
    }

    double Angle(const R3Vector& v) const;
    double angle(const R3Vector& v) const { return Angle(v); }

    double area(const R3Vector& v) const
    {
        return vectorProduct(v).length();
    }
    double Area(const R3Vector& v) const { return area(v); }

    static double signedVolume(    // Determinant
        const R3Vector& a,
        const R3Vector& b,
        const R3Vector& c
    );
};

inline R3Vector operator*(double c, const R3Vector& v) {
    return R3Vector(c*v.x, c*v.y, c*v.z);
}

class R2Matrix
{
public:
    double m[2][2];

public:
    R2Matrix()
    {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                m[i][j] = 0.;
    }

    R2Matrix(const R2Matrix& M)
    {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                m[i][j] = M.m[i][j];
    }

    R2Matrix(double a11, double a12, double a21, double a22)
    {
        m[0][0] = a11; m[0][1] = a12;
        m[1][0] = a21; m[1][1] = a22;
    }

    R2Matrix(double alpha)      // Rotation matrix
    {
        double cosAlpha = cos(alpha);
        double sinAlpha = sin(alpha);
        m[0][0] = cosAlpha; m[0][1] = (-sinAlpha);
        m[1][0] = sinAlpha; m[1][1] = cosAlpha;
    }

    ~R2Matrix() {}

    const R2Matrix& operator=(const R2Matrix& M)
    {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                m[i][j] = M.m[i][j];
        return *this;
    }

    double& operator()(int i, int j)
    {
        return m[i][j];
    }

    R2Vector operator*(const R2Vector& v) const
    {
        return R2Vector(
            m[0][0] * v.x + m[0][1] * v.y,
            m[1][0] * v.x + m[1][1] * v.y
        );
    }

    R2Matrix operator*( const R2Matrix& M ) const
    {
        return R2Matrix(
            m[0][0] * M.m[0][0] + m[0][1] * M.m[1][0],
            m[0][0] * M.m[0][1] + m[0][1] * M.m[1][1],
            m[1][0] * M.m[0][0] + m[1][1] * M.m[1][0],
            m[1][0] * M.m[0][1] + m[1][1] * M.m[1][1]
        );
    }

    double* operator[](int i) { return m[i]; }
    const double* operator[](int i) const { return m[i]; }

    double& ElementAt(int i, int j) { return m[i][j]; }
    double& at(int i, int j) { return ElementAt(i, j); }
    const double& ElementAt(int i, int j) const { return m[i][j]; }
    const double& at(int i, int j) const { return ElementAt(i, j); }
    double GetAt(int i, int j) const { return m[i][j]; }
    void SetAt(int i, int j, double v) { m[i][j] = v; }

    double Determinant() const
    {
        return (m[0][0] * m[1][1] - m[1][0] * m[0][1]);
    }

    bool Invert(R2Matrix& inv) const
    {
        double det = Determinant();
        //... if (det == 0.)
        if (fabs(det) <= R2_EPSILON)
            return false;

        inv.m[0][0] = m[1][1] / det;
        inv.m[1][0] = (-m[1][0]) / det;
        inv.m[0][1] = (-m[0][1]) / det;
        inv.m[1][1] = m[0][0] / det;

        return true;
    }

    static R2Matrix Zero()
    {
        return R2Matrix(
            0., 0.,
            0., 0.
        );
    }

    static R2Matrix Unit()
    {
        return R2Matrix(
            1., 0.,
            0., 1.
        );
    }

    void Erase()
    {
        m[0][0] = 0.;
        m[0][1] = 0.;
        m[1][0] = 0.;
        m[1][1] = 0.;
    }

    void SetZero() { Erase(); }
    void erase() { Erase(); }

    void SetUnit()
    {
        m[0][0] = 1.;
        m[0][1] = 0.;
        m[1][0] = 0.;
        m[1][1] = 1.;
    }

    void GaussMethod(
        int& rank, double& det, double eps = 1e-7
    );
};

inline R2Vector& R2Vector::operator*=(const R2Matrix& m)
{
    double x1, y1;
    x1 = x * m.m[0][0] + y * m.m[1][0];
    y1 = x * m.m[0][1] + y * m.m[1][1];
    x = x1; y = y1;
    return *this;
}

class R3Matrix
{
public:
    double m[3][3];

public:
    R3Matrix()
    {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m[i][j] = 0.;
    }

    R3Matrix(const R3Matrix& M)
    {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m[i][j] = M.m[i][j];
    }

    R3Matrix(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33
    )
    {
        m[0][0] = a11; m[0][1] = a12; m[0][2] = a13;
        m[1][0] = a21; m[1][1] = a22; m[1][2] = a23;
        m[2][0] = a31; m[2][1] = a32; m[2][2] = a33;
    }

    R3Matrix(const R3Vector v[3])
    {
        m[0][0] = v[0].x; m[0][1] = v[1].x; m[0][2] = v[2].x;
        m[1][0] = v[0].y; m[1][1] = v[1].y; m[1][2] = v[2].y;
        m[2][0] = v[0].z; m[2][1] = v[1].z; m[2][2] = v[2].z;
    }

    R3Matrix(
        const R3Vector& v0,
        const R3Vector& v1,
        const R3Vector& v2
    )
    {
        m[0][0] = v0.x; m[0][1] = v1.x; m[0][2] = v2.x;
        m[1][0] = v0.y; m[1][1] = v1.y; m[1][2] = v2.y;
        m[2][0] = v0.z; m[2][1] = v1.z; m[2][2] = v2.z;
    }

    R3Matrix(
        const R3Vector& rotationAxis,
        double alpha
    );

    static R3Matrix RotationMatrix(
        const R3Vector& rotationAxis,
        double alpha
    );

    ~R3Matrix() {}

    double& ElementAt(int i, int j) { return m[i][j]; }
    double& at(int i, int j) { return ElementAt(i, j); }
    const double& ElementAt(int i, int j) const { return m[i][j]; }
    const double& at(int i, int j) const { return ElementAt(i, j); }
    double GetAt(int i, int j) const { return m[i][j]; }
    void SetAt(int i, int j, double v) { m[i][j] = v; }
    double* operator[](int i) { return m[i]; }
    const double* operator[](int i) const { return m[i]; }

    const R3Matrix& operator=(const R3Matrix& M)
    {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m[i][j] = M.m[i][j];
        return *this;
    }

    double& operator() (int i, int j)
    {
        return m[i][j];
    }

    R3Vector operator*( const R3Vector& v ) const
    {
        return R3Vector(
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z
        );
    }

    ByteR3Vector operator*( const ByteR3Vector& v ) const;

    R3Matrix operator*( const R3Matrix& M ) const
    {
        R3Matrix r;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                r.m[i][j] = 0;
                for (int k = 0; k < 3; k++)
                {
                    r.m[i][j] += m[i][k] * M.m[k][j];
                }
            }
        }
        return r;
    }

    R3Matrix& operator*=( const R3Matrix& M )
    {
        R3Matrix r(*this);
        *this = r * M;
        return *this;
    }

    double Determinant() const
    {
        return (
            m[0][0] * m[1][1] * m[2][2] +
            m[1][0] * m[2][1] * m[0][2] +
            m[2][0] * m[0][1] * m[1][2] -

            m[0][2] * m[1][1] * m[2][0] -
            m[1][0] * m[2][2] * m[0][1] -
            m[2][1] * m[1][2] * m[0][0]
        );
    }

    double determinant() const { return Determinant(); }

    bool Invert(R3Matrix& inv) const
    {
        double det = Determinant();
        //... if (det == 0.)
        if (fabs(det) <= R2_EPSILON)
            return false;

        inv.m[0][0] =     //  A11/det
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
        inv.m[1][0] =     // -A01/det
            (-m[1][0] * m[2][2] + m[2][0] * m[1][2]) / det;
        inv.m[2][0] =     //  A02/det
            (m[1][0] * m[2][1] - m[2][0] * m[1][1]) / det;

        inv.m[0][1] =     // -A10/det
            (-m[0][1] * m[2][2] + m[2][1] * m[0][2]) / det;
        inv.m[1][1] =     // A11/det
            (m[0][0] * m[2][2] - m[2][0] * m[0][2]) / det;
        inv.m[2][1] =     // -A12/det
            (-m[0][0] * m[2][1] + m[2][0] * m[0][1]) / det;

        inv.m[0][2] =     //  A20/det
            (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / det;
        inv.m[1][2] =     // -A21/det
            (-m[0][0] * m[1][2] + m[1][0] * m[0][2]) / det;
        inv.m[2][2] =     //  A22/det
            (m[0][0] * m[1][1] - m[1][0] * m[0][1]) / det;

        return true;
    }

    void SetZero()
    {
        m[0][0] = 0.; m[0][1] = 0.; m[0][2] = 0.;
        m[1][0] = 0.; m[1][1] = 0.; m[1][2] = 0.;
        m[2][0] = 0.; m[2][1] = 0.; m[2][2] = 0.;
    }
    void Erase() { SetZero(); }
    void erase() { SetZero(); }

    void SetUnit()
    {
        m[0][0] = 1.; m[0][1] = 0.; m[0][2] = 0.;
        m[1][0] = 0.; m[1][1] = 1.; m[1][2] = 0.;
        m[2][0] = 0.; m[2][1] = 0.; m[2][2] = 1.;
    }

    static R3Matrix Zero()
    {
        return R3Matrix(
            0., 0., 0.,
            0., 0., 0.,
            0., 0., 0.
        );
    }

    static R3Matrix Unit()
    {
        return R3Matrix(
            1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.
        );
    }

    bool IsRestrictedRotation(double maxAngle) const;

    void GaussMethod(
        int& rank, double& det, double eps = 1e-7
    );
};

inline double R3Vector::signedVolume(    // Determinant
    const R3Vector& a,
    const R3Vector& b,
    const R3Vector& c
)
{
    return R3Matrix(a, b, c).determinant();
}

inline R3Vector& R3Vector::operator*=(const R3Matrix& m)
{
    double x1, y1, z1;
    x1 = x * m.m[0][0] + y * m.m[1][0] + z * m.m[2][0];
    y1 = x * m.m[0][1] + y * m.m[1][1] + z * m.m[2][1];
    z1 = x * m.m[0][2] + y * m.m[1][2] + z * m.m[2][2];
    x = x1; y = y1; z = z1;
    return *this;
}

class ByteR3Vector
{
public:
    signed char x, y, z;

public:
    ByteR3Vector() : x(0), y(0), z(0) {}

    ByteR3Vector( int X, int Y, int Z )
        :
            x(char(X)),
            y(char(Y)),
            z(char(Z))
    {}

    ByteR3Vector( double X, double Y, double Z )
        :
            x(char(int(X))),
            y(char(int(Y))),
            z(char(int(Z)))
    {}

    ByteR3Vector( const ByteR3Vector& v ) : x( v.x ), y( v.y ), z( v.z ) {}

    ByteR3Vector( const R3Vector& v )
        :
            x(char(int(v.x))),
            y(char(int(v.y))),
            z(char(int(v.z)))
    {}

    ~ByteR3Vector() {}

    ByteR3Vector& Normal64(R3Vector v)
    {
        // Make the lengh to be RNormalLength=120 approximately
        double l = (double)v.x * (double)v.x +
            (double)v.y * (double)v.y + (double)v.z * (double)v.z;
        if (l > 0.)
        {
            l = RNormalLength / sqrt(l);
            x = char(int( v.x * l ));
            y = char(int( v.y * l ));
            z = char(int( v.z * l ));
        }
        else
        {
            x = char(int(v.x));
            y = char(int(v.y));
            z = char(int(v.z));
        }
        return *this;
    }

    ByteR3Vector& Normal64(double X, double Y, double Z)
    {
        // Make the lengh to be RNormalLength=120 approximately
        double l = X*X + Y*Y + Z*Z;
        if (l > 0.)
        {
            l = RNormalLength / sqrt(l);
            x = char(int( X * l ));
            y = char(int( Y * l ));
            z = char(int( Z * l ));
        }
        else
        {
            x = char(int(X));
            y = char(int(Y));
            z = char(int(Z));
        }
        return *this;
    }

    ByteR3Vector& Normal64(int X, int Y, int Z)
    {
        // Make the lengh to be RNormalLength=120 approximately
        double l = (double)X*(double)X +
            (double)Y*(double)Y + (double)Z*(double)Z;
        if (l > 0.)
        {
            l = RNormalLength / sqrt(l);
            x = char(int( (double)X * l ));
            y = char(int( (double)Y * l ));
            z = char(int( (double)Z * l ));
        }
        else
        {
            x = char(X);
            y = char(Y);
            z = char(Z);
        }
        return *this;
    }

    ByteR3Vector& Normalize64()
    {
        double l = (double)x * (double)x +
            (double)y * (double)y + (double)z * (double)z;
        if (l <= R3_EPSILON)
            return *this;
        l = RNormalLength / sqrt(l);
        x = (char) ((double) x * l);
        y = (char) ((double) y * l);
        z = (char) ((double) z * l);
        return *this;
    }

    const ByteR3Vector& operator=( const ByteR3Vector& v )
    {
        x = v.x; y = v.y; z = v.z; return *this;
    }

    const ByteR3Vector& operator=( const R3Vector& v )
    {
        x = char(int(v.x)); y = char(int(v.y)); z = char(int(v.z));
        return *this;
    }

    ByteR3Vector operator+( const ByteR3Vector& v ) const
    {
        return ByteR3Vector( x + v.x, y + v.y, z + v.z );
    }

    ByteR3Vector operator-( const ByteR3Vector& v ) const
    {
        return ByteR3Vector( x - v.x, y - v.y, z - v.z );
    }

    ByteR3Vector& operator+=( const ByteR3Vector& v )
    {
        x += v.x; y += v.y; z += v.z; return *this;
    }

    ByteR3Vector& operator-=( const ByteR3Vector& v )
    {
        x -= v.x; y -= v.y; z -= v.z; return *this;
    }

    ByteR3Vector operator*( double r ) const
    {
        return ByteR3Vector( r * double(x),  r * double(y), r * double(z) );
    }

    ByteR3Vector& operator*=( double r )
    {
        x = char(int(r * double(x)));
        y = char(int(r * double(y)));
        z = char(int(r * double(y)));
        return *this;
    }

    // Scalar product
    double operator*( const ByteR3Vector& v ) const
    {
        return (
            double(x) * double(v.x) +
            double(y) * double(v.y) +
            double(z) * double(v.z)
        );
    }

    // Vector product
    ByteR3Vector operator&(const ByteR3Vector& v) const
    {
        return ByteR3Vector(
            int(y) * int(v.z) - int(z) * int(v.y),
            -(int(x) * int(v.z)) + int(z) * int(v.x),
            int(x) * int(v.y) - int(y) * int(v.x)
        );
    }

    ByteR3Vector VectorProduct(const ByteR3Vector& v) const
    {
        return operator&(v);
    }

    double Length() const
    {
        return sqrt(double(
            int(x)*int(x) + int(y)*int(y) + int(z)*int(z)
        ));
    }

    ByteR3Vector& operator*=(const R3Matrix& m)
    {
        double x1, y1, z1;
        x1 = double(x) * m.m[0][0] + double(y) * m.m[1][0] +
            double(z) * m.m[2][0];
        y1 = double(x) * m.m[0][1] + double(y) * m.m[1][1] +
            double(z) * m.m[2][1];
        z1 = double(x) * m.m[0][2] + double(y) * m.m[1][2] +
            double(z) * m.m[2][2];
        x = char(int(x1)); y = char(int(y1)); z = char(int(z1));
        return *this;
    }
};

inline ByteR3Vector R3Matrix::operator*(const ByteR3Vector& v) const
{
    return ByteR3Vector(
        m[0][0] * double(v.x) + m[0][1] * double(v.y) + m[0][2] * double(v.z),
        m[1][0] * double(v.x) + m[1][1] * double(v.y) + m[1][2] * double(v.z),
        m[2][0] * double(v.x) + m[2][1] * double(v.y) + m[2][2] * double(v.z)
    );
}

inline R3Vector::R3Vector( const ByteR3Vector& v ):
    x(double(int(v.x))),
    y(double(int(v.y))),
    z(double(int(v.z)))
{
}

inline R3Vector& R3Vector::operator=( const ByteR3Vector& v )
{
    x = double(int(v.x));
    y = double(int(v.y));
    z = double(int(v.z));
    return *this;
}

class R3Point
{
public:
    double x, y, z;

public:
    R3Point(): x( 0. ), y( 0. ), z( 0. ) {}
    R3Point(double X, double Y, double Z): x( X ), y( Y ), z( Z ) {}
    R3Point(const R3Point& v): x( v.x ), y( v.y ), z( v.z ) {}
    R3Point(const R3Vector& v): x( v.x ), y( v.y ), z( v.z ) {}
    ~R3Point() {}

    R3Point& operator=(const R3Point& v)
    {
        x = v.x; y = v.y; z = v.z; return *this;
    }

    bool operator==(const R3Point& v) const
    {
        /*...
        return (
            x == v.x && y == v.y && z == v.z
        );
        ...*/
        return (
            fabs(x - v.x) <= R3_EPSILON &&
            fabs(y - v.y) <= R3_EPSILON &&
            fabs(z - v.z) <= R3_EPSILON
        );
    }
    bool closeTo(const R3Point& v, double eps = R3_EPSILON) const {
        return (
            fabs(x - v.x) <= eps &&
            fabs(y - v.y) <= eps &&
            fabs(z - v.z) <= eps
        );
    }
    bool operator!=(const R3Point& v) const { return !operator==(v); }
    bool operator>=(const R3Point& v) const {
        return (
            x > v.x ||
            (x >= v.x && y > v.y) ||
            (x >= v.x && y >= v.y && z >= v.z)
        );
    }
    bool operator>(const R3Point& v) const {
        return (
            x > v.x ||
            (x >= v.x && y > v.y) ||
            (x >= v.x && y >= v.y && z > v.z)
        );
    }
    bool operator<(const R3Point& v) const { return !operator>=(v); }
    bool operator<=(const R3Point& v) const { return !operator>(v); }

    R3Point operator+(const R3Vector& v) const
    {
        return R3Point(x + v.x, y + v.y, z + v.z);
    }

    R3Point operator-(const R3Vector& v) const
    {
        return R3Point(x - v.x, y - v.y, z - v.z);
    }

    R3Vector operator-(const R3Point& p) const
    {
        return R3Vector(x - p.x, y - p.y, z - p.z);
    }

    const R3Point& operator+=(const R3Vector& v)
    {
        x += v.x; y += v.y; z += v.z; return *this;
    }

    const R3Point& operator-=(const R3Vector& v)
    {
        x -= v.x; y -= v.y; z -= v.z; return *this;
    }

    double signedSolidAngle(
        const R3Point& a, const R3Point& b, const R3Point& c
    ) const
    {
        R3Vector va = a - *this;
        R3Vector vb = b - *this;
        R3Vector vc = c - *this;
        double det = R3Matrix(va, vb, vc).determinant();
        double s;
        if (det < 0.) {
            s = (-1.);
            det = (-det);
        } else {
            s = 1.;
        }
        if (det <= R3_EPSILON)
            return 0.;

        double al = va.length();
        double bl = vb.length();
        double cl = vc.length();

        double div = al*bl*cl +
            (va*vb)*cl +    // Scalar product of vectors va and vb
            (va*vc)*bl +
            (vb*vc)*al;
        double at = atan2(fabs(det), div);
        if (at < 0.)
            at += PI;
        return (s * 2. * at);
    }

    double signedSolidAngle(const R3Point t[3]) const
    {
        return signedSolidAngle(t[0], t[1], t[2]);
    }

    double solidAngle(  // always non-negative
        const R3Point& a, const R3Point& b, const R3Point& c
    ) const
    {
        return fabs(signedSolidAngle(a, b, c));
    }

    double solidAngle(  // always non-negative
        const R3Point t[3]
    ) const
    {
        return fabs(signedSolidAngle(t[0], t[1], t[2]));
    }

    // Area of triangle
    double area(const R3Point& b, const R3Point& c) const
    {
        return ((b - *this).area(c - *this)) * 0.5;
    }

    static double area(
        const R3Point& a, const R3Point& b, const R3Point& c
    )
    {
        return a.area(b, c);
    }
};

// Fight against MFC
#include <vector>
typedef std::vector<R3Point> R3PointArray;

// Global functions

bool IntersectR2Lines(
    const R2Point& p0, const R2Vector& v0,
    const R2Point& p1, const R2Vector& v1,
    R2Point& p,
    double* t0 = NULL, double* t1 = NULL
);

bool IntersectR2LineAndLineSegment(
    const R2Point& p0, const R2Vector& v0,
    const R2Point& p1, const R2Point& q1,
    R2Point& p
);

bool IntersectR2LineSegments(
    const R2Point& p0, R2Point& q0,
    const R2Point& p1, R2Point& q1,
    R2Point& p
);

//
// Intersect a plane, defined by a point p0 and normal vector v0,
// and a circumference,
// defined by a center p1, radius r, and normal vector v1.
// Results are two points q1 and q2.
// Return value: true, if vectors v0 and v1 are not collinear,
// and an intersection is not empty
//
BOOL IntersectPlaneAndCircle(
    const R3Point& p0, const R3Vector& v0,              // Plane
    const R3Point& p1, const R3Vector& v1, double r,    // Circumference
    R3Point& q1, R3Point& q2
);

// MFC Container classes

//MFC typedef CArray <R3Point, R3Point&> R3PointArray;

// Helper functions
/*...
template <> void AFXAPI ConstructElements <R3Point>
    (R3Point* pElements, int nCount);
template <> void AFXAPI DestructElements <R3Point>
    (R3Point* pElements, int nCount);
...*/

BOOL IntersectPlaneAndCylinder(             // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cylinder
    double diameter, double height,

    R3PointArray& section
);

BOOL IntersectAxialPlaneAndCylinder(        // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cylinder
    double diameter, double height,

    R3PointArray& section
);

void DrawCylinder3DShape(
    McDrawable* pDraw, McGC& gc,
    double shiftX, double shiftY,
    double zoomX, double zoomY,

    const R3Point& planeOrigin,
    const R3Vector& plane_ex,
    const R3Vector& plane_ey,
    const R3Vector& plane_ez,

    const R3Point& implantPosition,
    const R3Vector& implantAxis,
    double diameter,
    double height,

    BOOL lowerHalf = FALSE,
    BOOL upperHalf = FALSE
);

void DrawSharpCylinder3DShape(
    McDrawable* pDraw, McGC& gc,
    double shiftX, double shiftY,
    double zoomX, double zoomY,

    const R3Point& planeOrigin,
    const R3Vector& plane_ex,
    const R3Vector& plane_ey,
    const R3Vector& plane_ez,

    const R3Point& implantPosition,
    const R3Vector& implantAxis,
    double diameter,
    double height,              // Total height
    double coneHeight,           // included in height

    BOOL upperHalf = FALSE
);

void DrawCone3DShape(
    McDrawable* pDraw, McGC& gc,
    double shiftX, double shiftY,
    double zoomX, double zoomY,

    const R3Point& planeOrigin,
    const R3Vector& plane_ex,
    const R3Vector& plane_ey,
    const R3Vector& plane_ez,

    const R3Point& implantPosition,
    const R3Vector& implantAxis,
    double upper_diameter,
    double lower_diameter,
    double height,

    BOOL lowerHalf = FALSE,
    BOOL upperHalf = FALSE
);

BOOL IntersectPlaneAndCone(         // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cone center, axis
    double upper_diameter, double lower_diameter,
    double height,

    R3PointArray& section
);

BOOL IntersectAxialPlaneAndCone(        // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cylinder
    double upper_diameter, double lower_diameter,
    double height,

    R3PointArray& section
);

// Globals
const int ELLIPSE_STEPS = 72;   // 5 degrees

// Angle step for ellipse drawing
extern double ELLIPSE_SIN[ELLIPSE_STEPS];
extern double ELLIPSE_COS[ELLIPSE_STEPS];
extern BOOL ellipseSinCosCalculated;
extern void CalculateEllipseSinCos();

// Intersection of plane and ...
bool IntersectPlanes(
    const R3Point& p0, const R3Vector& n0,
    const R3Point& p1, const R3Vector& n1,
    R3Point& p, R3Vector& v
);

inline bool IntersectPlaneAndLine(
    const R3Point& p0, const R3Vector& n0,
    const R3Point& p1, const R3Vector& v,
    double& t   // Point internal coordinate in a straight line
)
{
    if (fabs(n0 * v) <= R3_EPSILON)
        return false;
    // p = p1 + v*t,    t is a number
    // (p1 + v*t - p0, n0) = 0
    // t = (p1 - p0, n0) / (v, n0)
    t = ((p1 - p0)*n0) / (v*n0);
    return true;
}

inline bool IntersectPlaneAndLine(
    const R3Point& p0, const R3Vector& n0,
    const R3Point& p1, const R3Vector& v,
    R3Point& t
)
{
    double x;
    if (!IntersectPlaneAndLine(p0, n0, p1, v, x))
        return false;
    t = p1 + (v*x);
}

// Returns number of points in intersection
// (when line segment is in the plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p, const R3Vector& n,
    const R3Point& s0, const R3Point& s1,
    R3Point& t0, R3Point& t1
);

// Returns number of points in intersection
// (when line segment is in the plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Point& s0, const R3Point& s1,
    R2Point& t0, R2Point& t1    // Points in internal coordinate system
);

// Returns number of points in intersection
// (when line segment is in the plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Vector& ez,         // Normal to the plane == ex.VectorProduct(ey)
    const R3Point& s0, const R3Point& s1,
    R2Point& t0, R2Point& t1    // Points in internal coordinate system
);

inline bool IntersectPlaneAndLineSegment(
    const R3Point& p, const R3Vector& n,
    const R3Point& s0, const R3Vector& v,
    double t0, double t1,   // Ends of line segment in internal coord.
    double& t
)
{
    if (!IntersectPlaneAndLine(p, n, s0, v, t))
        return false;
    else
        return(t0 <= t && t <= t1);
}

// Returns the number of points in intersection with
// the border of triangle
int IntersectPlaneAndTriangle(
    const R3Point& p, const R3Vector& n,
    const R3Point triangeVertices[3],
    R3Point intersection[2]
);

// Returns the number of points in intersection with
// the border of triangle
int IntersectPlaneAndTriangle(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
);

int IntersectPlaneAndTriangle(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Vector& ez,         // Normal to the plane == ex.VectorProduct(ey)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
);

int IntersectPlaneAndOrientedTriangle(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Vector& ez,         // Normal to the plane == ex.VectorProduct(ey)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
);

bool PointInPlane(
    const R3Point& p0,
    const R3Point& p1, const R3Vector& n
);

bool IntersectR3Rectangles(
    const R3Point& leftTop0,
    const R3Vector& row0,
    const R3Vector& column0,
    const R3Point& leftTop1,
    const R3Vector& row1,
    const R3Vector& column1,
    R3Point& p0, R3Point& p1
);

//
// Vectors v1, v2, axis must be in the same plane
//
inline double SphereSectorVolume(
    const R3Vector& axis,
    const R3Vector& v1,
    const R3Vector& v2,
    double alpha
)
{
    R3Vector vy = v1.VectorProduct(v2);
    double len = vy.Length();
    if (len <= R3_EPSILON)
        return 0.;

    ASSERT(vy.VectorProduct(axis).Length() <= R3_EPSILON);

    vy *= (1./len);
    double s = len / 2.;

    R3Vector vx = vy.VectorProduct(axis);
    vx.Normalize();
    double r1 = v1.ScalarProduct(vx);
    double r2 = v2.ScalarProduct(vx);
    return (
        s * (r1 + r2) * alpha / 3.
    );
}

// Matrix of arbitrary size
class Matrix {
public:
    int rows;
    int cols;
    double* elements;

    void fill(double v) {
        int n = rows * cols;
        for (int i = 0; i < n; ++i) {
            elements[i] = v;
        }
    }

    void erase() {
        fill(0.);
    }

    Matrix():
        rows(1),
        cols(1),
        elements(new double[1])
    {
        elements[0] = 1.;
    }

    Matrix(int n):
        rows(n),
        cols(n),
        elements(new double[n*n])
    {
        erase();
    }

    Matrix(int r, int c):
        rows(r),
        cols(c),
        elements(new double[r*c])
    {
        erase();
    }

    Matrix(const Matrix& m):
        rows(m.rows),
        cols(m.cols),
        elements(new double[rows * cols])
    {
        int n = rows * cols;
        for (int i = 0; i < n; ++i) {
            elements[i] = m.elements[i];
        }
    }

    Matrix& operator=(const Matrix& m) {
        int n = m.rows * m.cols;
        if (n != rows * cols) {
            delete[] elements;
            rows = m.rows;
            cols = m.cols;
            elements = new double[rows * cols];
        }
        for (int i = 0; i < n; ++i) {
            elements[i] = m.elements[i];
        }
        return *this;
    }

    ~Matrix() { delete[] elements; }

    void setDimensions(int numRows, int numCols);

    double& at(int i, int j) {
        ASSERT(0 <= i && i < rows);
        ASSERT(0 <= j && j < cols);
        return elements[i*cols + j];
    }

    const double& at(int i, int j) const {
        ASSERT(0 <= i && i < rows);
        ASSERT(0 <= j && j < cols);
        return elements[i*cols + j];
    }

    double* operator[](int i) {
        ASSERT(0 <= i && i < rows);
        return (elements + i*cols);
    }

    const double* operator[](int i) const {
        ASSERT(0 <= i && i < rows);
        return (elements + i*cols);
    }

    void setUnit() {
        erase();
        int n = rows;
        if (n > cols)
            n = cols;
        for (int i = 0; i < n; ++i) {
            at(i, i) = 1.;
        }
    }

    void GaussMethod(
        int& rank, double& det, double eps = 1e-7
    );

    void printMatrix() const;

    Matrix& mirrorX();

    void zoom(double coeff, Matrix& zoomedMatrix) const;
    void resize(Matrix& zoomedMatrix) const;    // stretch / shrink
    void stretch(Matrix& stretchedMatrix) const;
    void shrink(Matrix& shrinkedMatrix) const;
};

class R3Box {
public:
    R3Point origin;
    R3Vector size;

    R3Box():
        origin(),
        size()
    {}

    R3Box(const R3Point& o, const R3Vector& s):
        origin(o),
        size(s)
    {}

    R3Box(const R3Point& o):    // Zero size box in point o
        origin(o),
        size()
    {}

    R3Box& operator=(const R3Point& o) {
        size = R3Vector(0., 0., 0.);
        origin = o;
        return *this;
    }

    R3Point& getOrigin() { return origin; }
    const R3Point& getOrigin() const { return origin; }
    R3Vector& getSize() { return size; }
    const R3Vector& getSize() const { return size; }

    void setOrigin(const R3Point& o) { origin = o; }
    void setMinimalCorner(const R3Point& o) { setOrigin(o); }
    void setMaximalCorner(const R3Point& c) {
        size = c - origin;
    }

    R3Box& shift(const R3Vector& v) {
        origin += v;
        return *this;
    }

    R3Box& extend(const R3Vector& v) {
        size += v;
        return *this;
    }

    R3Box& extendToPoint(const R3Point& p) {
        if (p.x < origin.x)
            origin.x = p.x;
        if (p.y < origin.y)
            origin.y = p.y;
        if (p.z < origin.z)
            origin.z = p.z;
        R3Vector v = p - origin;
        if (size.x < v.x)
            size.x = v.x;
        if (size.y < v.y)
            size.y = v.y;
        if (size.z < v.z)
            size.z = v.z;
        return *this;
    }

    bool contains(const R3Point& p) const {
        R3Vector v = p - origin;
        return (
            v.x >= 0. && v.y >= 0. && v.z >= 0. &&
            v.x <= size.x && v.y <= size.y && v.z <= size.z
        );
    }

    double aspect() const {
        double minSize = size.x;
        double maxSize = size.x;
        if (size.y < minSize)
            minSize = size.y;
        else if (size.y > maxSize)
            maxSize = size.y;

        if (size.z < minSize)
            minSize = size.z;
        else if (size.z > maxSize)
            maxSize = size.z;

        if (minSize > 0.)
            return (maxSize / minSize);
        else
            return 1e+10;   // Very large number
    }
};

inline double SignedDistanceToPlane(
    const R3Point& p,
    const R3Point& planeOrigin,
    const R3Vector& planeNormal     // Normalized
) {
    return (p - planeOrigin).ScalarProduct(planeNormal);
}

inline double DistanceToPlane(
    const R3Point& p,
    const R3Point& planeOrigin,
    const R3Vector& planeNormal     // Normalized
) {
    return fabs(SignedDistanceToPlane(p, planeOrigin, planeNormal));
}

bool CheckEllipseBy8Rays(
    const double* rayDistances,     // Array of size 8
    double dx, double dy,           // Size of a pixel
    double& minAxis,
    double& maxAxis
);

inline void ConvertR3Coordinates(
    const R3Point& srcOrigin,
    const R3Vector& srcEX,
    const R3Vector& srcEY,
    const R3Vector& srcEZ,

    const R3Point& srcPoint,

    const R3Point& dstOrigin,
    const R3Vector& dstEX,
    const R3Vector& dstEY,
    const R3Vector& dstEZ,

    R3Point& dstPoint
) {
    // Coordinate systems are orthonormal
    R3Vector v =
        (srcOrigin +
        srcEX*srcPoint.x +
        srcEY*srcPoint.y +
        srcEZ*srcPoint.z) -
        dstOrigin;
    dstPoint.x = v*dstEX;
    dstPoint.y = v*dstEY;
    dstPoint.z = v*dstEZ;
}

inline void defineR3Coordinates(
    const R3Point& srcPoint,

    const R3Point& origin,
    const R3Vector& ex,
    const R3Vector& ey,
    const R3Vector& ez,

    R3Point& dstPoint
) {
    // Coordinate system is orthonormal
    R3Vector v = srcPoint - origin;
    dstPoint.x = v*ex;
    dstPoint.y = v*ey;
    dstPoint.z = v*ez;
}

inline void mapR3PointToPlane(
    const R3Point& srcPoint,

    const R3Point& origin,
    const R3Vector& ex,
    const R3Vector& ey,

    R2Point& dstPoint
) {
    // Coordinates system is orthonormal
    R3Vector v = srcPoint - origin;
    dstPoint.x = v*ex;
    dstPoint.y = v*ey;
}

// For a point inside a triangle,
// calculate its bariocentric coordinates.
// Return true iff a point is inside a triangle
bool computeBariocentricCoordinates(
    const R2Point triangle[3],
    const R2Point& p,
    double coords[3]
);

#endif
