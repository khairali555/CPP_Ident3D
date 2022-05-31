// #include "stdafx.h"
#include <stdio.h>
#include "matrix.h"
#include "PI.h"
#include "graph.h"
#include <limits.h>

#ifdef ThreeDWK
#include "mgdilib.h"
#include "mdrawab.h"
#endif

double R2LineSegment::precision = R2_EPSILON;
double R2OrientedLineSegment::precision = R2_EPSILON;

void R2LineSegment::setPrecision(double prec) {
    precision = prec;
}

void R2OrientedLineSegment::setPrecision(double prec) {
    precision = prec;
}

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
) {
    R3Vector v = v1.VectorProduct(v0);
    double len = v.Length();
    if (len == 0.)
        return FALSE;
    v *= (1./len);      // Normalize

    R3Vector n = v1.VectorProduct(v);
    n.Normalize();

    double t = ((p0 - p1) * v0) / (n * v0);

    if (fabs(t) > r)
        return FALSE;

    double l = sqrt(r*r - t*t);
    R3Point q = p1 + n * t;

    q1 = q + v * l;
    q2 = q - v * l;
    return TRUE;
}

/*...
template <> void AFXAPI
    ConstructElements <R3Point> (R3Point* pNew, int nCount)
{
    for (int i = 0; i < nCount; i++, pNew++)
    {
        // call default constructor directly
        new (pNew) R3Point();
    }
}

template <> void AFXAPI
    DestructElements <R3Point> (R3Point* pElements, int nCount)
{
    for (int i = 0; i < nCount; i++, pElements++)
    {
        // Call default constructor directly
        pElements->R3Point::~R3Point();
    }
}
...*/

static const int CYLINDER_MAX_STEPS = 256;
static const double CYLINDER_STEP_FACTOR = 8.;

BOOL IntersectPlaneAndCylinder(             // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cylinder
    double diameter, double height,

    R3PointArray& section
) {
    //MFC section.RemoveAll();
    section.clear();
    R3PointArray mirror;

    double r = diameter / 2.;
    R3Vector axis = v1;
    axis.Normalize();

    // Define a step
    R3Vector normal = v0;
    normal.Normalize();     // Normal to plane

    // Alpha is an angle between a cylinder axis and a plane
    double sinAlpha = axis * normal;
    double cosAlpha = sqrt(1. - sinAlpha * sinAlpha);

    int numSteps;
    if (cosAlpha == 0.)
        numSteps = CYLINDER_MAX_STEPS;
    else
        numSteps = (int) (CYLINDER_STEP_FACTOR * height / cosAlpha);
    if (numSteps > CYLINDER_MAX_STEPS)
        numSteps = CYLINDER_MAX_STEPS;

    double step;
    if (numSteps < 1) {
        numSteps = 1; step = height;
    } else {
        step = height / (double) numSteps;
    }
    // End of step definition

    R3Point p = p1 - axis * (height / 2.);
    R3Vector advance = axis * step;
    int i;
    R3Point q1, q2;
    for (i = 0; i <= numSteps; i++) {
        if (IntersectPlaneAndCircle(
            p0, v0,                     // Plane
            p, axis, r,                 // Circle
            q1, q2
        )) {
            //MFC section.Add(q1); mirror.Add(q2);
            section.push_back(q1); mirror.push_back(q2);
        }
        p += advance;
    }

    //MFC if (section.GetSize() == 0) {
    if (section.size() == 0) {
        return FALSE;       // Intersection is empty
    }

    //MFC for (i = (int) mirror.GetSize() - 1; i >= 0; i--) {
    for (i = (int) mirror.size() - 1; i >= 0; i--) {
        //MFC section.Add(mirror.ElementAt(i));
        section.push_back(mirror.at(i));
    }
    return TRUE;
}

BOOL IntersectPlaneAndCone(         // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cone center, axis
    double upper_diameter, double lower_diameter,
    double height,

    R3PointArray& section
) {
    //MFC section.RemoveAll();
    section.clear();
    R3PointArray mirror;

    R3Vector axis = v1;
    axis.Normalize();

    // Define a step
    R3Vector normal = v0;
    normal.Normalize();     // Normal to plane

    // Alpha is an angle between a cylinder axis and a plane
    double sinAlpha = axis * normal;
    double cosAlpha = sqrt(1. - sinAlpha * sinAlpha);

    int numSteps;
    if (cosAlpha == 0.)
        numSteps = CYLINDER_MAX_STEPS;
    else
        numSteps = (int) (CYLINDER_STEP_FACTOR * height / cosAlpha);
    if (numSteps > CYLINDER_MAX_STEPS)
        numSteps = CYLINDER_MAX_STEPS;

    double step;
    if (numSteps < 1) {
        numSteps = 1; step = height;
    } else {
        step = height / (double) numSteps;
    }
    // End of step definition

    double lowerRadius = lower_diameter / 2.;
    double upperRadius = upper_diameter / 2.;
    double radiusStep = (upperRadius - lowerRadius) /
        (double) numSteps;
    double radius = lowerRadius;

    R3Point p = p1 - axis * (height / 2.);
    R3Vector advance = axis * step;
    int i;
    R3Point q1, q2;
    for (i = 0; i <= numSteps; i++) {
        if (IntersectPlaneAndCircle(
            p0, v0,                     // Plane
            p, axis, radius,            // Circle
            q1, q2
        )) {
            //MFC section.Add(q1); mirror.Add(q2);
            section.push_back(q1); mirror.push_back(q2);
        }
        p += advance;
        radius += radiusStep;
    }

    //MFC if (section.GetSize() == 0) {
    if (section.size() == 0) {
        return FALSE;       // Intersection is empty
    }

    //MFC for (i = (int) mirror.GetSize() - 1; i >= 0; i--) {
    for (i = (int) mirror.size() - 1; i >= 0; i--) {
        //MFC section.Add(mirror.ElementAt(i));
        section.push_back(mirror.at(i));
    }
    return TRUE;
}

// Angle step for ellipse drawing
double ELLIPSE_SIN[ELLIPSE_STEPS];
double ELLIPSE_COS[ELLIPSE_STEPS];
BOOL ellipseSinCosCalculated = FALSE;

void CalculateEllipseSinCos() {
    if (ellipseSinCosCalculated)
        return;
    int i;
    double phi = 0.;
    double dphi = 2. * PI / (double) ELLIPSE_STEPS;
    for (i = 0; i < ELLIPSE_STEPS; i++) {
        ELLIPSE_SIN[i] = sin(phi);
        ELLIPSE_COS[i] = cos(phi);
        phi += dphi;
    }
    ellipseSinCosCalculated = TRUE;
}

BOOL IntersectAxialPlaneAndCylinder(        // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cylinder
    double diameter, double height,

    R3PointArray& section
) {
    //MFC section.RemoveAll();
    section.clear();
    R3Vector vv0 = v0; vv0.Normalize();
    R3Vector vv1 = v1; vv1.Normalize();
    double r = diameter / 2.;

    double cosinus = vv0 * vv1;
    //... if (cosinus == 0.)
    //...     return FALSE;
    if (fabs(cosinus) <= 0.99)  // This value should be adjusted.
    {
        return IntersectPlaneAndCylinder(
            p0, v0,             // Plane
            p1, v1,             // Cylinder
            diameter, height,
            section
        );
    }

    //... double t = (vv0 * (p0 - p1)) / cosinus;
    double t = vv1 * (p0 - p1);
    if (fabs(t) > height / 2.)
        return FALSE;

    R3Point q = p1 + vv1 * t;
    R3Vector e1, e2;
    if (fabs(cosinus) == 1.) {
        // Intersection is a circle
        // find any 2 perpendicular vectors in a plane
        if (vv0.x != 0.) {
            e1.y = 1.; e1.z = 0.;
            e1.x = (-(e1.y * vv0.y + e1.y * vv0.y) / vv0.x);
        } else if (vv0.y != 0.) {
            e1.x = 1.; e1.z = 0.;
            e1.y = (-(e1.x * vv0.x + e1.z * vv0.z) / vv0.y);
        } else {
            ASSERT(vv0.z != 0);
            e1.x = 1.; e1.y = 0.;
            e1.z = (-(e1.x * vv0.x + e1.y * vv0.y) / vv0.z);
        }
        e1.Normalize();
        e2 = vv0.VectorProduct(e1);     // Length(e2) == 1.
    } else {
        e1 = vv0.VectorProduct(vv1);
        ASSERT(e1.Length() > 0.);
        e1.Normalize();
        e2 = vv0.VectorProduct(e1);     // Length(e2) == 1.
    }
    e1 *= r;
    e2 *= (r / cosinus);

    // Draw ellipse
    CalculateEllipseSinCos();
    int i;
    for (i = 0; i < ELLIPSE_STEPS; i++) {
        //MFC section.Add(
        section.push_back(
            q + e1 * ELLIPSE_COS[i] + e2 * ELLIPSE_SIN[i]
        );
    }
    return TRUE;
}

BOOL IntersectAxialPlaneAndCone(        // TRUE, if nonempty
    const R3Point& p0, const R3Vector& v0,  // Plane

    const R3Point& p1, const R3Vector& v1,  // Cone
    double upper_diameter, double lower_diameter,
    double height,

    R3PointArray& section
) {
    //MFC section.RemoveAll();
    section.clear();
    R3Vector vv0 = v0; vv0.Normalize();
    R3Vector vv1 = v1; vv1.Normalize();
    double upper_r = upper_diameter / 2.;
    double lower_r = lower_diameter / 2.;

    double cosinus = vv0 * vv1;
    //... if (cosinus == 0.)
    //...     return FALSE;
    if (fabs(cosinus) <= 0.99)  // This value should be adjusted.
    {
        return IntersectPlaneAndCone(
            p0, v0,             // Plane
            p1, v1,             // Cone
            upper_diameter, lower_diameter, height,
            section
        );
    }

    //... double t = (vv0 * (p0 - p1)) / cosinus;
    double t = vv1 * (p0 - p1);
    if (fabs(t) > height / 2.)
        return FALSE;
    double r =
        (upper_r + lower_r) / 2. +  // radius in the center of a cone
        ((upper_r - lower_r) * t) /
        height;

    R3Point q = p1 + vv1 * t;
    R3Vector e1, e2;
    if (fabs(cosinus) == 1.) {
        // Intersection is a circle
        // find any 2 perpendicular vectors in a plane
        if (vv0.x != 0.) {
            e1.y = 1.; e1.z = 0.;
            e1.x = (-(e1.y * vv0.y + e1.y * vv0.y) / vv0.x);
        } else if (vv0.y != 0.) {
            e1.x = 1.; e1.z = 0.;
            e1.y = (-(e1.x * vv0.x + e1.z * vv0.z) / vv0.y);
        } else {
            ASSERT(vv0.z != 0);
            e1.x = 1.; e1.y = 0.;
            e1.z = (-(e1.x * vv0.x + e1.y * vv0.y) / vv0.z);
        }
        e1.Normalize();
        e2 = vv0.VectorProduct(e1);     // Length(e2) == 1.
    } else {
        e1 = vv0.VectorProduct(vv1);
        ASSERT(e1.Length() > 0.);
        e1.Normalize();
        e2 = vv0.VectorProduct(e1);     // Length(e2) == 1.
    }
    e1 *= r;
    e2 *= (r / cosinus);

    // Draw ellipse
    CalculateEllipseSinCos();
    int i;
    for (i = 0; i < ELLIPSE_STEPS; i++) {
        //MFC section.Add(
        section.push_back(
            q + e1 * ELLIPSE_COS[i] + e2 * ELLIPSE_SIN[i]
        );
    }
    return TRUE;
}

#ifdef ThreeDWK
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

    BOOL lowerHalf /* = FALSE */,
    BOOL upperHalf /* = FALSE */
)
{
    R3Vector planeX(plane_ex); planeX.Normalize();
    R3Vector planeY(plane_ey); planeY.Normalize();
    R3Vector planeZ(plane_ez); planeZ.Normalize();

    R3Vector implAxis(implantAxis);
    implAxis.Normalize();
    BOOL axisDown = FALSE;
    if (implAxis * planeZ < 0.)
    {
        implAxis *= (-1.);
        axisDown = TRUE;
    }

    R3Vector e1 = planeZ.VectorProduct(implAxis);
    e1.Normalize();                     // First axis of ellipse
    R3Vector e2 = implAxis.VectorProduct(e1);
    e2.Normalize();                     // Second axis of ellipse
    if (e2 * planeZ < 0.)
        e2 *= (-1.);

    // Upper ellipse center
    R3Point o1 = implantPosition + implAxis * (height / 2.);

    // Lower ellipse center
    R3Point o2 = implantPosition + implAxis * (-height / 2.);

    R3Point p1 = o1 + e1 * (diameter / 2.);
    R3Point p2 = o1 + e1 * (-diameter / 2.);
    R3Point p3 = o2 + e1 * (diameter / 2.);
    R3Point p4 = o2 + e1 * (-diameter / 2.);

    // Drawing
    // 1. Lower half-ellipse
    RPoint o2_2D(
        ((o2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RVect e1_2D(
        (e1 * planeX) * (diameter / 2.) * zoomX,
        (e1 * planeY) * (diameter / 2.) *zoomY
    );
    RVect e2_2D(
        (e2 * planeX) * (diameter / 2.) * zoomX,
        (e2 * planeY) * (diameter / 2.) * zoomY
    );

    CalculateEllipseSinCos();

    // Draw half-ellipse
    RPoint p(o2_2D + e1_2D);
    pDraw->MoveTo((int) p.x, (int) p.y);
    int i;
    for (i = 1; i < (ELLIPSE_STEPS/2) + 1; i++)
    {
        p = o2_2D + e1_2D * ELLIPSE_COS[i] + e2_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }

    // 2. Upper complete ellipse
    BOOL drawHalf = FALSE;
    if (
        (!axisDown && upperHalf) ||
        (axisDown && lowerHalf)
    )
        drawHalf = TRUE;

    RPoint o1_2D(
        ((o1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    p = o1_2D + e1_2D;
    RPoint p0(p);

    int numSteps = ELLIPSE_STEPS;
    if (drawHalf)
        numSteps = (ELLIPSE_STEPS/2) + 1;

    pDraw->MoveTo((int) p.x, (int) p.y);
    for (i = 1; i < numSteps; i++)
    {
        p = o1_2D + e1_2D * ELLIPSE_COS[i] + e2_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }
    if (!drawHalf)
        pDraw->LineTo((int) p0.x, (int) p0.y);

    // 3. Lateral side of cylinder
    RPoint p1_2D(
        ((p1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p2_2D(
        ((p2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p3_2D(
        ((p3 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p3 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p4_2D(
        ((p4 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p4 - planeOrigin) * planeY) * zoomY + shiftY
    );

    pDraw->MoveTo((int) p3_2D.x, (int) p3_2D.y);
    pDraw->LineTo((int) p1_2D.x, (int) p1_2D.y);

    pDraw->MoveTo((int) p4_2D.x, (int) p4_2D.y);
    pDraw->LineTo((int) p2_2D.x, (int) p2_2D.y);
}

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

    BOOL lowerHalf /* = FALSE */,
    BOOL upperHalf /* = FALSE */
)
{
    double upper_diam = upper_diameter;
    double lower_diam = lower_diameter;

    R3Vector planeX(plane_ex); planeX.Normalize();
    R3Vector planeY(plane_ey); planeY.Normalize();
    R3Vector planeZ(plane_ez); planeZ.Normalize();

    R3Vector implAxis(implantAxis);
    implAxis.Normalize();
    BOOL axisDown = FALSE;
    if (implAxis * planeZ < 0.)
    {
        // Dirty trick!
        implAxis *= (-1.);
        lower_diam = upper_diameter;
        upper_diam = lower_diameter;
        axisDown = TRUE;
    }

    R3Vector e1 = planeZ.VectorProduct(implAxis);
    e1.Normalize();                     // First axis of ellipse
    R3Vector e2 = implAxis.VectorProduct(e1);
    e2.Normalize();                     // Second axis of ellipse
    if (e2 * planeZ < 0.)
        e2 *= (-1.);

    // Upper ellipse center
    R3Point o1 = implantPosition + implAxis * (height / 2.);

    // Lower ellipse center
    R3Point o2 = implantPosition + implAxis * (-height / 2.);

    R3Point p1 = o1 + e1 * (upper_diam / 2.);
    R3Point p2 = o1 + e1 * (-upper_diam / 2.);
    R3Point p3 = o2 + e1 * (lower_diam / 2.);
    R3Point p4 = o2 + e1 * (-lower_diam / 2.);

    // Drawing
    // 1. Lower half-ellipse
    RPoint o2_2D(
        ((o2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RVect e1_lower_2D(
        (e1 * planeX) * (lower_diam / 2.) * zoomX,
        (e1 * planeY) * (lower_diam / 2.) *zoomY
    );
    RVect e2_lower_2D(
        (e2 * planeX) * (lower_diam / 2.) * zoomX,
        (e2 * planeY) * (lower_diam / 2.) * zoomY
    );
    RVect e1_upper_2D(
        (e1 * planeX) * (upper_diam / 2.) * zoomX,
        (e1 * planeY) * (upper_diam / 2.) *zoomY
    );
    RVect e2_upper_2D(
        (e2 * planeX) * (upper_diam / 2.) * zoomX,
        (e2 * planeY) * (upper_diam / 2.) * zoomY
    );

    CalculateEllipseSinCos();

    // Draw half-ellipse
    RPoint p(o2_2D + e1_lower_2D);
    pDraw->MoveTo((int) p.x, (int) p.y);
    int i;
    for (i = 1; i < (ELLIPSE_STEPS/2) + 1; i++)
    {
        p = o2_2D + e1_lower_2D * ELLIPSE_COS[i] +
            e2_lower_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }

    // 2. Upper complete ellipse
    BOOL drawHalf = FALSE;
    if (
        (!axisDown && upperHalf) ||
        (axisDown && lowerHalf)
    )
        drawHalf = TRUE;

    RPoint o1_2D(
        ((o1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    p = o1_2D + e1_upper_2D;
    RPoint p0(p);

    int numSteps = ELLIPSE_STEPS;
    if (drawHalf)
        numSteps = (ELLIPSE_STEPS/2) + 1;

    pDraw->MoveTo((int) p.x, (int) p.y);
    for (i = 1; i < numSteps; i++)
    {
        p = o1_2D + e1_upper_2D * ELLIPSE_COS[i] +
            e2_upper_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }
    if (!drawHalf)
        pDraw->LineTo((int) p0.x, (int) p0.y);

    // 3. Lateral side of cylinder
    RPoint p1_2D(
        ((p1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p2_2D(
        ((p2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p3_2D(
        ((p3 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p3 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p4_2D(
        ((p4 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p4 - planeOrigin) * planeY) * zoomY + shiftY
    );

    pDraw->MoveTo((int) p3_2D.x, (int) p3_2D.y);
    pDraw->LineTo((int) p1_2D.x, (int) p1_2D.y);

    pDraw->MoveTo((int) p4_2D.x, (int) p4_2D.y);
    pDraw->LineTo((int) p2_2D.x, (int) p2_2D.y);
}
#endif // ThreeDWK

// Intersection of plane and ...
bool IntersectPlanes(
    const R3Point& p0, const R3Vector& n0,
    const R3Point& p1, const R3Vector& n1,
    R3Point& p, R3Vector& v
)
{
    if (n0.Parallel(n1))
        return false;
    v = n0.VectorProduct(n1);
    R3Vector u = n0.VectorProduct(v);   // Vector in first plane
    // p = p0 + u*x,     x is a number
    // (p - p1, n1) = 0
    // (p0 + u*x - p1, n1) = 0
    // (p0 - p1, n1) + x*(u, n1) = 0
    // x = (p1 - p0, n1) / (u, n1)
    double x = ((p1 - p0)*n1) / (u*n1);
    p = p0 + (u*x);
    return true;
}

//... Now inline function
//... bool IntersectPlaneAndLine(
//...     const R3Point& p0, const R3Vector& n0,
//...     const R3Point& p1, const R3Vector& v,
//...     R3Point& t
//... )
//... {
//...     if (fabs(n0 * v) <= R3_EPSILON)
//...         return false;
//...     // t = p1 + v*x,    x is a number
//...     // (p1 + v*x - p0, n0) = 0
//...     // x = (p1 - p0, n0) / (v, n0)
//...     double x = ((p1 - p0)*n0) / (v*n0);
//...     t = p1 + (v*x);
//...     return true;
//... }

// Returns number of points in intersection
// (when line segment is in a plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p, const R3Vector& n,
    const R3Point& s0, const R3Point& s1,
    R3Point& t0, R3Point& t1
)
{
    R3Vector v = s1 - s0;
    if (fabs(v*n) <= R3_EPSILON)
    {
        if (
            PointInPlane(s0, p, v) &&
            PointInPlane(s1, p, v)
        )
        {
            t0 = s0; t1 = s1;
            return 2;
        }
        else
        {
            return 0;
        }
    }

    double x0 = (s0 - p)*n;
    double x1 = (s1 - p)*n;

    if (
        (x0 > R3_EPSILON && x1 > R3_EPSILON) ||
        (x0 < (-R3_EPSILON) && x1 < (-R3_EPSILON))
    )
        return 0;

    // t0 = s0 + v*x
    // (t0 - p, n) = 0
    // (s0 + v*x - p, n) = 0
    // x = (p - s0, n) / (v, n)
    double x = ((p - s0)*n) / (v*n);
    t0 = s0 + v*x;
    return 1;
}

// Returns number of points in intersection
// (when line segment is in the plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p,       // Point in a plane (coordinate origin)
    const R3Vector& ex,     // Coordinate system in a plane
    const R3Vector& ey,     // (unit vectors)
    const R3Point& s0, const R3Point& s1,
    R2Point& t0, R2Point& t1    // Points in internal coordinate system
)
{
    R3Vector ez = ex.VectorProduct(ey);
    return IntersectPlaneAndLineSegment(
        p,
        ex, ey, ez,
        s0, s1,
        t0, t1
    );
}

// Returns the number of points in intersection
// (when line segment is in the plane, returns vertices of line segment)
int IntersectPlaneAndLineSegment(
    const R3Point& p,       // Point in a plane (coordinate origin)
    const R3Vector& ex,     // Coordinate system in a plane
    const R3Vector& ey,     // (unit vectors)
    const R3Vector& ez,     // Normal to the plane == ex.VectorProduct(ey)
    const R3Point& s0, const R3Point& s1,
    R2Point& t0, R2Point& t1    // Points in internal coordinate system
)
{
    R3Point p0, p1;
    int res = IntersectPlaneAndLineSegment(
        p, ez,
        s0, s1,
        p0, p1
    );
    if (res >= 1)
    {
        R3Vector v = p0 - p;
        t0.x = v * ex;
        t0.y = v * ey;
    }
    if (res >= 2)
    {
        R3Vector v = p1 - p;
        t1.x = v * ex;
        t1.y = v * ey;
    }
    return res;
}

bool PointInPlane(
    const R3Point& p0,
    const R3Point& p1, const R3Vector& n
)
{
    return (fabs((p0 - p1)*n) <= R3_EPSILON);
}

double R2Vector::Angle(const R2Vector& v) const
{
    if (Length() <= R2_EPSILON || v.Length() <= R2_EPSILON)
        return 0.;
    R2Vector n = Normal();
    return atan2(v * n, v * (*this));
}

double R3Vector::Angle(const R3Vector& v) const
{
    double len = Length();
    if (fabs(len) <= R3_EPSILON)
        return 0.;
    R3Vector e3 = VectorProduct(v);
    double lenProd = e3.Length();
    if (lenProd <= R3_EPSILON)
    {
        double s = (*this) * v;
        if (s >= 0.)
            return 0.;
        else
            return PI;
    }

    e3 *= len / lenProd;
    R3Vector e2 = e3.VectorProduct(*this);
    e2 *= (1. / len);
    // Assert: length(*this) == length(e2) == len
    return atan2(v * e2, v * (*this));
}

R3Matrix::R3Matrix(
    const R3Vector& rotationAxis,
    double alpha
)
{
    if (rotationAxis.Length() <= R3_EPSILON)
    {
        m[0][0] = 1.; m[0][1] = 0.; m[0][2] = 0.;
        m[1][0] = 0.; m[1][1] = 1.; m[1][2] = 0.;
        m[2][0] = 0.; m[2][1] = 0.; m[2][2] = 1.;
        return;
    }

    R3Vector e3 = rotationAxis;
    e3.Normalize();
    R3Vector e1 = rotationAxis.Normal().Normalize();
    R3Vector e2 = e3.VectorProduct(e1);

    // In basis (e1, e2, e3) the rotation matrix is
    // ( cos(alpha)    -sin(alpha)    0 )
    // ( sin(alpha)    cos(alpha)     0 )
    // ( 0             0              1 )

    double cosAlpha = cos(alpha);
    double sinAlpha = sin(alpha);

    R3Vector ex(1., 0., 0.);
    R3Vector ey(0., 1., 0.);
    R3Vector ez(0., 0., 1.);

    double ex_x = ex * e1;
    double ex_y = ex * e2;
    double ex_z = ex * e3;

    // Rotate ex
    double exRot_x = cosAlpha * ex_x - sinAlpha * ex_y;
    double exRot_y = sinAlpha * ex_x + cosAlpha * ex_y;
    double exRot_z = ex_z;

    R3Vector exRot = e1 * exRot_x + e2 * exRot_y + e3 * exRot_z;

    double ey_x = ey * e1;
    double ey_y = ey * e2;
    double ey_z = ey * e3;

    // Rotate ey
    double eyRot_x = cosAlpha * ey_x - sinAlpha * ey_y;
    double eyRot_y = sinAlpha * ey_x + cosAlpha * ey_y;
    double eyRot_z = ey_z;

    R3Vector eyRot = e1 * eyRot_x + e2 * eyRot_y + e3 * eyRot_z;

    double ez_x = ez * e1;
    double ez_y = ez * e2;
    double ez_z = ez * e3;

    // Rotate ey
    double ezRot_x = cosAlpha * ez_x - sinAlpha * ez_y;
    double ezRot_y = sinAlpha * ez_x + cosAlpha * ez_y;
    double ezRot_z = ez_z;

    R3Vector ezRot = e1 * ezRot_x + e2 * ezRot_y + e3 * ezRot_z;

    m[0][0] = exRot.x; m[0][1] = eyRot.x; m[0][2] = ezRot.x;
    m[1][0] = exRot.y; m[1][1] = eyRot.y; m[1][2] = ezRot.y;
    m[2][0] = exRot.z; m[2][1] = eyRot.z; m[2][2] = ezRot.z;
}

R3Matrix R3Matrix::RotationMatrix(
    const R3Vector& rotationAxis,
    double alpha
)
{
    return R3Matrix(rotationAxis, alpha);
}

// Returns true, if for everz axis (x, y, z)
// the angle between new and old axes is less or equal
//
bool R3Matrix::IsRestrictedRotation(double maxAngle) const
{
    R3Vector e(1., 0., 0.);
    R3Vector e1 = (*this) * e;
    if (fabs(e.Angle(e1)) > maxAngle)
        return false;

    e.x = 0.; e.y = 1.;
    e1 = (*this) * e;
    if (fabs(e.Angle(e1)) > maxAngle)
        return false;

    e.y = 0.; e.z = 1.;
    e1 = (*this) * e;
    if (fabs(e.Angle(e1)) > maxAngle)
        return false;

    return true;
}

bool IntersectR3Rectangles(
    const R3Point& leftTop0,
    const R3Vector& row0,
    const R3Vector& column0,
    const R3Point& leftTop1,
    const R3Vector& row1,
    const R3Vector& column1,
    R3Point& p0, R3Point& p1
)
{
    R3Point p;
    R3Vector v;
    if (!IntersectPlanes(
        leftTop0, row0.VectorProduct(column0),
        leftTop1, row1.VectorProduct(column1),
        p, v
    ))
        return false;

    // Consider plane of first rectangle
    R3Vector ex0 = row0; ex0.Normalize();
    R3Vector ey0 = column0; ey0.Normalize();
    double lenX0 = row0.Length();
    double lenY0 = column0.Length();

    RPoint o(0., 0.);
    RVect dx(1., 0.);
    RVect dy(0., 1.);

    RPoint pp = o +
        dx * (ex0 * (p - leftTop0)) +
        dy * (ey0 * (p - leftTop0));
    RVect vv =
        dx * (ex0 * v) +
        dy * (ey0 * v);

    RRect r0(RPoint(0., 0.), lenX0, lenY0);
    RPoint c0, c1;

    if (!Clip(pp, vv, r0, c0, c1))
        return FALSE;

    // Go back to R3 space
    R3Point cc0 = leftTop0 +
        ex0 * c0.x + ey0 * c0.y;
    R3Point cc1 = leftTop0 +
        ex0 * c1.x + ey0 * c1.y;

    // Consider plane of second rectangle
    R3Vector ex1 = row1; ex1.Normalize();
    R3Vector ey1 = column1; ey1.Normalize();
    double lenX1 = row1.Length();
    double lenY1 = column1.Length();
    RPoint ccc0 = o +
        dx * (ex1 * (cc0 - leftTop1)) +
        dy * (ey1 * (cc0 - leftTop1));
    RPoint ccc1 = o +
        dx * (ex1 * (cc1 - leftTop1)) +
        dy * (ey1 * (cc1 - leftTop1));

    RRect r1(RPoint(0., 0.), lenX1, lenY1);
    RPoint cccc0, cccc1;

    if (!Clip(ccc0, ccc1, r1, cccc0, cccc1))
        return FALSE;

    // Finally, go back to R3 space and return a result
    p0 = leftTop1 +
        ex1 * cccc0.x + ey1 * cccc0.y;
    p1 = leftTop1 +
        ex1 * cccc1.x + ey1 * cccc1.y;
    return TRUE;
}

#ifdef ThreeDWK
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
    double coneHeight,          // included in height

    BOOL upperHalf /* = FALSE */
)
{
    R3Vector implAxis(implantAxis);
    implAxis.Normalize();
    R3Vector axis(implAxis);

    double cylHeight = height - coneHeight;
    R3Point cylPosition(
        implantPosition + (axis * (coneHeight/2.))
    );

    R3Vector planeX(plane_ex); planeX.Normalize();
    R3Vector planeY(plane_ey); planeY.Normalize();
    R3Vector planeZ(plane_ez); planeZ.Normalize();

    BOOL axisUp = TRUE;
    if (implAxis * planeZ < 0.)
    {
        implAxis *= (-1.);
        axisUp = FALSE;
    }

    R3Vector e1 = planeZ.VectorProduct(implAxis);
    e1.Normalize();                     // First axis of ellipse
    R3Vector e2 = implAxis.VectorProduct(e1);
    e2.Normalize();                     // Second axis of ellipse
    if (e2 * planeZ < 0.)
        e2 *= (-1.);

    // Upper ellipse center
    R3Point o1 = cylPosition + implAxis * (cylHeight / 2.);

    // Lower ellipse center
    R3Point o2 = cylPosition + implAxis * (-cylHeight / 2.);

    R3Point p1 = o1 + e1 * (diameter / 2.);
    R3Point p2 = o1 + e1 * (-diameter / 2.);
    R3Point p3 = o2 + e1 * (diameter / 2.);
    R3Point p4 = o2 + e1 * (-diameter / 2.);

    // Drawing
    RPoint o2_2D(
        ((o2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RVect e1_2D(
        (e1 * planeX) * (diameter / 2.) * zoomX,
        (e1 * planeY) * (diameter / 2.) *zoomY
    );
    RVect e2_2D(
        (e2 * planeX) * (diameter / 2.) * zoomX,
        (e2 * planeY) * (diameter / 2.) * zoomY
    );

    // Lateral side of cylinder
    RPoint p1_2D(
        ((p1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p2_2D(
        ((p2 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p2 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p3_2D(
        ((p3 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p3 - planeOrigin) * planeY) * zoomY + shiftY
    );
    RPoint p4_2D(
        ((p4 - planeOrigin) * planeX) * zoomX + shiftX,
        ((p4 - planeOrigin) * planeY) * zoomY + shiftY
    );

    // Draw cone
    R3Point coneCenter(o1);
    RPoint coneP1_2D(p1_2D);
    RPoint coneP2_2D(p2_2D);
    if (axisUp)
    {
        coneCenter = o2;
        coneP1_2D = p3_2D;
        coneP2_2D = p4_2D;
    }
    coneCenter -= axis * coneHeight;
    RPoint coneCenter_2D(
        ((coneCenter - planeOrigin) * planeX) * zoomX + shiftX,
        ((coneCenter - planeOrigin) * planeY) * zoomY + shiftY
    );

    RVect coneDH_2D = coneCenter_2D - (
        coneP1_2D +
        ((coneP2_2D - coneP1_2D) * 0.5)
    );

    BOOL coneDrawn = FALSE;
    if (coneDH_2D.Length() >= e2_2D.Length() * 1.1)
        coneDrawn = TRUE;
    // Now, we draw an arrow always!
    pDraw->MoveTo((int) coneP1_2D.x, (int) coneP1_2D.y);
    pDraw->LineTo((int) coneCenter_2D.x, (int) coneCenter_2D.y);
    pDraw->LineTo((int) coneP2_2D.x, (int) coneP2_2D.y);

    CalculateEllipseSinCos();

    // Lower half-ellipse
    RPoint p(o2_2D + e1_2D);
    pDraw->MoveTo((int) p.x, (int) p.y);
    int i;
    for (i = 1; i < (ELLIPSE_STEPS/2) + 1; i++)
    {
        p = o2_2D + e1_2D * ELLIPSE_COS[i] + e2_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }

    // Upper complete ellipse
    RPoint o1_2D(
        ((o1 - planeOrigin) * planeX) * zoomX + shiftX,
        ((o1 - planeOrigin) * planeY) * zoomY + shiftY
    );
    p = o1_2D + e1_2D;
    RPoint p0(p);
    int numSteps = ELLIPSE_STEPS;
    BOOL drawHalf = FALSE;
    if (
        (!axisUp && coneDrawn) ||
        (axisUp && upperHalf)
    )
    {
        // If cone is drawn, draw only a half of ellipse
        numSteps = (ELLIPSE_STEPS/2) + 1;
        drawHalf = TRUE;
    }

    pDraw->MoveTo((int) p.x, (int) p.y);
    for (i = 1; i < numSteps; i++)
    {
        p = o1_2D + e1_2D * ELLIPSE_COS[i] + e2_2D * ELLIPSE_SIN[i];
        pDraw->LineTo((int) p.x, (int) p.y);
    }
    if (!drawHalf)
        pDraw->LineTo((int) p0.x, (int) p0.y);

    pDraw->MoveTo((int) p3_2D.x, (int) p3_2D.y);
    pDraw->LineTo((int) p1_2D.x, (int) p1_2D.y);

    pDraw->MoveTo((int) p4_2D.x, (int) p4_2D.y);
    pDraw->LineTo((int) p2_2D.x, (int) p2_2D.y);
}

#endif /* ThreeDWK */

void Matrix::GaussMethod(
    int& rank, double& det, double eps /* = 1e-7 */
) {
    // printf("Gauss Method:\n");
    // printMatrix();

    int i = 0;
    int j = 0;
    while (i < rows && j < cols) {
        // Invariant: minor (0..i-1, 0..j-1) has a staircase form
        //            minor (i..rows-1, 0..cols-1) is zero

        // Look for the maximal element in column j
        int k = i;
        double m = fabs(at(i, j));
        for (int l = i+1; l < rows; ++l) {
            double t = fabs(at(l, j));
            if (t > m) {
                m = t; k = l;
            }
        }
        if (m <= eps) {
            // Zero column
            for (int l = i; l < rows; ++l) {
                at(l, j) = 0.;
            }
            ++j;
        } else {
            // Swap rows i and k
            if (i != k) {
                for (int l = j; l < cols; ++l) {
                    double t = at(i, l);
                    at(i, l) = at(k, l);
                    at(k, l) = (-t); // Minus to preserve determinant
                }
            }
            // Make zero j-column, starting form element at row i+1
            for (int l = i + 1; l < rows; ++l) {
                double r = at(l, j) / at(i, j);
                at(l, j) = 0.;
                for (int s = j + 1; s < cols; ++s) {
                    at(l, s) -= r * at(i, s);
                }
            }
            ++i; ++j;
        }
    }
    rank = i;
    det = 1.;
    int n = rows;
    if (cols < rows)
        n = cols;
    for (int l = 0; l < n; ++l) {
        det *= at(l, l);
    }

    // printMatrix();
}

void Matrix::printMatrix() const {
    printf("Matrix:\n");
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (j > 0)
                printf("\t");
            printf("%.3lf", at(i, j));
        }
        printf("\n");
    }
}

void R2Matrix::GaussMethod(
    int& rank, double& det, double eps /* = 1e-7 */
) {
    // printf("Gauss Method:\n");
    // printMatrix();

    int i = 0;
    int j = 0;
    int rows = 2;
    int cols = 2;
    while (i < rows && j < cols) {
        // Invariant: minor (0..i-1, 0..j-1) has a staircase form
        //            minor (i..rows-1, 0..cols-1) is zero

        // Look for the maximal element in column j
        int k = i;
        double m = fabs(at(i, j));
        for (int l = i+1; l < rows; ++l) {
            double t = fabs(at(l, j));
            if (t > m) {
                m = t; k = l;
            }
        }
        if (m <= eps) {
            // Zero column
            for (int l = i; l < rows; ++l) {
                at(l, j) = 0.;
            }
            ++j;
        } else {
            // Swap rows i and k
            if (i != k) {
                for (int l = j; l < cols; ++l) {
                    double t = at(i, l);
                    at(i, l) = at(k, l);
                    at(k, l) = (-t); // Minus to preserve determinant
                }
            }
            // Make zero j-column, starting form element at row i+1
            for (int l = i + 1; l < rows; ++l) {
                double r = at(l, j) / at(i, j);
                at(l, j) = 0.;
                for (int s = j + 1; s < cols; ++s) {
                    at(l, s) -= r * at(i, s);
                }
            }
            ++i; ++j;
        }
    }
    rank = i;
    det = 1.;
    int n = rows;
    if (cols < rows)
        n = cols;
    for (int l = 0; l < n; ++l) {
        det *= at(l, l);
    }

    // printMatrix();
}

void R3Matrix::GaussMethod(
    int& rank, double& det, double eps /* = 1e-7 */
) {
    // printf("Gauss Method:\n");
    // printMatrix();

    int i = 0;
    int j = 0;
    int rows = 3;
    int cols = 3;
    while (i < rows && j < cols) {
        // Invariant: minor (0..i-1, 0..j-1) has a staircase form
        //            minor (i..rows-1, 0..cols-1) is zero

        // Look for the maximal element in column j
        int k = i;
        double m = fabs(at(i, j));
        for (int l = i+1; l < rows; ++l) {
            double t = fabs(at(l, j));
            if (t > m) {
                m = t; k = l;
            }
        }
        if (m <= eps) {
            // Zero column
            for (int l = i; l < rows; ++l) {
                at(l, j) = 0.;
            }
            ++j;
        } else {
            // Swap rows i and k
            if (i != k) {
                for (int l = j; l < cols; ++l) {
                    double t = at(i, l);
                    at(i, l) = at(k, l);
                    at(k, l) = (-t); // Minus to preserve determinant
                }
            }
            // Make zero j-column, starting form element at row i+1
            for (int l = i + 1; l < rows; ++l) {
                double r = at(l, j) / at(i, j);
                at(l, j) = 0.;
                for (int s = j + 1; s < cols; ++s) {
                    at(l, s) -= r * at(i, s);
                }
            }
            ++i; ++j;
        }
    }
    rank = i;
    det = 1.;
    int n = rows;
    if (cols < rows)
        n = cols;
    for (int l = 0; l < n; ++l) {
        det *= at(l, l);
    }

    // printMatrix();
}

double R2Point::distanceToLine(
    const R2Point& p, const R2Vector& v
) const {
    R2Vector n = v.normal();
    n.normalize();
    return fabs(((*this) - p)*n);
}

bool IntersectR2Lines(
    const R2Point& p0, const R2Vector& v0,
    const R2Point& p1, const R2Vector& v1,
    R2Point& p,
    double* t0 /* = NULL */, double* t1 /* = NULL */
) {
    if (v0.Parallel(v1)) {
        double d = p0.distanceToLine(p1, v1);
        if (d <= R2_EPSILON) {
            p = p0;
            if (t0 != NULL)
                *t0 = 0.;
            if (t1 != NULL)
                *t1 = 0.;
            return TRUE;
        } else {
            return false;
        }
    }

    R2Vector n = v1.normal();
    // p = p0 + v0 * t;
    // (p - p1) * n == 0
    // ((p0 - p1) + v0 * t) * n == 0
    // (v0 * n) * t == (p1 - p0) * n
    // t == (p1 - p0) * n / (v0 * n)

    double t = ((p1 - p0) * n) / (v0 * n);
    p = p0 + (v0 * t);
    if (t0 != NULL) {
        *t0 = t;
        if (t1 != NULL) {
            R2Vector w = p - p1;
            //... ASSERT(w.parallel(v1));
            *t1 = (v1 * w) / (v1 * v1);
        }
    }
    return true;
}

bool IntersectR2LineAndLineSegment(
    const R2Point& p0, const R2Vector& v0,
    const R2Point& p1, const R2Point& q1,
    R2Point& p
) {
    if (q1 == p1) {
        double d = p1.distanceToLine(p0, v0);
        if (d > R2_EPSILON)
            return false;
        p = p1;
        return true;
    }

    R2Vector v1 = q1 - p1;

    if (v0.Parallel(v1)) {
        double d = p0.distanceToLine(p1, q1);
        if (d > R2_EPSILON) {
            return false;
        }
        p = p1;
        return true;
    }

    R2Point s;
    bool intersects = IntersectR2Lines(
        p0, v0,
        p1, q1 - p1,
        s
    );
    ASSERT(intersects);
    if (!intersects)
        return false;
    if (s.between(p1, q1)) {
        p = s;
        return true;
    }
    return false;
}

bool IntersectR2LineSegments(
    const R2Point& p0, const R2Point& q0,
    const R2Point& p1, const R2Point& q1,
    R2Point& p
) {
    if (p0 == q0) {
        if (p0.between(p1, q1)) {
            p = p0;
            return true;
        }
        return false;
    }
    if (p1 == q1) {
        if (p1.between(p0, q0)) {
            p = p1;
            return true;
        }
        return false;
    }
    R2Vector v0 = (q0 - p0);
    R2Vector v1 = (q1 - p1);
    if (v0.parallel(v1)) {
        double d = p0.distanceToLine(p1, v1);
        if (d <= R2_EPSILON) {
            if (p0.between(p1, q1)) {
                p = p0;
                return true;
            }
            if (q0.between(p1, q1)) {
                p = q0;
                return true;
            }
            if (p1.between(p0, q0)) {
                p = p1;
                return true;
            }
            if (q1.between(p0, q0)) {
                p = q1;
                return true;
            }
            return false;
        }
        return false;
    }

    R2Point s;
    bool intersects = IntersectR2Lines(
        p0, v0,
        p1, v1,
        s
    );
    ASSERT(intersects);
    if (!intersects)
        return false;
    if (
        (s - p0)*(s - q0) <= 0. &&
        (s - p1)*(s - q1) <= 0.
    ) {
        p = s;
        return true;
    }
    return false;
}

// Returns the number of points in intersection with
// the border of triangle
int IntersectPlaneAndTriangle(
    const R3Point& p, const R3Vector& n,
    const R3Point triangeVertices[3],
    R3Point intersection[2]
) {
    int numPoints = 0;
    int k = IntersectPlaneAndLineSegment(
        p, n,
        triangeVertices[0], triangeVertices[1],
        intersection[0], intersection[1]
    );
    if (k == 2) {
        return 2;
    } else if (k == 1) {
        ++numPoints;
    }

    R3Point p0, p1;
    k = IntersectPlaneAndLineSegment(
        p, n,
        triangeVertices[1], triangeVertices[2],
        p0, p1
    );
    if (k == 1) {
        if (numPoints == 0 || p0 != intersection[0]) {
            intersection[numPoints] = p0;
            ++numPoints;
            if (numPoints == 2)
                return 2;
        }
    } else if (k == 2) {
        if (numPoints == 0) {
            intersection[0] = p0;
            intersection[1] = p1;
            return 2;
        } else if (numPoints == 1) {
            // Should not come here
            if (p0 != intersection[0]) {
                intersection[1] = p0;
                return 2;
            } else if (p1 != intersection[0]) {
                intersection[1] = p1;
                return 2;
            } else {
                // All 3 points are equal!
                ASSERT(false);
            }
        }
    }

    k = IntersectPlaneAndLineSegment(
        p, n,
        triangeVertices[2], triangeVertices[0],
        p0, p1
    );
    if (k == 1) {
        if (numPoints == 0 || p0 != intersection[0]) {
            intersection[numPoints] = p0;
            ++numPoints;
            if (numPoints == 2)
                return 2;
        }
    } else if (k == 2) {
        if (numPoints == 0) {
            intersection[0] = p0;
            intersection[1] = p1;
            return 2;
        } else if (numPoints == 1) {
            // Should not come here
            if (p0 != intersection[0]) {
                intersection[1] = p0;
                return 2;
            } else if (p1 != intersection[0]) {
                intersection[1] = p1;
                return 2;
            } else {
                // All 3 points are equal!
                ASSERT(false);
            }
        }
    }

    return numPoints;
}

// Returns the number of points in intersection with
// the border of triangle
int IntersectPlaneAndTriangle(
    const R3Point& p,       // Point in a plane (coordinate origin)
    const R3Vector& ex,     // Coordinate system in a plane
    const R3Vector& ey,     // (unit vectors)
    const R3Vector& ez,     // Normal to the plane == ex.VectorProduct(ey)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
)
{
    R3Point intersect3D[2];
    int res = IntersectPlaneAndTriangle(
        p, ez,
        triangeVertices,
        intersect3D
    );
    if (res >= 1)
    {
        R3Vector v = intersect3D[0] - p;
        intersection[0].x = v * ex;
        intersection[0].y = v * ey;
    }
    if (res >= 2)
    {
        R3Vector v = intersect3D[1] - p;
        intersection[1].x = v * ex;
        intersection[1].y = v * ey;
    }
    return res;
}

int IntersectPlaneAndTriangle(
    const R3Point& p,       // Point in a plane (coordinate origin)
    const R3Vector& ex,     // Coordinate system in a plane
    const R3Vector& ey,     // (unit vectors)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
)
{
    R3Vector ez = ex.VectorProduct(ey);
    return IntersectPlaneAndTriangle(
        p,
        ex, ey, ez,
        triangeVertices,
        intersection
    );
}

int IntersectPlaneAndOrientedTriangle(
    const R3Point& p,           // Point in a plane (coordinate origin)
    const R3Vector& ex,         // Coordinate system in a plane
    const R3Vector& ey,         // (unit vectors)
    const R3Vector& ez,         // Normal to the plane == ex.VectorProduct(ey)
    const R3Point triangeVertices[3],
    R2Point intersection[2]
)
{
    R3Point intersection3D[2];
    int res = IntersectPlaneAndTriangle(
        p, ez,              // Plane
        triangeVertices,    // Triangle
        intersection3D
    );
    if (res == 0)
        return 0;

    if (res == 2 && intersection3D ) {
        // Select the direction of resulting line that corresponds
        // to the orientation of triangle
        R3Vector triangleNormal =
            (triangeVertices[1] - triangeVertices[0]).VectorProduct(
            (triangeVertices[2] - triangeVertices[0])
        );
        R3Vector intersectionDirection =
            ez.VectorProduct(triangleNormal);
        if (
            (intersection3D[1] - intersection3D[0]) *
            intersectionDirection
            >= 0.
        ) {
            // Direction is correct
            R3Vector v = intersection3D[0] - p;
            intersection[0].x = v * ex;
            intersection[0].y = v * ey;

            v = intersection3D[1] - p;
            intersection[1].x = v * ex;
            intersection[1].y = v * ey;
        } else {
            // Swap resulting points
            R3Vector v = intersection3D[1] - p;
            intersection[0].x = v * ex;
            intersection[0].y = v * ey;

            v = intersection3D[0] - p;
            intersection[1].x = v * ex;
            intersection[1].y = v * ey;
        }

#       ifdef _DEBUG
        R3Vector lineDir(
            (p + ex*intersection[1].x + ey*intersection[1].y) -
            (p + ex*intersection[0].x + ey*intersection[0].y)
        );
        double det = R3Matrix(
            lineDir,
            ez,
            triangleNormal
        ).determinant();
        ASSERT(det >= 0.);
#       endif

        return 2;
    } else {
        R3Vector v = intersection3D[0] - p;
        intersection[0].x = v * ex;
        intersection[0].y = v * ey;
        return 1;
    }
}

static const R2Vector RAY_DIRECTIONS[8] = {
    R2Vector(1., 0.),
    R2Vector(SQRT2, SQRT2),
    R2Vector(0., 1.),
    R2Vector(-SQRT2, SQRT2),
    R2Vector(-1., 0.),
    R2Vector(-SQRT2, -SQRT2),
    R2Vector(0., -1.),
    R2Vector(SQRT2, -SQRT2)
};

bool CheckEllipseBy8Rays(
    const double* rayDistances, // Array of size 8
    double dx, double dy,       // Size of a pixel
    double& minAxis,
    double& maxAxis
) {
    // Calculate the gravity center of plain figure
    R2Point vertices[8];
    R2Vector ray_dirs[8];
    const R2Vector* rayDirs;
    int i;
    R2Point o(0., 0.);

    //... (fabs(dx - 1.) <= R2_EPSILON);
    ASSERT(fabs(dx) > R2_EPSILON && fabs(dy) > R2_EPSILON);

    if (fabs(dy - 1.) > R2_EPSILON || fabs(dx - 1.) > R2_EPSILON) {
        // Calculate directions of length 1
        double diag = sqrt(dx*dx + dy*dy);

        ray_dirs[0].x = 1.;         ray_dirs[0].y = 0.;         // (1., 0.)
        ray_dirs[1].x = dx/diag;    ray_dirs[1].y = dy/diag;    // (1., 1.)
        ray_dirs[2].x = 0.;         ray_dirs[2].y = 1.;         // (0., 1.)
        ray_dirs[3].x = (-dx/diag); ray_dirs[3].y = dx/diag;    // (-1., 1.)
        ray_dirs[4].x = (-1.);      ray_dirs[4].y = 0.;         // (-1., 0.)
        ray_dirs[5].x = (-dx/diag); ray_dirs[5].y = (-dy/diag); // (-1., -1.)
        ray_dirs[6].x = 0.;         ray_dirs[6].y = (-1.);      // (0., -1.)
        ray_dirs[7].x = dx/diag;    ray_dirs[7].y = (-dy/diag); // (1., -1.)

        rayDirs = ray_dirs;
    } else {
        rayDirs = RAY_DIRECTIONS;
    }

    for (i = 0; i < 8; ++i) {
        vertices[i] = o +
            rayDirs[i] * rayDistances[i];
    }

    R2Vector shiftToCenter(0., 0.);
    double signedArea = 0.;

    for (i = 1; i < 7; i++)
    {
        R2Vector v0 = vertices[i] - vertices[0];
        R2Vector v1 = vertices[i+1] - vertices[0];
        R2Vector c = (v0 + v1) * (1./3.);

        // Area of triangle
        double a = R2Vector::signed_area(v0, v1);

        shiftToCenter += c*a;
        signedArea += a;
    }

    if (fabs(signedArea) <= R2_EPSILON)
        return false;

    shiftToCenter *= (1./signedArea);
    R2Point center = vertices[0] + shiftToCenter;

    // Define minimal and maximal distances to vertices
    double minD2 = (double) SHRT_MAX;
    double maxD2 = 0.;

    for (i = 0; i < 8; ++i) {
        double dist2 = center.distanceSquare(
            vertices[i]
        );
        if (dist2 < minD2)
            minD2 = dist2;
        if (dist2 > maxD2)
            maxD2 = dist2;
    }
    minAxis = sqrt(minD2);
    maxAxis = sqrt(maxD2);
    return true;
}

void Matrix::setDimensions(int numRows, int numCols) {
    if (
        numRows > rows || numRows < rows / 2 ||
        numCols > cols || numCols < cols / 2
    ) {
        delete[] elements; elements = NULL;
    }
    rows = numRows;
    cols = numCols;
    if (elements == NULL) {
        elements = new double[rows * cols];
    }
    erase();
}

Matrix& Matrix::mirrorX() {
    int y;
    for (y = 0; y < rows; ++y) {
        int x0 = 0; int x1 = cols - 1;
        while (x0 < x1) {
            double tmp = at(y, x0);
            at(y, x0) = at(y, x1);
            at(y, x1) = tmp;
            ++x0; --x1;
        }
    }
    return *this;
}

void Matrix::zoom(double coeff, Matrix& zoomedMatrix) const {
    if (coeff == 1.) {
        zoomedMatrix = *this;
        return;
    }
    int c = (int)((double) cols * coeff + 0.49);
    int r = (int)((double) rows * coeff + 0.49);
    zoomedMatrix.setDimensions(r, c);
    if (coeff >= 1.)
        stretch(zoomedMatrix);
    else
        shrink(zoomedMatrix);
}

void Matrix::resize(Matrix& zoomedMatrix) const {   // stretch / shrink
    if (zoomedMatrix.cols >= cols)
        stretch(zoomedMatrix);
    else
        shrink(zoomedMatrix);
}

//
// Bilinear stretch algorithm
//
void Matrix::stretch(Matrix& m) const {
    ASSERT(rows > 0 && cols > 0);
    ASSERT(m.rows > 0 && m.cols > 0);
    if (
        rows <= 0 || cols <= 0 ||
        m.rows <= 0 || m.cols <= 0
    )
        return;

    double zoomY = (double) m.rows / (double) rows;
    double zoomX = (double) m.cols / (double) cols;
    ASSERT(zoomX >= 1.);
    double dx = 1. / zoomX;
    int x, y;   // coordinates in destination
    y = 0;
    while (y < m.rows) {
        double y_src  = (double) y / zoomY;
        int y0 = (int) y_src;
        int y1 = y0 + 1;
        double w0y = 1. - (y_src - (double) y0);
        double w1y = 1. - w0y;
        if (y1 >= rows) {
            y1 = y0;
            w0y = 1.; w1y = 0.;
        }

        x = 0;
        double x_src = 0.;
        int x0 = 0;
        int x1 = 1;
        if (x1 >= cols)
            x1 = x0;
        double w0x = 1.;
        double w1x = 0.;
        while (x < m.cols) {
            double v0 = at(y0, x0) * w0x + at(y0, x1) * w1x;
            double v1 = at(y1, x0) * w0x + at(y1, x1) * w1x;
            double v = v0 * w0y + v1 * w1y;
            m.at(y, x) = v;

            ++x;
            x_src += dx;
            w0x -= dx;
            w1x += dx;
            if (w0x < 0.) {
                ++x0;
                ++x1;
                w0x += 1.;
                w1x = 1. - w0x;
                if (x1 >= cols) {
                    x1 = x0;
                    w0x = 1.; w1x = 0;
                }
            }
        } // end while (x...

        ++y;
    } // end while (y...
}

void Matrix::shrink(Matrix& m) const {
    if (m.rows == rows && m.cols == cols) {
        m = *this;
        return;
    }

    ASSERT(
        m.rows > 0 && m.cols > 0 &&
        rows > 0 && cols > 0
    );
    ASSERT(m.cols <= cols);
    double zoomX = (double) m.cols / (double) cols;
    double zoomY = (double) m.rows / (double) rows;
    double dx = 1./zoomX;
    double dy = 1./zoomY;

    double y_src = 0.;
    int y_dst;
    for (y_dst = 0; y_dst < m.rows; ++y_dst, y_src += dy) {
        double x_src = 0.;
        int x_dst;
        for (x_dst = 0; x_dst < m.cols; ++x_dst, x_src += dx) {
            double s = 0.;
            double v = 0.;
            double yy0 = y_src;
            while (yy0 < y_src  + dy - R2_EPSILON) {
                double yy1 = ceil(yy0);
                if (fabs(yy1 - yy0) <= R2_EPSILON) { // yy0 is integer
                    yy1 = yy0 + 1.;
                }
                if (yy1 > y_src  + dy)
                    yy1 = y_src  + dy;

                double dyy = yy1 - yy0;
                int yy = (int)(yy0 + R2_EPSILON);

                double xx0 = x_src;
                while (xx0 < x_src  + dx - R2_EPSILON) {
                    double xx1 = ceil(xx0);
                    if (fabs(xx1 - xx0) <= R2_EPSILON) { // xx0 is integer
                        xx1 = xx0 + 1.;
                    }
                    if (xx1 > x_src  + dx)
                        xx1 = x_src  + dx;

                    double dxx = xx1 - xx0;
                    int xx = (int)(xx0 + R2_EPSILON);

                    double ds = dxx * dyy;
                    s += ds;
                    v += at(yy, xx) * ds;

                    xx0 = xx1;
                } // end while (xx0...

                yy0 = yy1;
            } // end while (yy0...

            ASSERT(s > 0.);
            if (fabs(s) > R2_EPSILON)
                v /= s;

            m.at(y_dst, x_dst) = v;

        } // end for (x_dst...
    } // end for (y_dst...
}

// For a point inside a triangle,
// calculate its bariocentric coordinates.
// Return true iff a point is inside a triangle
bool computeBariocentricCoordinates(
    const R2Point triangle[3],
    const R2Point& p,
    double coord[3]
) {
    if (p == triangle[2]) {
        coord[0] = 0.;
        coord[1] = 0.;
        coord[2] = 1.;
        return true;
    }
    R2Vector v = p - triangle[2];
    double t1, t2;

    R2Vector w = triangle[1] - triangle[0];
    R2Point q;

    bool res = IntersectR2Lines(
        triangle[2], v,
        triangle[0], w,
        q,
        &t2, &t1
    );

    if (!res || t2 <= 0.) {
        coord[0] = 0.;
        coord[1] = 0.;
        coord[2] = 1.;
        return (fabs(t2) <= R2_EPSILON);
    }

    if (t2 < 1.)
    {
        res = false;
        t2 = 1.;
    }
    double t12 = 1. / t2;
    coord[2] = 1. - t12;
    coord[0] = (1. - t1) * t12;
    coord[1] = t1 * t12;

    // Test
    ASSERT(
        fabs(coord[0] + coord[1] + coord[2] - 1.) <=
        4. * R2_EPSILON
    );
    ASSERT(
        !res ||
        R2Point::distance(
            p,
            triangle[0] * coord[0] +
            triangle[1] * coord[1] +
            triangle[2] * coord[2]
        ) <= 4. * R2_EPSILON
    );

    return (
        res &&
        -R2_EPSILON <= coord[0] && coord[0] <= 1. + R2_EPSILON &&
        -R2_EPSILON <= coord[1] && coord[1] <= 1. + R2_EPSILON
    );
}
