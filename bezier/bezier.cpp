/*
 * @file test_creationclass.cpp
 */

/*****************************************************************************
**  $Id: test_creationclass.cpp 8865 2008-02-04 18:54:02Z andrew $
**
**  This is part of the dxflib library
**  Copyright (C) 2001 Andrew Mustun
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU Library General Public License as
**  published by the Free Software Foundation.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public License
**  along with this program; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/

#include "bezier.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>

using namespace std;

const double oneThird = 1.0 / 3.0;
const double piByThree = M_PI / 3.0;

inline double qCbrt(double x) {
  if (x > 0.0)
    return pow(x, oneThird);
  else if (x < 0.0)
    return -pow(-x, oneThird);
  else
    return 0.0;
}

bezierSegment::bezierSegment() {
  x1 = 0.0;
  y1 = 0.0;
  x2 = 0.0;
  y2 = 0.0;
  x3 = 0.0;
  y3 = 0.0;
  x4 = 0.0;
  y4 = 0.0;
  id = 0;
  lseg = SEGMENTDEEL;
  tx = 0.0;
  ty = 0.0;
  px = x1;
  py = y1;
  t0 = 0.0;
  //  dirX = getDirection(X);
  //  dirY = getDirection(Y);
}

bezierSegment::bezierSegment(Bezier *_B, double _x1, double _y1, double _x2,
                             double _y2, double _x3, double _y3, double _x4,
                             double _y4) {
  x1 = _x1;
  y1 = _y1;
  x2 = _x2;
  y2 = _y2;
  x3 = _x3;
  y3 = _y3;
  x4 = _x4;
  y4 = _y4;
  id = 0;
  lseg = SEGMENTDEEL;
  tx = 0.0;
  ty = 0.0;
  px = x1;
  py = y1;
  t0 = 0.0;
  B = _B;
  Xcoeffs[0] = x4 - x1 + 3.0 * (x2 - x3);
  Xcoeffs[1] = 3.0 * (x1 - 2.0 * x2 + x3);
  Xcoeffs[2] = 3.0 * (x2 - x1);
  Xcoeffs[3] = x1;
  Ycoeffs[0] = y4 - y1 + 3.0 * (y2 - y3);
  Ycoeffs[1] = 3.0 * (y1 - 2.0 * y2 + y3);
  Ycoeffs[2] = 3.0 * (y2 - y1);
  Ycoeffs[3] = y1;
  //  dirX = getDirection(X);
  //  dirY = getDirection(Y);
}

int bezierSegment::print() {
  cout.precision(6);
  cout << fixed << setw(7) << x1 << "," << fixed << setw(7) << y1 << "   ";
  cout << fixed << setw(7) << x2 << "," << fixed << setw(7) << y2 << "   ";
  cout << fixed << setw(7) << x3 << "," << fixed << setw(7) << y3 << "   ";
  cout << fixed << setw(7) << x4 << "," << fixed << setw(7) << y4 << endl;
  return 0;
}

double bezierSegment::lengte() {
  // double l; // resulterende lengte
  // double px1, px2, py1, py2;
  // double t, t2, tn, tn2;
  // px1 = x1;
  // py1 = y1;
  // l = 0.0;
  // for (int i = 1; i <= lseg; i++) {
  //   t = 1.0 * i / lseg;
  //   t2 = t * t;
  //   tn = 1.0 - t;
  //   tn2 = tn * tn;
  //   px2 = (t * t2 * x4) + (3.0 * t2 * tn * x3) + (3.0 * t * tn2 * x2) +
  //         (tn * tn2 * x1);
  //   py2 = (t * t2 * y4) + (3.0 * t2 * tn * y3) + (3.0 * t * tn2 * y2) +
  //         (tn * tn2 * y1);
  //   l += sqrt(((px2 - px1) * (px2 - px1)) + ((py2 - py1) * (py2 - py1)));
  //   px1 = px2;
  //   py1 = py2;
  // }
  // cout << endl;
  curveLengte = ttLengte(0.0, 1.0);
  return curveLengte;
}

// transleer over lijnstuk d1 -> d4, interpoleer d2 en de op 1/3 en 2/3

int bezierSegment::transform(double dx1, double dy1, double dx4, double dy4) {
  double dx2, dy2, dx3, dy3;
  dx2 = (dx4 - dx1) / 3.0;
  dy2 = (dy4 - dy1) / 3.0;
  dx3 = dx1 + 2.0 * dx2;
  dy3 = dy1 + 2.0 * dy2;
  dx2 += dx1;
  dy2 += dy1;
  x1 += dx1;
  y1 += dy1;
  x2 += dx2;
  y2 += dy2;
  x3 += dx3;
  y3 += dy3;
  x4 += dx4;
  y4 += dy4;
  return 0;
}

int bezierSegment::scale(double ox, double oy, double scalefac) {
  x1 = ((x1 - ox) * scalefac);
  y1 = ((y1 - oy) * scalefac);
  x2 = ((x2 - ox) * scalefac);
  y2 = ((y2 - oy) * scalefac);
  x3 = ((x3 - ox) * scalefac);
  y3 = ((y3 - oy) * scalefac);
  x4 = ((x4 - ox) * scalefac);
  y4 = ((y4 - oy) * scalefac);
  px = x1;
  py = y1;
}

/******************************************************************************
/ bereken de t bij de volgende X.X is afhankelijk van de richting van de curve

/ xo = huidige x positie   t = t waarde van punt xo
/ tx is t waarde van nieuwe x coordinaat (x0 +/- stepSize)
******************************************************************************/

// ph is huidige positie
// pvB is volgende positie boven ph
// pvO is volgende positie inder ph

double bezierSegment::next(int axis) {
  double p1, p4;
  int dir;
  if (axis == X) {
    p1 = x1;
    p4 = x4;
    dir = dirX;
  } else {
    p1 = y1;
    p4 = y4;
    dir = dirY;
  }

  switch (dir) {
  case RCMIN: {
    double psize = (p1 - p4) / 7.0;
    for (double p = p4; p <= p1; p += psize) {
      cout << cubicRoots(axis, p, root) << "|";
      if (inRange(root[0]))
        cout << " b " << root[0];
      if (inRange(root[1]))
        cout << " b " << root[1];
      if (inRange(root[2]))
        cout << " b " << root[2];
      cout << endl;
    }
    break;
  }
  case RCPLUS: {
    double psize = (p4 - p1) / 7.0;
    for (double p = p1; p <= p4; p += psize) {
      cout << cubicRoots(axis, p, root) << "|";
      if (inRange(root[0]))
        cout << " a " << root[0];
      if (inRange(root[1]))
        cout << " a " << root[1];
      if (inRange(root[2]))
        cout << " a " << root[2];
      cout << endl;
    }
    break;
  }
  case RCHOL: {
    cout << "**********************************  Hol" << endl;
    break;
  }
  case RCBOL: {
    cout << "**********************************  Bol" << endl;
    break;
  }
  }
  return 0.0;
}

double bezierSegment::next2(int axis, double *cut, double *t) {
  //  double pb = *cut + stepSize
  double stepSize = 4.03723587505742E-05;
  double pb = *cut + stepSize;
  double po = *cut - stepSize;
  int numberOfRoots = cubicRoots(axis, pb, root);
  int n = 0;
  for (int i = 0; i < numberOfRoots; i++) {
    //  cout << i << " r " << root[i] << " ";
    if (inRange(root[i])) {
      cout << "   up root " << i << " " << numberOfRoots << " "
           << " -" << root[i];
      n++;
    }
  }
  numberOfRoots = cubicRoots(axis, po, root);
  n = 0;
  for (int i = 0; i < numberOfRoots; i++) {
    //  cout << i << " r " << root[i] << " ";
    if (inRange(root[i])) {
      cout << " down root " << i << " " << numberOfRoots << " "
           << " -" << root[i];
      n++;
    }
  }
  // if (n == 0)
  // cout << " n= " << n << " roots# " << numberOfRoots;
  return 0;
}

int bezierSegment::cubicRoots(int axis, double cut, double *roots) {
  double a, b, c, d;
  roots[0] = -0.1;
  roots[1] = -0.1;
  roots[2] = -0.1;

  double p1, p2, p3, p4;
  if (axis == X) {
    p1 = x1 - cut;
    p2 = x2 - cut;
    p3 = x3 - cut;
    p4 = x4 - cut;
  } else {
    p1 = y1 - cut;
    p2 = y2 - cut;
    p3 = y3 - cut;
    p4 = y4 - cut;
  }

  a = p4 - p1 + 3.0 * (p2 - p3);
  b = 3.0 * (p1 - 2.0 * p2 + p3);
  c = 3.0 * (p2 - p1);
  d = p1;
  int numberOfRoots = 1;

  // Simple cases with linear, quadratic or invalid equations
  if (almostZero(a)) {
    if (almostZero(b)) {
      if (almostZero(c))
        return 0;
      roots[0] = -d / c;
      numberOfRoots = 1;
      return numberOfRoots;
    }
    const float discriminant = c * c - 4.0 * b * d;
    if (discriminant < 0.0) {
      numberOfRoots = 0;
      return numberOfRoots;
    }

    if (almostZero(discriminant)) {
      roots[0] = -c / (2.0 * b);
      numberOfRoots = 1;
      return numberOfRoots;
    }
    roots[0] = (-c + sqrt(discriminant)) / (2.0 * b);
    roots[1] = (-c - sqrt(discriminant)) / (2.0 * b);
    numberOfRoots = 2;
    return numberOfRoots;
  }
  //  s = -a/3.0 + ()
  // See
  // https://en.wikipedia.org/wiki/Cubic_function#General_solution_to_the_cubic_equation_with_real_coefficients
  // for a description. We depress the general cubic to a form that can more
  // easily be solved. Solve it and then
  // substitue the results back to get the roots of the original cubic.

  // Put cubic into normal format: x^3 + Ax^2 + Bx + C = 0
  const double A = b / a;
  const double B = c / a;
  const double C = d / a;
  // Substitute x = y - A/3 to eliminate quadratic term (depressed form):
  // x^3 + px + q = 0
  const double Asq = A * A;
  const double p = oneThird * (-oneThird * Asq + B);
  const double q = 1.0 / 2.0 * (2.0 / 27.0 * A * Asq - oneThird * A * B + C);
  // Use Cardano's formula
  const double pCubed = p * p * p;
  const double discriminant = q * q + pCubed;
  if (almostZero(discriminant)) {
    if (q == 0.0) {
      // One repeated triple root
      roots[0] = 0.0;
      numberOfRoots = 1;
    } else {
      // One single and one double root
      double u = pow(-q, oneThird);
      roots[0] = 2.0 * u;
      roots[1] = -u;
      numberOfRoots = 2;
    }
  } else if (discriminant < 0) {
    // Three real solutions
    double phi = oneThird * acos(-q / sqrt(-pCubed));
    double t = 2.0 * sqrt(-p);
    roots[0] = t * cos(phi);
    roots[1] = -t * cos(phi + piByThree);
    roots[2] = -t * cos(phi - piByThree);
    numberOfRoots = 3;
  } else {
    // One real solution
    double sqrtDisc = sqrt(discriminant);
    double u = qCbrt(sqrtDisc - q);
    double v = -qCbrt(sqrtDisc + q);
    roots[0] = u + v;
    numberOfRoots = 1;
  }
  // Substitute back in
  const double sub = oneThird * A;
  for (int i = 0; i < numberOfRoots; ++i) {
    roots[i] -= sub;
    // Take care of cases where we are close to zero or one
    if (almostZero(roots[i]))
      roots[i] = 0.f;
    if (almostZero(roots[i] - 1.0))
      roots[i] = 1.0;
  }
  return numberOfRoots;
}

double bezierSegment::ttLengte(double ts, double te) { // t start; t eind
  double l;                                            // resulterende lengte
  double px1, px2, py1, py2;
  double t, t2, tn, tn2;

  t = ts;
  t2 = t * t;
  tn = 1.0 - t;
  tn2 = tn * tn;
  px1 = (t * t2 * x4) + (3.0 * t2 * tn * x3) + (3.0 * t * tn2 * x2) +
        (tn * tn2 * x1);
  py1 = (t * t2 * y4) + (3.0 * t2 * tn * y3) + (3.0 * t * tn2 * y2) +
        (tn * tn2 * y1);
  l = 0.0;
  int tseg = lseg * (te - ts);
  double dt = 1.0 / lseg;
  for (int i = 0; i < tseg; i++) {
    t += dt;
    t2 = t * t;
    tn = 1.0 - t;
    tn2 = tn * tn;
    px2 = (t * t2 * x4) + (3.0 * t2 * tn * x3) + (3.0 * t * tn2 * x2) +
          (tn * tn2 * x1);
    py2 = (t * t2 * y4) + (3.0 * t2 * tn * y3) + (3.0 * t * tn2 * y2) +
          (tn * tn2 * y1);
    l += sqrt(((px2 - px1) * (px2 - px1)) + ((py2 - py1) * (py2 - py1)));
    px1 = px2;
    py1 = py2;
  }
  //  cout << " l = " << l;
  return l;
}

double bezierSegment::t_atLengthD(double ts, double length) {
  double l; // resulterende lengte
  double px1, px2, py1, py2;
  double t, t2, tn, tn2;
  t = ts;
  t2 = t * t;
  tn = 1.0 - t;
  tn2 = tn * tn;
  px1 = (t * t2 * x4) + (3.0 * t2 * tn * x3) + (3.0 * t * tn2 * x2) +
        (tn * tn2 * x1);
  py1 = (t * t2 * y4) + (3.0 * t2 * tn * y3) + (3.0 * t * tn2 * y2) +
        (tn * tn2 * y1);
  l = 0.0;
  double dt = 1.0 / lseg;
  double dl;
  while (l < length) {
    t += dt;
    if (t > 1.0)
      return -ttLengte(ts, 1.0);
    t2 = t * t;
    tn = 1.0 - t;
    tn2 = tn * tn;
    px2 = (t * t2 * x4) + (3.0 * t2 * tn * x3) + (3.0 * t * tn2 * x2) +
          (tn * tn2 * x1);
    py2 = (t * t2 * y4) + (3.0 * t2 * tn * y3) + (3.0 * t * tn2 * y2) +
          (tn * tn2 * y1);

    dl = sqrt(((px2 - px1) * (px2 - px1)) + ((py2 - py1) * (py2 - py1)));
    l += dl;
    px1 = px2;
    py1 = py2;
  }
  t = t-((l-length)/dl*dt);

  return t;
}

int bezierSegment::print(int axis) {
  if (axis == X) {
    cout << setw(7) << x1 << "  " << setw(7) << x2 << "  " << setw(7) << x3
         << "  " << setw(7) << x4 << "  ";
  } else {
    cout << y1 << "  " << y2 << "  " << y3 << "  " << y4 << "  ";
  }
  return 0;
}


int bezierSegment::getDirection(int axis) {
  double p1, p2, p3, p4;
  if (axis == X) {
    p1 = x1;
    p2 = x2;
    p3 = x3;
    p4 = x4;
  } else {
    p1 = y1;
    p2 = y2;
    p3 = y3;
    p4 = y4;
  }
  double pt0;
  double pt1;
  cout.precision(15);
  //  cout << p1 << "  " << p2 << "  " << p3 << "  " << p4 << "  ";
  double a, b, c, d;
  a = (-p1 + (3.0 * p2) - (3.0 * p3) + p4) * 3.0;
  //  cout << a << "  ";
  b = 6.0 * (p1 - 2.0 * p2 + p3);
  //  cout << b << "  ";
  c = 3.0 * (p1 - p2);
  //  cout << c << "  ";

  pt0 = 3.0 * (p2 - p1);
  pt1 = 3.0 * (p4 - p3);
  d = (b * b) - (4 * a * c);
  // cout << d << "  ";
  if (d > 0.0) {
    double sqrtd = sqrt(d);
    double sol1 = (-b - sqrtd) / (2.0 * a);
    double sol2 = (-b + sqrtd) / (2.0 * a);

    if ((sol1 > 0.0) && (sol1 < 1.0))
      if (pt0 > 0.0) {
        direction = RCBOL;
      } else {
        direction = RCHOL;
      }

    if ((sol2 > 0.0) && (sol2 < 1.0))
      if (pt0 > 0.0) {
        direction = RCBOL;
      } else {
        direction = RCHOL;
      }
  } else {
    if (pt0 * pt1 > 0)
      if (pt0 > 0) {
        direction = RCPLUS;
      } else {
        direction = RCMIN;
      }
  }

  return direction;
}
/***********************************************
 * definities van bezier class functies
 **********************************************/
Bezier::Bezier() {
  knooppunten = 0;
  id = 1;
  minX = 1000000.0;
  maxX = -1000000.0;
  minY = 1000000.0;
  maxY = -1000000.0;
  stepSize = 2.2429e-7; // grootte van een stap van de motor.
}

void Bezier::addControlPoint(const DL_ControlPointData &data) {
  if (data.x > maxX)
    maxX = data.x;
  if (data.x < minX)
    minX = data.x;
  if (data.y > maxY)
    maxY = data.y;
  if (data.y < minY)
    minY = data.y;
  if (knooppunten == 0) {
    segment.x1 = data.x;
    segment.y1 = data.y;
    knooppunten++;
  } else {
    if (knooppunten % 3 == 0) {
      segment.x4 = data.x;
      segment.y4 = data.y;
      bezierSegment *seg = new bezierSegment(
          this, segment.x1, segment.y1, segment.x2, segment.y2, segment.x3,
          segment.y3, segment.x4, segment.y4);
      seg->id = id++;
      segments.push_back(*seg);
      segment.x1 = segment.x4;
      segment.y1 = segment.y4;

    } else if (knooppunten % 3 == 1) {
      segment.x2 = data.x;
      segment.y2 = data.y;
    }
    if (knooppunten % 3 == 2) {
      segment.x3 = data.x;
      segment.y3 = data.y;
    }
    knooppunten++;
  }
}

int Bezier::normalize() {
  midX = (maxX + minX) / 2.0;
  midY = (maxY + minY) / 2.0;
  double scale = 1.0 / (maxX - minX);
  list<bezierSegment>::iterator it;
  for (it = segments.begin(); it != segments.end(); it++) {
    it->scale(midX, midY, scale);
  }
  minX = (minX - midX) * scale;
  maxX = (maxX - midX) * scale;
  minY = (minY - midY) * scale;
  maxY = (maxY - midY) * scale;
  midX = 0.0;
  midY = 0.0;
  return 0;
}

int Bezier::normalize(double lengte) {
  midX = (maxX + minX) / 2.0;
  midY = (maxY + minY) / 2.0;
  double scale = lengte / (maxX - minX);
  list<bezierSegment>::iterator it;
  for (it = segments.begin(); it != segments.end(); it++) {
    it->scale(midX, midY, scale);
  }
  minX = (minX - midX) * scale;
  maxX = (maxX - midX) * scale;
  minY = (minY - midY) * scale;
  maxY = (maxY - midY) * scale;
  midX = 0.0;
  midY = 0.0;
  return 0;
}

int Bezier::correctCurve() {
  it = segments.begin();
  list<bezierSegment>::iterator itn = it;
  itn++;
  do {
    it->x4 = (it->x3 + itn->x2) / 2.0;
    itn->x1 = it->x4;
    it->y4 = (it->y3 + itn->y2) / 2.0;
    itn->y1 = it->y4;
    it = itn;
    itn++;
  } while (itn != segments.end());
  return 0;
}

double Bezier::calcLengte() {
  lengte = 0.0;
  for (it = segments.begin(); it != segments.end(); it++) {
    lengte += it->ttLengte(0.0, 1.0);
  }
  textLengte = lengte;
  return lengte;
}

int Bezier::segmentDirections() {
  for (it = segments.begin(); it != segments.end(); it++) {
    it->getDirection(X);
  }
  return 0;
}

int Bezier::devideOnT(
    int n, int aantal) { // n = startpunt in list, aantal = extra punten
  double fac;
  for (int i = 0; i < aantal - 1; i++) {
    fac = 1.0 / (aantal - i);
    splitOnT(i + n, fac);
  }
  knooppunten += (aantal - 1) * 3;
  return aantal;
}

int Bezier::splitOnT(int n, double t) {
  double x2, y2;
  double x3, y3;
  double x4, y4;
  double x5, y5;
  double x6, y6;
  double x7, y7;
  list<bezierSegment>::iterator it;
  it = segments.begin();
  advance(it, n);
  // it->print();

  //  1 - 2 - 3 - 4
  x2 = (it->x2 - it->x1) * t + it->x1;
  y2 = (it->y2 - it->y1) * t + it->y1;
  x7 = (it->x3 - it->x2) * t + it->x2;
  y7 = (it->y3 - it->y2) * t + it->y2;
  x6 = (it->x4 - it->x3) * t + it->x3;
  y6 = (it->y4 - it->y3) * t + it->y3;
  x3 = (x7 - x2) * t + x2;
  y3 = (y7 - y2) * t + y2;
  x5 = (x6 - x7) * t + x7;
  y5 = (y6 - y7) * t + y7;
  x4 = (x5 - x3) * t + x3;
  y4 = (y5 - y3) * t + y3;
  bezierSegment *seg =
      new bezierSegment(this, it->x1, it->y1, x2, y2, x3, y3, x4, y4);
  seg->id = it->id * 1000 + n;
  // seg -> print();
  it->x1 = x4;
  it->y1 = y4;
  it->x2 = x5;
  it->y2 = y5;
  it->x3 = x6;
  it->y3 = y6;
  // it->print();
  segments.insert(it, *seg);

  return 0;
}

int Bezier::splitOnD(int n, double d) { return 0; }

int Bezier::setStepSize(double _stepSize) {
  stepSize = _stepSize;
  return 0;
}

int Bezier::calcIntervals(int axis) {
  int n = 0;
  double tijd = 0.0;
  double segStep;
  it = segments.begin();
  if (axis == X) {
    ph = it->x1;
  } else {
    ph = it->y1;
  }

  for (it = segments.begin(); it != segments.end(); it++) {
    segLen = 0.0;
    while (segLen < 10.0) {
      //  segLen = (it->next(Y))*3600.0/textLengte;
      segLen = (it->next(axis));
      if (segLen < 10.0) {
        tijd += segLen * 3600.0 / textLengte;
        //    cout << int(segLen * 200000000) << endl;
      } else {
        segStep = segLen - 10.0;
        tijd += segStep * 3600.0 / textLengte;
        //    cout << int(segStep * 200000000) << endl;
      }
      n++;
    }
    // cout << "nextX = " << ts << " Nr " << i++ << endl;
    it->t0 = 0;
  }
  cout << "tijd = " << tijd;
  return n;
}

int Bezier::divideToLine(int aantal){
  const double timeSeg = textLengte / aantal;
  double timeSegTMP = timeSeg;
  int n = 0;
  double t;
  for (it = segments.begin(); it != segments.end(); it++) {
    t = 0.0;
    cout << "--> ";
    while ((t = it->t_atLengthD(t, timeSegTMP)) > 0) {
      cout << t << " ";
      timeSegTMP = timeSeg;
      n++;
    }
    timeSegTMP = timeSeg + t;
  }
  return n;
}

/***********************************************
 * definities van trail class functies
 **********************************************/
trail::trail() {
  aantal = 3600;
  x = new double[aantal];
  y = new double[aantal];
  end = 0;
}

trail::trail(int _aantal) {
  aantal = _aantal;
  x = new double[aantal];
  y = new double[aantal];
  end = 0;
}

int trail::add(double _x, double _y){
  if(end >= aantal) return -1;
  x[end] = _x;
  y[end] = _y;
  end++;
  return end;
}

trail::~trail() {
  delete[] x;
  delete[] y;
}

double map(double pc, double p0, double p1, double t0, double t1) {

  return (((pc - p0) / (p1 - p0)) + (t0 / (t1 - t0))) * (t1 - t0);
}

// EOF

bool almostZero(double value) {
  // 1e-3 might seem excessively fuzzy, but any smaller value will make the
  // factors a, b, and c large enough to knock out the cubic solver.
  double threshold = 1e-6;
  return value > -threshold && value < threshold;
}

bool inRange(double x) { return x >= 0.0 && x <= 1.0; }
