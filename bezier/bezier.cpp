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


/***********************************************
 * definities van bezier class functies
 **********************************************/
Bezier::Bezier(unsigned long _aantal) {
  knooppunten = 0;
  id = 1;
  minX = 1000000.0;
  maxX = -1000000.0;
  minY = 1000000.0;
  maxY = -1000000.0;
  aantal = _aantal;
  starTrail =  new trail(aantal);
  textTrail =  new trail(aantal);
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
// int Bezier::calcTextTrail(){
//   for(unsigned long i)
// }
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


int Bezier::divideToLine(int aantal){
  const double timeSeg = textLengte / aantal;
  double timeSegTMP = timeSeg;
  int n = 0;
  double t;
  for (it = segments.begin(); it != segments.end(); it++) {
    t = 0.0;
    while ((t = it->t_atLengthD(t, timeSegTMP)) > 0) {
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

trail::trail(long _aantal) {
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
