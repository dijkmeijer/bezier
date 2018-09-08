/**
 * @file test_creationclass.h
 */

/*****************************************************************************
**
**   class bezier is gebaseerd op dl_creationadapter van de dxflib
**   bezier::addControlPoint laad een dxf file met een B-spline:
**   een curve bestaande uit cubische beziercurvesegmenten (4 controlPoints).
**   Elk punt van de segmenten wordt opgeslagen in een <list> obect als
**   bezierSegment.
**   met bezier splitOnT kan elk segment in kleinere beziercurven worden
**   verdeeld om de afwijking na vervormen (op een startrail) te minimaliseren.
**   splitOnD = idem op afstand D.
**   met normalize() wordt de curve veschaald naar lengte = 1 en veplaats met
**   het middelpunt in de oorsprong.
**   normaize(double) = idem met lengte = double.
**
**   bezierSegment = hulp classe voor bewerkingen op een B-spline segmenten.
**
**   class trail heeft de punten van de startrail (moet worden berekend
**   met novas)
**
**   na vervormen van de B-spline met de trail kan de staptijd worden per
**   as worden berekend. t(n+1) op x+1 of x-1 => T = lengte (bezier(t(n+1)-t)
**   nb lengte beziercurve is maat voor tijd.
**
******************************************************************************/

#ifndef BEZIER_H
#define BEZIER_H

#include "dl_creationadapter.h"
#include <iterator>
#include <list>
using namespace std;

/**
 * This class takes care of the entities read from the file.
 * Usually such a class would probably store the entities.
 * this one just prints some information about them to stdout.
 *
 * @author Andrew Mustun
 */

#define RCPLUS 1 // positive richtincoefficient, stijgend
#define RCMIN 2  // negatieve richtincoefficient, dalend
#define RCBOL 3  // richtings verandering van stijgend naar dalend
#define RCHOL 4  // richtings verandering van dalend naar stijgend

#define X 1
#define Y 2

#define SEGMENTDEEL 10000 // aantal delen per segement om de lengte te berekenen
extern double stepSize;
extern double ph;

class Bezier;
bool inRange(double);

class bezierSegment {
public:
  bezierSegment();
  bezierSegment(Bezier *, double, double, double, double, double, double,
                double, double);
  double x1, x2, x3, x4;
  double y1, y2, y3, y4;
  Bezier *B;
  int id;
  double lengte();
  double next(int);                      // start x(n), dx => return t;
  double next2(int, double *, double *); // start x(n), dx => return t;
  int transform(double, double, double, double);
  int scale(double, double, double);
  double ttLengte(double, double); // lengte bezier curve tussen twee t's
  int lseg; // ondervedeling voor de berekening van de lengte; default 1000;
  double tDevider; // ondervedeling voor de berkening van t
  int print();
  double type;
  int direction;
  int getDirection(int);
  double d; // beziercoefficient
  double t0;
  int cubicRoots(int, double, double *);
  double root[3];
  int print(int);
  double t_atLengthD(double, double);

private:
  double curveLengte;
  double Xcoeffs[4]; // beziercoefficient
  double Ycoeffs[4]; // beziercoefficient
  int dirX;
  int dirY;
  double tx;
  double px;
  int ix;
  double ty;
  double py;
  int iy;
};

class Bezier : public DL_CreationAdapter {
public:
  Bezier();
  int knooppunten;
  double lengte;
  double calcLengte();
  int id;
  double sx; // stapgrootte in de X richting
  double sy; // stapgrootte in de y richting
  double stepSize;
  double ph;
  double minX;
  double maxX;
  double minY;
  double maxY;
  double midX;
  double midY;
  double textLengte;
  list<bezierSegment> segments;
  bezierSegment segment;
  bezierSegment segmentToT; // shegment tot t
  double segLen;
  int segmentDirections(); // bereken de richting van de segmenten
  void addControlPoint(const DL_ControlPointData &data);
  int correctCurve(); // p4(i) = p1(i-1) tussen p2(i)en p3(i-1)
  int normalize();
  int normalize(double);
  int splitOnT(int, double); // splits segment (i) op t, voeg deelsegment aan
                             // list en verander restsegment.
  int devideOnT(int, int);
  int splitOnD(int, double); // splits segment (i) op afstand, voeg deelsegment
  int setStepSize(double);
  int calcIntervals(int);
  int divideToLine(int);
  // aan list en verander restsegment.
private:
  list<bezierSegment>::iterator it;
};

class trail {
public:
  double *x;
  double *y;
  int aantal;
  int add(double, double);
  int end;
  trail();
  trail(int);
  ~trail();
};

double map(double, double, double, double, double);
bool almostZero(double);

#endif
