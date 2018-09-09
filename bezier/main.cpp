/*
 * @file main.cpp
 */

/*****************************************************************************
**  $Id: main.cpp 3591 2006-10-18 21:23:25Z andrew $
**
**  This is part of the dxflib library
**  Copyright (C) 2000-2001 Andrew Mustun
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

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "bezier.h"
#include "dl_dxf.h"

void usage();
void ReadingDXF(char *, int);

using namespace std;

int main(int argc, char **argv) {

  // Check given arguments:
  // if (argc < 2) {
  //   usage();
  //   return 0;
  // }
  int out;

  if (argv[1][0] == 'l')
    out = 2;
  else
    out = 1;
  ReadingDXF(argv[2], out);
  // testBezier();
  return 0;
}

/*
 * @brief Prints error message if file name not specified as command
 * line argument.
 */
void usage() { cout << "\nUsage: test <DXF file>\n\n"; }

void ReadingDXF(char *file, int out) {
  // Load DXF file into memory:
  //  cout << "Reading file " << file << "...\n";
  Bezier *bezier = new Bezier(10000);
  DL_Dxf *dxf = new DL_Dxf();
  if (!dxf->in(file, bezier)) { // if file open failed
    cerr << file << " could not be opened.\n";
    return;
  }

  bezier->normalize(1500.0);

  list<bezierSegment>::iterator it;

  if (out == 1)
    cout << "Totale lengte = " << bezier->calcLengte() << endl;

  //  bereken de trail
  double xLengte = bezier->minX - bezier->maxX;
  double xFac = xLengte / bezier->aantal;
  double halfAantal = bezier->aantal / 2.0;
  for (int i = 0; i < bezier->aantal; i++) {
    bezier->starTrail->add((i - halfAantal) * xFac, (i - halfAantal) * xFac);
  }

  
  if (out == 2) {

    for (it = bezier->segments.begin(); it != bezier->segments.end(); it++) {
      cout << setw(7) << it->x1 << " " << setw(7) << it->y1 - 0.1 << "   ";
      cout << setw(7) << it->x2 << " " << setw(7) << it->y2 - 0.1 << "   ";
      cout << setw(7) << it->x3 << " " << setw(7) << it->y3 - 0.1 << "   ";
      cout << setw(7) << it->x4 << " " << setw(7) << it->y4 - 0.1 << endl;
    }
  }
  if (out == 1) {
    bezier->setStepSize(4.03723587505742E-05);
    const double timeSeg = bezier->textLengte / bezier->aantal;

    cout << "\naantal lijnen = " << bezier->divideToLine(bezier->aantal) << endl;
  }
  delete dxf;
  delete bezier;
}
