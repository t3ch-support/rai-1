/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "../Core/thread.h"
#include "../Geo/geo.h"

#include <math.h>

struct Depth2PointCloud : Thread {
  //inputs
  Var<floatA> depth;
  Var<rai::Transformation> pose;
  //outputs
  Var<arr> points;

  float fx, fy, px, py;
  floatA _depth;
  arr _points;

  Depth2PointCloud(Var<floatA>& _depth, float _fx=NAN, float _fy=NAN, float _px=NAN, float _py=NAN);
  Depth2PointCloud(Var<floatA>& _depth, const arr& FxyCxy);
  virtual ~Depth2PointCloud();

  void open() {}
  void step();
  void close() {}
};

void depthData2point(double* pt, double* FxyCxy);
void depthData2point(arr& pt, const arr& FxyCxy);
void depthData2pointCloud(arr& pts, const floatA& depth, float fx, float fy, float px, float py);
void depthData2pointCloud(arr& pts, const floatA& depth, const arr& FxyCxy);
void point2depthData(double* pt, double* FxyCxy);
