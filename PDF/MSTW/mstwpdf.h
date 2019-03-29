/////////////////////////////////////////////////////////////////////
// Class header file for MSTW PDFs.
// Original C++ version by Jeppe Andersen.
// Modified by Graeme Watt <watt(at)hep.ucl.ac.uk>.
//
// Usage:
// Initialising an instance of the class:
//   c_mstwpdf pdf("filename.dat")
// with "filename.dat" being the data file to interpolate.
//
// A call to the method 
//   c_mstwpdf::update(x,q)
// updates the parton content to the values at (x,q^2).
// The parton contents can then be accessed in
//   c_mstwpdf::cont.upv etc.
// This method is faster if all the flavours need to be evaluated.
//
// A call to the method 
//   c_mstwpdf::parton(f,x,q)
// will return the value of the PDF with flavour 'f' at (x,q^2).
// Now use PDG notation for parton flavour:
//  f =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
//    = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
// Can also get valence quarks directly:
//  f =  7, 8, 9,10,11,12
//    = dv,uv,sv,cv,bv,tv.
// Photon: f = 13.
// This method is faster if only one flavour needs to be evaluated.
/////////////////////////////////////////////////////////////////////

#ifndef _MSTWPDF_H_INCLUDED_
#define _MSTWPDF_H_INCLUDED_
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <math.h>

using namespace std;

class c_mstwpdf {
 private:
  const static int np=12;   // number of parton flavours
  const static int nx=64;   // number of x grid points
  const static int nq=48;   // number of q^2 grid points
  const static int nqc0=4;  // number of q^2 bins below charm threshold
  const static int nqb0=14; // number of q^2 bins below bottom threshold
  static double xx[nx+1];   // grid points in x
  static double qq[nq+1];   // grid points in q^2
  double c[np+1][nx][nq][5][5]; // coefficients used for interpolation
  double parton_interpolate(int flavour,double xxx,double qqq);
  double parton_extrapolate(int flavour,double xxx,double qqq);
  bool warn,fatal;
 public:
  struct s_partoncontent {
    double upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot;
  } cont;
  double xmin,xmax,qsqmin,qsqmax;
  double mCharm,mBottom,alphaSQ0,alphaSMZ;
  double distance,tolerance;
  int alphaSorder,alphaSnfmax;
  void update(double x,double q); // update cont
  double parton(int flavour,double x,double q);
  // The constructor (initialises the functions):
  c_mstwpdf(string filename,bool warn=false,bool fatal=true);
};

#endif
