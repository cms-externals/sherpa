///////////////////////////////////////////////////////////////////////////////
// File: siscone.cpp                                                         //
// Description: source file for the main SISCone class                       //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 219                                                          $//
// $Date:: 2008-04-03 15:21:00 +0200 (Thu, 03 Apr 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "config.h"

#include "ranlux.h"
#include "momentum.h"
#include "defines.h"
#include "siscone.h"
#include "siscone_error.h"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace siscone{
using namespace std;

/***************************************************************
 * Csiscone implementation                                     *
 * final class: gather everything to compute the jet contents. *
 *                                                             *
 * This is the class user should use.                          *
 * It computes the jet contents of a list of particles         *
 * given a cone radius and a threshold for splitting/merging.  *
 ***************************************************************/

// default ctor
//--------------
Csiscone::Csiscone(){
  rerun_allowed = false;
}

// default dtor
//--------------
Csiscone::~Csiscone(){
  rerun_allowed = false;
}

thread_local bool Csiscone::init_done=false;

/*
 * compute the jets from a given particle set doing multiple passes
 * such pass N looks for jets among all particles not put into jets
 * during previous passes.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _f           shared energy threshold for splitting&merging
 *  - _n_pass_max  maximum number of runs
 *  - _ptmin       minimum pT of the protojets
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *    NOTE: using pt leads to IR unsafety for some events with momentum
 *          conservation. So we strongly advise not to change the default
 *          value.
 * return the number of jets found.
 **********************************************************************/
int Csiscone::compute_jets(vector<Cmomentum> &_particles, double _radius, double _f, 
			   int _n_pass_max, double _ptmin,
			   Esplit_merge_scale _split_merge_scale){
  // initialise random number generator
  if (!init_done){
    // initialise random number generator
    ranlux_init();

    // print the banner
    cout << "#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
    cout << "#                    SISCone   version " << setw(28) << left << siscone_version() << "o" << endl;
    cout << "#              http://projects.hepforge.org/siscone                o" << endl;
    cout << "#                                                                  o" << endl;
    cout << "# This is SISCone: the Seedless Infrared Safe Cone Jet Algorithm   o" << endl;
    cout << "# SISCone was written by Gavin Salam and Gregory Soyez             o" << endl;
    cout << "# It is released under the terms of the GNU General Public License o" << endl;
    cout << "#                                                                  o" << endl;
    cout << "# A description of the algorithm is available in the publication   o" << endl;
    cout << "# JHEP 05 (2007) 086 [arXiv:0704.0292 (hep-ph)].                   o" << endl;
    cout << "# Please cite it if you use SISCone.                               o" << endl;
    cout << "#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
    cout << endl;
    // do not do this again
    init_done=true;
  }

  // run some general safety tests (NB: f will be checked in split-merge)
  if (_radius <= 0.0 || _radius >= 0.5*M_PI) {
    ostringstream message;
    message << "Illegal value for cone radius, R = " << _radius 
            << " (legal values are 0<R<pi/2)";
    throw Csiscone_error(message.str());
  }



  ptcomparison.split_merge_scale = _split_merge_scale;
  partial_clear(); // make sure some things are initialised properly

  // init the split_merge algorithm with the initial list of particles
  // this initialises particle list p_left of remaining particles to deal with
  init_particles(_particles);

  bool finished = false;

  rerun_allowed = false;
  protocones_list.clear();

  do{
    // initialise stable_cone finder
    // here we use the list of remaining particles
    // AFTER COLLINEAR CLUSTERING !!!!!!
    Cstable_cones::init(p_uncol_hard);

    // get stable cones
    if (get_stable_cones(_radius)){
      // we have some new protocones.
      // add them to candidates
      protocones_list.push_back(protocones);
      add_protocones(&protocones, R2, _ptmin);
    } else {
      // no new protocone: leave
      finished=true;
    }

    _n_pass_max--;
  } while ((!finished) && (n_left>0) && (_n_pass_max!=0));

  rerun_allowed = true;

  // split & merge
  return perform(_f, _ptmin);
}

/*
 * recompute the jets with a different overlap parameter.
 * we use the same particles and R as in the preceeding call.
 *  - _f           shared energy threshold for splitting&merging
 *  - _ptmin       minimum pT of the protojets
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *    NOTE: using pt leads to IR unsafety for some events with momentum
 *          conservation. So we strongly advise not to change the default
 *          value.
 * return the number of jets found, -1 if recomputation not allowed.
 ********************************************************************/
int Csiscone::recompute_jets(double _f, double _ptmin,
			     Esplit_merge_scale _split_merge_scale){
  if (!rerun_allowed)
    return -1;

  ptcomparison.split_merge_scale = _split_merge_scale;

  // restore particle list
  partial_clear();
  init_pleft();

  // initialise split/merge algorithm
  unsigned int i;
  for (i=0;i<protocones_list.size();i++)
    add_protocones(&(protocones_list[i]), R2, _ptmin);

  // split & merge
  return perform(_f, _ptmin);
}  


// finally, a bunch of functions to access to 
// basic information (package name, version)
//---------------------------------------------

/* 
 * return SISCone package name.
 * This is nothing but "SISCone", it is a replacement to the
 * PACKAGE_NAME string defined in config.h and which is not
 * public by default.
 * return the SISCone name as a string
 */
string siscone_package_name(){
  return VERSION;
}

/* 
 * return SISCone version number.
 * return a string of the form "X.Y.Z" with possible additional tag
 *        (alpha, beta, devel) to mention stability status
 */
string siscone_version(){
  return VERSION;
}

}
