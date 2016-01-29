/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
 *
 *    Copyright (C) 2016 University of Southern California and
 *                       Meng Zhou
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SINGLE_FAMILY_HMM_HPP
#define SINGLE_FAMILY_HMM_HPP

#include <string>
#include <vector>
#include <fstream>

class SingleFamilyProfileHMM {
public:

  SingleFamilyProfileHMM() : a,b,c {}

  double
  ViterbiDecoding() const;

  double
  BaumWelchTraining() const;

  double
  PosteriorDecoding() const;
}

#endif
