///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
#include "TCMatrix.h"

ClassImp(TCMatrix)

TCMatrix::TCMatrix(int i, int j): _re(i,j), _im(i,j)
{}

void
TCMatrix::set(int i, int j, complex<double> c){
  _re[i][j]=c.real();
  _im[i][j]=c.imag();
}

TComplex
TCMatrix::get(int i, int j){
  TComplex c(_re[i][j],_im[i][j]);
  return c;
}

TComplex
TCMatrix::operator() (int i, int j){
  return this->get(i,j);
}
