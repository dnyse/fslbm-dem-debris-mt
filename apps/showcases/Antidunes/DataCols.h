//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file DataCols.h
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#pragma once

#include <map>
#include <string>

namespace walberla {
namespace antidunes {

class DataCols {
public:
   void add(const std::string& key, int val) {
      intCols[key] = val;
   }

   void add(const std::string& key, double val) {
      realCols[key] = val;
   }

   void add(const std::string& key, const std::string& val) {
      strCols[key] = val;
   }

   std::map< std::string, int > getIntCols() const { return intCols; }
   std::map< std::string, std::string > getStrCols() const { return strCols; }
   std::map< std::string, double > getRealCols() const { return realCols; }

private:
   std::map< std::string, int > intCols;
   std::map< std::string, std::string > strCols;
   std::map< std::string, double > realCols;  
};

} // namespace antidunes
} // namespace walberla