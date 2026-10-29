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
//! \file ClosedFormContactsInfra.h
//! \brief Shared return type (ContactResult) for every closed-form contact-detection function in
//! this split -- renamed from PaperAnalyticalContacts.h / namespace paper_analytical 2026-10-29.
//! See the sibling files: ClosedFormContactsBoxPoint.h, ClosedFormContactsCylinderPlaneGeometry.h,
//! ClosedFormContactsCylinderCylinder.h, all of which include this one.

#pragma once

#include <cstdlib>
#include <mesa_pd/data/DataTypes.h>

#include <core/math/Constants.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace closed_form_contacts
{
using mesa_pd::Mat3;
using mesa_pd::Vec3;

//! Common return type for every function in this file. Overlap follows the PAPERS' OWN convention
//! (positive = penetrating), which is the OPPOSITE sign from SolidCylinderContactDetection's
//! penetrationDepth (negative = penetrating) -- callers/tests comparing the two must negate one of
//! them. Kept this way deliberately so the numbers here can be checked directly against the papers'
//! own equations without an extra sign translation layer while reading this file.
struct ContactResult
{
   bool        contact  = false;
   real_t      overlap  = real_t(0);
   Vec3        normal   = Vec3(real_t(0));
   Vec3        point    = Vec3(real_t(0));
   const char* scenario = "none";
};

} // namespace closed_form_contacts
} // namespace antidunes
} // namespace walberla
