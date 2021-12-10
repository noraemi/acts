// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

Acts::SurfaceBounds::BoundsType Acts::CylinderBounds::type() const {
  return SurfaceBounds::eCylinder;
}

Acts::Vector2 Acts::CylinderBounds::shifted(
    const Acts::Vector2& lposition) const {
  return {Acts::detail::radian_sym((lposition[Acts::eBoundLoc0] / get(eR)) -
                                   get(eAveragePhi)),
          lposition[Acts::eBoundLoc1]};
}

Acts::ActsMatrix<2, 2> Acts::CylinderBounds::jacobian() const {
  ActsMatrix<2, 2> j;
  j(0, eBoundLoc0) = 1 / get(eR);
  j(0, eBoundLoc1) = 0;
  j(1, eBoundLoc0) = 0;
  j(1, eBoundLoc1) = 1;
  return j;
}

bool Acts::CylinderBounds::inside(const Vector2& lposition,
                                  const BoundaryCheck& bcheck) const {
  return bcheck.transformed(jacobian())
      .isInside(shifted(lposition),
                Vector2(-get(eHalfPhiSector), -get(eHalfLengthZ)),
                Vector2(get(eHalfPhiSector), get(eHalfLengthZ)));
}

bool Acts::CylinderBounds::inside3D(const Vector3& position,
                                    const BoundaryCheck& bcheck) const {
  // additional tolerance from the boundary check if configred
  bool checkAbsolute = bcheck.m_type == BoundaryCheck::Type::eAbsolute;

  // this fast check only applies to closed cylindrical bounds
  double addToleranceR =
      (checkAbsolute && m_closed) ? bcheck.m_tolerance[0] : 0.;
  double addToleranceZ = checkAbsolute ? bcheck.m_tolerance[1] : 0.;
  // check if the position compatible with the radius
  if ((s_onSurfaceTolerance + addToleranceR) <=
      std::abs(perp(position) - get(eR))) {
    return false;
  } else if (checkAbsolute && m_closed) {
    return ((s_onSurfaceTolerance + addToleranceZ + get(eHalfLengthZ)) >=
            std::abs(position.z()));
  }
  // detailed, but slower check
  Vector2 lpos(detail::radian_sym(phi(position) - get(eAveragePhi)),
               position.z());
  return bcheck.transformed(jacobian())
      .isInside(lpos, Vector2(-get(eHalfPhiSector), -get(eHalfLengthZ)),
                Vector2(get(eHalfPhiSector), get(eHalfLengthZ)));
}

std::ostream& Acts::CylinderBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eR) << ", " << get(eHalfLengthZ) << ", ";
  sl << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

std::vector<Acts::Vector3> Acts::CylinderBounds::createCircles(
    const Transform3 ctrans, size_t lseg) const {
  std::vector<Vector3> vertices;

  double avgPhi = get(eAveragePhi);
  double halfPhi = get(eHalfPhiSector);

  bool fullCylinder = coversFullAzimuth();

  // Get the phi segments from the helper - ensures extra points
  auto phiSegs = fullCylinder ? detail::VerticesHelper::phiSegments()
                              : detail::VerticesHelper::phiSegments(
                                    avgPhi - halfPhi, avgPhi + halfPhi,
                                    {static_cast<ActsScalar>(avgPhi)});

  // Write the two bows/circles on either side
  std::vector<int> sides = {-1, 1};
  for (auto& side : sides) {
    for (size_t iseg = 0; iseg < phiSegs.size() - 1; ++iseg) {
      int addon = (iseg == phiSegs.size() - 2 and not fullCylinder) ? 1 : 0;
      /// Helper method to create the segment
      detail::VerticesHelper::createSegment(
          vertices, {get(eR), get(eR)}, phiSegs[iseg], phiSegs[iseg + 1], lseg,
          addon, Vector3(0., 0., side * get(eHalfLengthZ)), ctrans);
    }
  }

  double bevelMinZ = get(eBevelMinZ);
  double bevelMaxZ = get(eBevelMaxZ);

  // Modify the vertices position if bevel is defined
  if ((bevelMinZ != 0. || bevelMaxZ != 0.) && vertices.size() % 2 == 0) {
    auto halfWay = vertices.end() - vertices.size() / 2;
    double mult{1};
    auto func = [&mult](Vector3& v) { v(2) += v(1) * mult; };
    if (bevelMinZ != 0.) {
      mult = std::tan(-bevelMinZ);
      std::for_each(vertices.begin(), halfWay, func);
    }
    if (bevelMaxZ != 0.) {
      mult = std::tan(bevelMaxZ);
      std::for_each(halfWay, vertices.end(), func);
    }
  }
  return vertices;
}
