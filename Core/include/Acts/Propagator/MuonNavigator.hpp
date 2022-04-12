// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

//#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <bitset>

namespace Acts {

/// A Material Collector struct
struct MuonNavigator {

  mutable bool writeGeometry = true;
  struct this_result {
    bool objectHit;
  };

  using result_type = this_result;



  float stationPhi(int iphi, int offset=0) const
  {
    float offsetPhi = offset ? M_PI / 8. : 0.;
    return  - (M_PI / 2. - M_PI / 4. * (iphi - 1) + offsetPhi);
  }



  void stationParameters(std::string station, int I, int Iz) 
  {
    return ;

  }

  enum muStation
  {
    Z, negZ, R, halfWidth, halfLength, thickness, dx, dy, dz, phiOffset, nTypes = 9
  };






  

  /// Collector action for the ActionList of the Propagator
  /// It checks if the state has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the result object to be filled
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    const auto& logger = state.options.logger;
    
    ACTS_VERBOSE("NORA Tihi I'm a new actor");


    auto position = stepper.position(state.stepping);
    float x = position(0); 
    float y = position(1); 
    float z = position(2);

    float r = sqrt(x*x + y*y);
    ACTS_VERBOSE("NORA State position (x, y, z), (r, z): (" << x << ", " << y << ", " << z << "), (" << r << ", " << z << ")");

    /// Toroid coil consist of eight pieces:
    /// two broad sides, two short shides, and four edge pieces.
    /// Each of these three categories have the same size in z
    /// Radius of them all are the same r: 550 mm
    /// Broad: z 11727 mm
    /// Short: z 3505 mm
    /// Edge:  z 1277 mm

    double radius(550.);
    double angle(M_PI / 8);
    double rotate(M_PI / 4);

    double broadz(11727. - radius * std::sin(angle) ); 
    double edgez(1277./2 - radius * std::sin(angle) );

    double shortz(3505./2  - radius * std::sin(angle) - 42.); // 3505 according to the xml, somewhere its wrong and its too long with 42 mm
    double bottomSideR(4720. + radius);
    double topSideR(10030. -  radius );
    double sidez(12630.);

    //Toroids
    std::vector<std::vector< std::shared_ptr<Acts::AbstractVolume> > > coils;
    //BIL :->
    //std::vector< std::shared_ptr<Acts::Layer> > BIL;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BIL;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BIS;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BML;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BMS;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BOL;
    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > BOS;

    std::map<std::string, std::shared_ptr<Acts::AbstractVolume> > MS_BARREL;

    const u_int32_t ncoils = 8;
    float coilAngle = M_PI / 4;
    float startAngle = 22.5 * M_PI / 180.;
    for (u_int32_t i = 0; i < ncoils; i++) 
    {

      std::vector<std::shared_ptr<Acts::AbstractVolume> > coil;
      
      float coilRotate = -(startAngle + coilAngle * i);
      Eigen::AngleAxisd coilPhi(coilRotate, Eigen::Vector3d(0, 0, 1));

      Acts::Translation3 translation{0., topSideR, 0.};
      auto pTransform = coilPhi * Acts::Transform3(translation);
      auto coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, broadz, M_PI, 0., angle, angle);
      auto cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));

      translation = {0., bottomSideR, 0.};
      pTransform = coilPhi * Acts::Transform3(translation);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, broadz, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));

      Vector3 c{0., bottomSideR + edgez * std::sin(rotate) , -(broadz + edgez * std::sin(rotate))};
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(rotate, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, edgez, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));

      c = {0., bottomSideR + edgez * std::sin(rotate), (broadz+edgez * std::sin(rotate))};
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(-rotate, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, edgez, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));

      c = {0., topSideR - edgez * std::sin(rotate), -(broadz + edgez * std::sin(rotate))};
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(rotate + M_PI/2, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, edgez, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));

      c = {0., topSideR - edgez * std::sin(rotate), (broadz + edgez * std::sin(rotate))};
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(rotate-M_PI, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, edgez, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));
   
      c = {0., bottomSideR + shortz + 2 * edgez * std::sin(rotate), -(broadz + 2 * edgez * std::cos(rotate) ) };
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, shortz, M_PI, 0., -angle, -angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));
   
      c = {0., bottomSideR + shortz + 2 * edgez * std::sin(rotate), (broadz + 2 * edgez * std::cos(rotate) ) };
      pTransform =  coilPhi * Eigen::Translation3d(c) * Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d(1, 0, 0));// * Eigen::Translation3d(-c);
      coilCylinderBounds =
          std::make_shared<CylinderVolumeBounds>(0., radius, shortz, M_PI, 0., angle, angle);
      cylinder = std::make_shared<Acts::AbstractVolume>(pTransform, coilCylinderBounds);
      coil.emplace_back(std::move(cylinder));
      
      coils.emplace_back(coil);
    }

    std::map<std::tuple<std::string, int, int, std::byte>, std::vector<float> > mapMS;

    if(writeGeometry)
    {
      ///////////////BIL 4,3 changed from -3251.051 to -3451.051 and -4531.261 to -4331.051 
      mapMS[std::make_tuple("BIL", 7, 1, std::byte{0b11100000}) ] = {       330., -1231.051, 4740.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 3, 1, std::byte{0b00010000}) ] = {       510., -1230.841, 4740.91, 2671.5 / 2,  720.841 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 7, 1, std::byte{0b00001000}) ] = {       560., -1461.051, 5256.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 8,-1, std::byte{0b00000010}) ] = {  -1231.261,       0.0, 4740.91, 2671.5 / 2, 1081.261 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 7, 1, std::byte{0b00000010}) ] = {       330., -1231.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 1, 2, std::byte{0b11111010}) ] = {      1250., -2331.261, 4790.91, 2671.5 / 2, 1081.261 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 5, 3, std::byte{0b11011010}) ] = {      2350., -3251.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 4, 3, std::byte{0b00100000}) ] = {      2350., -3431.261, 4790.91, 2671.5 / 2, 1081.261 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 6, 4, std::byte{0b11111000}) ] = {      3450., -4531.261, 4790.91, 2671.5 / 2, 1081.261 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 7, 4, std::byte{0b00000010}) ] = {      3450., -4351.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 2, 5, std::byte{0b11111000}) ] = {      4550., -5451.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 9, 5, std::byte{0b00000010}) ] = {      4380., -5271.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 4, 6, std::byte{0b11101000}) ] = {      5470., -6551.261, 4790.91, 2671.5 / 2, 1081.261 / 2, 252.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BIL", 5, 6, std::byte{0b00010010}) ] = {      5470., -6371.051, 4790.91, 2671.5 / 2,  901.051 / 2, 252.0, 0., 0., 0., 0. };

      ///////////////BIS
      mapMS[std::make_tuple("BIS", 3, 1, std::byte{0b11110010}) ] = {        10., -1091.261, 4407.93, 1615.0 / 2, 1062.400 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 2, std::byte{0b11110010}) ] = {      1110., -2011.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 3, std::byte{0b11110010}) ] = {      2030., -2931.051, 4407.93, 1615.0 / 2, 1062.400 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 6, 1, std::byte{0b00000100}) ] = {        10., -1091.261, 4407.93, 1615.0 / 2, 1062.400 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 5, 2, std::byte{0b00000100}) ] = {      1110., -2011.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 5, 3, std::byte{0b00000100}) ] = {      2030., -2931.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 3, 1, std::byte{0b00001001}) ] = {        10., -1091.261, 4422.93, 1615.0 / 2, 1062.400 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 2, std::byte{0b00001001}) ] = {      1110., -2011.051, 4422.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 3, std::byte{0b00001001}) ] = {      2030., -2931.051, 4422.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 4, std::byte{0b11100000}) ] = {      2950., -3851.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 5, std::byte{0b11100000}) ] = {      3870., -4771.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 6, std::byte{0b11100000}) ] = {      4790., -5691.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 4, std::byte{0b00001001}) ] = {      3130., -4031.051, 4422.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 5, std::byte{0b00001001}) ] = {      4050., -4951.051, 4422.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 6, std::byte{0b00001001}) ] = {      4970., -5871.051, 4422.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 4, std::byte{0b00010010}) ] = {      2950., -4031.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 5, std::byte{0b00010010}) ] = {      3870., -4951.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 2, 6, std::byte{0b00010010}) ] = {      4790., -5871.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 5, 4, std::byte{0b00000100}) ] = {      3130., -4031.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 5, 5, std::byte{0b00000100}) ] = {      4050., -4951.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 5, 6, std::byte{0b00000100}) ] = {      4970., -5871.051, 4407.93, 1615.0 / 2,  880.300 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS",  7, 7, std::byte{0b10000000}) ] = {    5716.3, -6791.051, 4432.30, 1660.0 / 2, 1162.700 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 13, 7, std::byte{0b01000000}) ] = {    5716.3, -6791.051, 4432.00, 1660.0 / 2, 1162.700 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS",  9, 7, std::byte{0b00100000}) ] = {    5716.3, -6791.051, 4432.00, 1660.0 / 2, 1162.700 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 10, 7, std::byte{0b00010000}) ] = {    5897.5, -6791.051, 4432.30, 1660.0 / 2,  981.500 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS",  8, 7, std::byte{0b00001000}) ] = {    5897.5, -6791.051, 4447.30, 1660.0 / 2,  981.500 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 11, 7, std::byte{0b00000100}) ] = {    5897.5, -6791.051, 4432.30, 1630.0 / 2,  981.500 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 12, 7, std::byte{0b00000010}) ] = {    5897.5, -6791.051, 4432.30, 1630.0 / 2,  981.500 / 2, 142.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BIS", 14, 7, std::byte{0b00000001}) ] = {    5897.5, -6791.051, 4447.00, 1660.0 / 2,  981.500 / 2, 142.0, 0., 0., 0., 1. };

      /////////////////////////BML
      mapMS[std::make_tuple("BML", 2, 1, std::byte{0b11000000}) ] = {       150.,       0.0, 6743.94, 3551.5 / 2, 1441.681 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 5,-1, std::byte{0b10000000}) ] = {  -1831.121,       0.0, 6743.94, 3551.5 / 2,  961.121 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 6,-1, std::byte{0b01000000}) ] = {  -1831.681,       0.0, 6743.94, 3551.5 / 2, 1441.681 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 1, 1, std::byte{0b00100000}) ] = {       390., -1830.000, 6743.94, 3551.5 / 2, 1441.681 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 7, 1, std::byte{0b00011000}) ] = {       870., -1830.000, 6743.94, 3551.5 / 2,  961.121 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 8,-1, std::byte{0b00001000}) ] = {  -1831.401,       0.0, 6743.94, 3551.5 / 2, 1201.401 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 3, 1, std::byte{0b00000101}) ] = {       630., -1830.000, 6743.94, 3551.5 / 2, 1201.401 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 9,-1, std::byte{0b00000010}) ] = {  -1831.961,       0.0, 6743.94, 3551.5 / 2, 1681.961 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 2, 2, std::byte{0b11111111}) ] = {   1850.000, -3530.000, 6743.94, 3551.5 / 2, 1681.961 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 2, 3, std::byte{0b11111111}) ] = {   3550.000, -5230.000, 6743.94, 3551.5 / 2, 1681.961 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 3, 4, std::byte{0b11111101}) ] = {   5250.000, -6450.000, 6743.94, 3551.5 / 2, 1201.401 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 3, 5, std::byte{0b11111101}) ] = {   6470.000, -7670.000, 6743.94, 3551.5 / 2, 1201.401 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 1, 6, std::byte{0b11111101}) ] = {   7690.000, -9130.000, 6743.94, 3551.5 / 2, 1441.681 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 4, 7, std::byte{0b11111101}) ] = {   9180.000, -9660.000, 7430.00, 3551.5 / 2,  961.121 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 3, 4, std::byte{0b00000010}) ] = {   6470.000, -7670.000, 6743.94, 3551.5 / 2, 1201.401 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 1, 5, std::byte{0b00000010}) ] = {   7690.000, -9130.000, 6743.94, 3551.5 / 2, 1441.681 / 2, 387.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BML", 4, 6, std::byte{0b00000010}) ] = {   9180.000, -9660.000, 7430.00, 3551.5 / 2,  961.121 / 2, 387.0, 0., 0., 0., 0. };

      ////////////////////////BMS
      mapMS[std::make_tuple("BMS", 3, 1, std::byte{0b11101001}) ] = {       150.,       0.0, 7768.94, 3071.5 / 2, 1681.961 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 3, 1, std::byte{0b00010000}) ] = {       150.,       0.0, 7768.94, 3071.5 / 2, 1681.961 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 3, 1, std::byte{0b11101001}) ] = {       150., -1830.000, 7768.94, 3071.5 / 2, 1681.961 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 7, 1, std::byte{0b00010000}) ] = {       150.,       0.0, 7768.94, 3071.5 / 2, 1681.961 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 8,-1, std::byte{0b00010000}) ] = {     -1830.,       0.0, 7768.94, 3071.5 / 2, 1681.961 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 5, 2, std::byte{0b11111001}) ] = {      1830., -3645.000, 7768.94, 3071.5 / 2, 1441.681 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 4, 3, std::byte{0b11111001}) ] = {      3680., -5120.000, 7768.94, 3071.5 / 2, 1441.681 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 6, 4, std::byte{0b11111001}) ] = {      5140., -6913.500, 7768.94, 3071.5 / 2, 1441.681 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 1, 5, std::byte{0b11111001}) ] = {      6925., -7885.000, 7768.94, 3071.5 / 2,  961.121 / 2, 414.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BMS", 2, 6, std::byte{0b11111001}) ] = {      7905., -9345.000, 7768.94, 3071.5 / 2, 1441.681 / 2, 414.0, 0., 0., 0., 1. };

      ////////////////////////BOL
      mapMS[std::make_tuple("BOL", 3,-1, std::byte{0b11000000}) ] = {  -2311.681,       0.0, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 1, 1, std::byte{0b11000000}) ] = {       150.,       0.0, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 3,-1, std::byte{0b00100000}) ] = {  -2311.681,       0.0, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. }; // was -3031.681 for Z
    //  mapMS[std::make_tuple("BOL", 3,-2, std::byte{0b00100000}) ] = {  -4491.681,       0.0, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. }; //This one is weird
      mapMS[std::make_tuple("BOL", 4, 1, std::byte{0b00100000}) ] = {       390.,       0.0, 9244.44, 4961.5 / 2, 1922.241 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 6,-1, std::byte{0b00001000}) ] = {  -2311.961,       0.0, 9244.44, 4961.5 / 2, 1681.961 / 2, 321.0, 0., 0., 0., 0. }; //wut
      mapMS[std::make_tuple("BOL", 5, 1, std::byte{0b00011000}) ] = {       870., -2310.000, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. }; 
      mapMS[std::make_tuple("BOL", 7,-1, std::byte{0b00000010}) ] = {  -2312.521,       0.0, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 2, 1, std::byte{0b00000101}) ] = {       630., -2310.000, 9244.44, 4961.5 / 2, 1681.961 / 2, 321.0, 0., 0., 0., 0. }; /// 00000111 to 00000101
      mapMS[std::make_tuple("BOL", 1, 2, std::byte{0b11111101}) ] = {      2330., -4490.521, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 5, 2, std::byte{0b00000010}) ] = {      2330., -3770.000, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 2, 3, std::byte{0b11111101}) ] = {      4510., -6190.000, 9244.44, 4961.5 / 2, 1681.961 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 5, 3, std::byte{0b00000010}) ] = {      3790., -5230.000, 9244.44, 4961.5 / 2, 1441.681 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 1, 4, std::byte{0b11111101}) ] = {      6210., -8370.521, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 1, 4, std::byte{0b00000010}) ] = {      6250., -8410.521, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 1, 5, std::byte{0b11111101}) ] = {      8390.,-10550.521, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 1, 5, std::byte{0b00000010}) ] = {      8430.,-10590.521, 9244.44, 4961.5 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 2, 6, std::byte{0b11111101}) ] = {     10570.,-12250.000, 9244.44, 4961.5 / 2, 1681.961 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 8, 7, std::byte{0b00000010}) ] = {      7205.,       0.0,12574.50, 3773.3 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 8,-7, std::byte{0b00000010}) ] = {   -9202.52,       0.0,12834.50, 3773.3 / 2, 2162.521 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 9, 8, std::byte{0b00000010}) ] = {      7205.,       0.0,13235.50, 4060.0 / 2, 1300.000 / 2, 321.0, 0., 0., 0., 0. };
      mapMS[std::make_tuple("BOL", 9,-8, std::byte{0b00000010}) ] = {   -8520.00,       0.0,13395.50, 4060.0 / 2, 1300.000 / 2, 321.0, 0., 0., 0., 0. };

      /////////////////////////BOS
      mapMS[std::make_tuple("BOS", 2, 1, std::byte{0b11101001}) ] = {        10., -2170.531,10183.94, 3773.3 / 2, 2162.521 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 3, 1, std::byte{0b00010000}) ] = {       730., -2170.000,10183.94, 3773.3 / 2, 1441.681 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 2, 2, std::byte{0b11111001}) ] = {      2190., -4350.521,10183.94, 3773.3 / 2, 2162.521 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 2, 3, std::byte{0b11111001}) ] = {      4370., -6530.521,10183.94, 3773.3 / 2, 2162.521 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 2, 4, std::byte{0b11111001}) ] = {      6550., -8710.521,10183.94, 3773.3 / 2, 2162.521 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 2, 5, std::byte{0b11111001}) ] = {      8730.,-10890.521,10183.94, 3773.3 / 2, 2162.521 / 2, 314.0, 0., 0., 0., 1. };
      mapMS[std::make_tuple("BOS", 1, 6, std::byte{0b11111001}) ] = {     10910.,-12830.241,10183.94, 3773.3 / 2, 1922.241 / 2, 314.0, 0., 0., 0., 1. };


      for( auto station : mapMS ) 
      {
        std::tuple<std::string, int, int, std::byte> stationTuple = station.first;
        std::vector<float> parameters = station.second;
        std::bitset bits = std::bitset<8>(std::to_integer<int>(std::get<3>(stationTuple)));
        for (int i = bits.size() - 1; i >= 0 ; i--) {
          if ( bits.test(i) ) {
            Acts::Translation3 translation{ 0. + parameters[muStation::dx], parameters[muStation::R] + parameters[muStation::thickness], parameters[muStation::Z] + parameters[muStation::halfLength] + parameters[muStation::dz] };
            auto pTransform = Eigen::AngleAxisd(stationPhi(i+1,(int)parameters[muStation::phiOffset]), Eigen::Vector3d(0, 0, 1)) * Acts::Transform3(translation);
            auto pBox = std::make_shared<const Acts::CuboidVolumeBounds>(parameters[muStation::halfWidth], parameters[muStation::thickness], parameters[muStation::halfLength]);
            auto pBoxVolume = std::make_shared<Acts::AbstractVolume>(pTransform, pBox);

            std::string stationName = std::get<0>(stationTuple) + "_I" + std::to_string(std::get<1>(stationTuple)) + "_Iz" + std::to_string(std::get<2>(stationTuple)) + "_Iphi" + std::to_string(i+1);
            MS_BARREL[stationName] = std::move(pBoxVolume);

            if(parameters[muStation::negZ] < 0.0) 
            {
              translation = { 0. + parameters[muStation::dx], parameters[muStation::R] + parameters[muStation::thickness], parameters[muStation::negZ] + parameters[muStation::halfLength] - parameters[muStation::dz] };
              pTransform = Eigen::AngleAxisd(stationPhi(i+1,(int)parameters[muStation::phiOffset]), Eigen::Vector3d(0, 0, 1)) * Acts::Transform3(translation);
              pBox = std::make_shared<const Acts::CuboidVolumeBounds>(parameters[muStation::halfWidth], parameters[muStation::thickness], parameters[muStation::halfLength]);
              pBoxVolume = std::make_shared<Acts::AbstractVolume>(pTransform, pBox);

              stationName = std::get<0>(stationTuple) + "_I" + std::to_string(std::get<1>(stationTuple)) + "_Iz-" + std::to_string(std::get<2>(stationTuple)) + "_Iphi" + std::to_string(i+1);
              MS_BARREL[stationName] = std::move(pBoxVolume);
            }

          } 
        } 
      }
     
      auto gctx = Acts::GeometryContext();
      std::stringstream cStream;
      Acts::ViewConfig vConfig = Acts::s_viewVolume;
      Acts::ViewConfig lConfig = Acts::s_viewLine;
      Acts::ViewConfig sConfig = Acts::s_viewSensitive;
      Acts::ViewConfig gConfig = Acts::s_viewGrid;
      //vConfig.triangulate = 1;

      Acts::ObjVisualization3D helper;
      int coilnumber = 1;
      for (auto coil: coils) 
      {
        int cylinderPart = 1;
        for (auto cylinder: coil)
        {
          GeometryView3D::drawVolume(helper, *cylinder, gctx, Transform3::Identity(),
                                  vConfig);


          helper.write("Coil" + std::to_string(coilnumber) + "_CylinderPart" + std::to_string(cylinderPart) );
          helper.write(cStream);
          helper.clear();
          cylinderPart++;
        }
        coilnumber++;
      }
          
      int stationnumber = 1;

      for (auto station : MS_BARREL)
      {
        GeometryView3D::drawVolume(helper, *station.second, gctx, Transform3::Identity(), vConfig);
        helper.write(station.first);
        helper.write(cStream);
        helper.clear();
        stationnumber++;
      }

      writeGeometry = false;
    }
   /* const CylinderBounds* bounds = &pCylinderLayer->bounds();
    if (bounds->inside3D(position, 0.1))
    {
        ACTS_VERBOSE("NORA I'm inside a cylinder \\o/");
        result.objectHit = true;
    } 
    else 
    {
      result.objectHit = false;
    }*/ 
  }
  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /*state*/) const {}
};  // namespace Acts
}  // namespace Acts
