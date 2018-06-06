// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief pair of ints for definition of a cell
struct DigitizationCell
{
  // identification and data
  size_t channel0 = 0;
  size_t channel1 = 1;
  float  data     = 0.;
  // connstruct them
  DigitizationCell(size_t ch0, size_t ch1, float d = 0.)
    : channel0(ch0), channel1(ch1), data(d)
  {
  }
  // copy them
  DigitizationCell(const DigitizationCell& dc)
    : channel0(dc.channel0), channel1(dc.channel1), data(dc.data)
  {
  }
  /// the deposited energy
  /// @param analogueReadout flag indicating if we have analgue readout
  /// @note this function is needed because possible derived classes may
  /// calculate the energy deposit differently. Furthermore this allows to apply
  /// an energy cut, because the energy deposit can also be stored for digital
  /// readout.
  virtual double
  depositedEnergy(bool /*analogueReadout*/) const
  {
    return data;
  }
};

/// @brief DigitizationStep for further handling
struct DigitizationStep
{
  double stepLength;   /// this is the path length within the cell
  double driftLength;  /// this is the path length of the setp center to the
                       /// readout surface
  DigitizationCell stepCell;      /// this is the cell identifier of the segment
  Vector3D         stepEntry;     /// this is the Entry point into the segment
  Vector3D         stepExit;      /// this is the Exit point from the segment
  Vector2D stepReadoutProjected;  /// this is the projected position at the
                                  /// readout surface
  Vector2D stepCellCenter;        /// this is the cell position

  /// Standard constructor
  DigitizationStep()
    : stepLength(0.)
    , driftLength(0.)
    , stepCell(0, 0)
    , stepEntry(0., 0., 0.)
    , stepExit(0., 0., 0.)
    , stepReadoutProjected(0., 0.)
    , stepCellCenter(0., 0.)
  {
  }

  /// Constructor with arguments
  ///
  /// @param sl step length of this step
  /// @param dl drift length of this step
  /// @param dc is the digitization zell (with indices)
  /// @param entryP is the entry position into the cell
  /// @param exitP is the exit position from the cell
  /// @param projectedPosition is the position on the readout surface
  /// @param cellPosition is the nominal position of the cell
  DigitizationStep(double           sl,
                   double           dl,
                   DigitizationCell dc,
                   const Vector3D&  entryP,
                   const Vector3D&  exitP,
                   const Vector2D&  projectedPosition,
                   const Vector2D&  cellPosition)
    : stepLength(sl)
    , driftLength(dl)
    , stepCell(dc)
    , stepEntry(entryP)
    , stepExit(exitP)
    , stepReadoutProjected(projectedPosition)
    , stepCellCenter(cellPosition)
  {
  }
};
}