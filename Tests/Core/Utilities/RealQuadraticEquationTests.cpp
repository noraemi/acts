// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE RealQuadraticEquation Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

//
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"
//
#include <limits>

using Acts::detail::RealQuadraticEquation;

// namespace bdata = boost::unit_test::data;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating RealQuadraticEquation object
  BOOST_AUTO_TEST_CASE(RealQuadraticEquationConstruction)
  {
    double a(1.0), b(-3.), c(2.);
    // test default construction: not deleted, not written
    // RealQuadraticEquation defaultConstructedRealQuadraticEquation;
    //
    /// Test construction with parameters
    BOOST_REQUIRE_NO_THROW(RealQuadraticEquation(a, b, c));
    //
    /// Copy constructor (implicit), void removes 'unused' compiler warning
    RealQuadraticEquation orig(a, b, c);
    BOOST_REQUIRE_NO_THROW(RealQuadraticEquation copied(orig); (void)copied);
  }
  /// Unit test for RealQuadraticEquation properties
  BOOST_AUTO_TEST_CASE(RealQuadraticEquationProperties)
  {
    double a(1.0), b(-3.), c(2.);
    //
    /// Test construction with parameters
    RealQuadraticEquation equation(a, b, c);
    //
    /// Test for solutions
    BOOST_TEST(equation.first == 2.);
    BOOST_TEST(equation.second == 1.);
    BOOST_TEST(equation.solutions == 2);
  }

  /// Unit test for testing RealQuadraticEquation assignment
  BOOST_AUTO_TEST_CASE(RealQuadraticEquationAssignment)
  {
    double                a(1.0), b(-3.), c(2.);
    RealQuadraticEquation realQuadraticEquationObject(a, b, c);
    // operator == not implemented in this class
    //
    /// Test assignment (implicit)
    RealQuadraticEquation assignedRealQuadraticEquationObject(9., -3.5, 6.7);
    assignedRealQuadraticEquationObject = realQuadraticEquationObject;
    BOOST_TEST(assignedRealQuadraticEquationObject.first == 2.);
    BOOST_TEST(assignedRealQuadraticEquationObject.second == 1.);
    BOOST_TEST(assignedRealQuadraticEquationObject.solutions == 2);
    /// equality not written and not implicit
    // BOOST_TEST(assignedRealQuadraticEquationObject ==
    // realQuadraticEquationObject);
  }
  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts