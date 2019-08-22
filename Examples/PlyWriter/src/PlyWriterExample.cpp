// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "detail/PlyWriterExampleBase.hpp"

/// @brief main executable
///
/// @param argc The argument count
/// @param argv The argument list
int
main(int argc, char* argv[])
{
  // --------------------------------------------------------------------------------
  GenericOptions  genericOptions;
  GenericGeometry genericGeometry;

  // now process it
  return propagationExample(argc, argv, genericOptions, genericGeometry);
}

/*

example run command:

./ACTFWPlyWriterExample                   (binary)
-n1                                       (number of runs)
--mat-input-type file                     (material description input type)
--mat-input-file mbp-40-1000.json         (mapped-binned-proto file)
-j1                                       (number of jobs)
--ply-output-file output.ply              (output file name)
--prop-eta-range -6 6                     (eta range for propagator)
--ply-colour-range 255 0 0 0 0 255        (global min and max colours - 2 RGB
values)

*/