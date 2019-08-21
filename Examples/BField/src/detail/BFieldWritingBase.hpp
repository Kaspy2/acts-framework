// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Plugins/Root/RootBFieldWriter.hpp"

namespace FW {
namespace BField {

  template <typename vmap_t, typename bfield_t>
  void
  writeField(vmap_t vm, std::shared_ptr<const bfield_t> bField)
  {
    // Write the interpolated magnetic field
    typename FW::Root::RootBFieldWriter<bfield_t>::Config writerConfig;
    if (vm["bf-out-rz"].template as<bool>())
      writerConfig.gridType = Root::GridType::rz;
    else
      writerConfig.gridType = Root::GridType::xyz;
    writerConfig.treeName = vm["bf-map-out"].template as<std::string>();
    writerConfig.fileName = vm["bf-file-out"].template as<std::string>();
    writerConfig.bField   = bField;
    std::cout << "setting rBounds" << std::endl;
    if (vm.count("bf-rRange") && vm.count("bf-zRange")) {
      auto rBounds = vm["bf-rRange"].template as<read_range>();
      auto zBounds = vm["bf-zRange"].template as<read_range>();
      writerConfig.rBounds
          = {{rBounds[0] * Acts::units::_mm, rBounds[1] * Acts::units::_mm}};
      writerConfig.zBounds
          = {{zBounds[0] * Acts::units::_mm, zBounds[1] * Acts::units::_mm}};
    }
    writerConfig.rBins   = vm["bf-rBins"].template as<size_t>();
    writerConfig.zBins   = vm["bf-ZBins"].template as<size_t>();
    writerConfig.phiBins = vm["bf-PhiBins"].template as<size_t>();

    FW::Root::RootBFieldWriter<bfield_t>::run(writerConfig);
  }
}  // namespace BField
}  // namespace FW
