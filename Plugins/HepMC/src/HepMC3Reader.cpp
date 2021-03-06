// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/HepMC/HepMC3Reader.hpp"

bool
FW::HepMC3ReaderAscii::readEvent(HepMC::ReaderAscii&              reader,
                                 std::shared_ptr<HepMC::GenEvent> event)
{
  // Read event and store it
  return reader.read_event(*event);
}

bool
FW::HepMC3ReaderAscii::status(HepMC::ReaderAscii& reader)
{
  return !reader.failed();
}
