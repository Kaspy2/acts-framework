/// @file
/// @date 2016-05-23
/// @author Andreas Salburger
/// @author Moritz Kiehnn <msmk@cern.ch>

#ifndef ACTFW_ISERVICE_H
#define ACTFW_ISERVICE_H

#include <string>

#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Interface for common services.
/// @todo Remove once all initializers and finalizers are gone
class IService
{
public:
  /// Virtual Destructor
  virtual ~IService() = default;

  /// Framework name() method
  virtual std::string
  name() const = 0;

  /// Framework intialize method
  /// @todo Remove once all custom implementations are gone
  virtual ProcessCode
  initialize() { return ProcessCode::SUCCESS; }

  /// Framework finalize method
  /// @todo Remove once all custom implementations are gone
  virtual ProcessCode
  finalize() { return ProcessCode::SUCCESS; }
};

}  // namespace FW

#endif  // ACTFW_ISERVICE_H
