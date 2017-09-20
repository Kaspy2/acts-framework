///////////////////////////////////////////////////////////////////
// GeometryService.hpp
///////////////////////////////////////////////////////////////////

#ifndef DD4HEP_GEOMETRYSERVICE_H
#define DD4HEP_GEOMETRYSERVICE_H 1

#include <memory>
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/GeometryInterfaces/IDD4hepService.hpp"
#include "ACTFW/GeometryInterfaces/ITGeoService.hpp"
#include "ACTFW/GeometryInterfaces/ITrackingGeometryService.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "TGeoNode.h"

namespace FW {

namespace DD4hep {

/// @class GeometryService
///
/// @brief service creating geometries from dd4hep input
///
/// The GeometryService creates the DD4hep, the TGeo and the ACTS
/// TrackingGeometry
/// from DD4hep xml input. The geometries are created only on demand.

class GeometryService : public FW::IDD4hepService,
                        public FW::ITGeoService,
                        public FW::ITrackingGeometryService
{
public:
  /// @class Config
  /// nested config file of the GeometryService
  class Config
  {
  public:
    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;
    /// XML-file with the detector description
    std::string xmlFileName;
    /// Logger for the geometry transformation
    Acts::Logging::Level lvl;
    /// The name of the service
    std::string name;
    /// Binningtype in phi
    Acts::BinningType bTypePhi;
    /// Binningtype in r
    Acts::BinningType bTypeR;
    /// Binningtype in z
    Acts::BinningType bTypeZ;
    /// The tolerance added to the geometrical extension in r
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing the
    /// layers (e.g. barrel, endcap volumes) have no specific shape (assemblies)
    double envelopeR;
    /// The tolerance added to the geometrical extension in z
    /// of the layers contained to build the volume envelope around
    /// @note this parameter only needs to be set if the volumes containing the
    /// layers (e.g. barrel, endcap volumes) have no specific shape (assemblies)
    double envelopeZ;
    ///  Flag indicating if the
    /// Acts::DigitizationModule (needed for Acts geometric digitization) will
    /// be
    /// build for every single sensitive DD4hep DetElement translating directly
    /// the
    /// DD4hep Segmentation.
    /// @attention Turning on this flag can be very time and memory consuming!
    /// If
    /// different modules are sharing the same segmentation (which will be the
    /// case
    /// most of the times) please use the
    /// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
    /// constructor.
    /// More information on the usage can be found
    /// in the description of the Acts::ActsExtension class.
    bool buildDigitizationModules;

    Config(const std::string&   lname = "GeometryService",
           Acts::Logging::Level level = Acts::Logging::INFO)
      : logger(Acts::getDefaultLogger(lname, level))
      , xmlFileName("")
      , lvl(level)
      , name(lname)
      , bTypePhi(Acts::equidistant)
      , bTypeR(Acts::equidistant)
      , bTypeZ(Acts::equidistant)
      , envelopeR(0.)
      , envelopeZ(0.)
      , buildDigitizationModules(false)
    {
    }
  };
  /// Constructor
  GeometryService(const Config& cfg);

  /// Virtual destructor
  ~GeometryService() override;

  /// Framework name() method
  std::string
  name() const final override;

  /// Interface method to access the DD4hep geometry
  /// @return The world DD4hep DetElement
  dd4hep::DetElement
  dd4hepGeometry() final override;

  /// Interface method to Access the TGeo geometry
  /// @return The world TGeoNode (physical volume)
  TGeoNode*
  tgeoGeometry() final override;

  /// Interface method to access to the interface of the DD4hep geometry
  dd4hep::Detector*
  lcdd() final override;

  /// Interface method to access the ACTS TrackingGeometry
  std::unique_ptr<const Acts::TrackingGeometry>
  trackingGeometry() final override;

private:
  /// Private method to initiate building of the DD4hep geometry
  FW::ProcessCode
  buildDD4hepGeometry();

  /// Private method to initiate building of the ACTS tracking geometry
  FW::ProcessCode
  buildTrackingGeometry();

  /// The config class
  Config m_cfg;
  /// Pointer to the interface to the DD4hep geometry
  dd4hep::Detector* m_lcdd;
  /// The world DD4hep DetElement
  dd4hep::DetElement m_dd4hepGeometry;
  /// The ACTS TrackingGeometry
  std::unique_ptr<const Acts::TrackingGeometry> m_trackingGeometry;

  /// Private access to the logging instance
  const Acts::Logger&
  logger() const
  {
    return *m_cfg.logger;
  }
};

} // namespace DD4hep
} // namespace FW

#endif  // DD4HEP_GEOMETRYSERVICE_H
