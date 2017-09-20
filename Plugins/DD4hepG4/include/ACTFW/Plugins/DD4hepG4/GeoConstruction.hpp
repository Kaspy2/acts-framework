#ifndef DD4HEPG4_GEOCONSTRUCTION_H
#define DD4HEPG4_GEOCONSTRUCTION_H

// DD4hep
#include "DDG4/Geant4GeometryInfo.h"

// Geant4
#include "G4VUserDetectorConstruction.hh"

namespace dd4hep {
        class Detector;
}
/// Temporary borrowed from FCCSW -> will be replaced later
/** @class GeoConstruction DetectorDescription/DetDesServices/src/GeoConstruction.h GeoConstruction.h
 *
 *  Class to create Geant4 detector geometry from TGeo representation
 *  On demand (ie. when calling "Construct") the DD4hep geometry is converted
 *  to Geant4 with all volumes, assemblies, shapes, materials etc.
 *
 *  @author Markus Frank
 *  @author Anna Zaborowska
 */

namespace FW {

namespace DD4hepG4 {

    class GeoConstruction : public G4VUserDetectorConstruction {
    public:
        /// Constructor
        GeoConstruction(dd4hep::Detector& lcdd);
        /// Default destructor
        virtual ~GeoConstruction() = default;
        /// Geometry construction callback: Invoke the conversion to Geant4
        /// All volumes (including world) are deleted in ~G4PhysicalVolumeStore()
        virtual G4VPhysicalVolume* Construct() final;
    private:
        /// Reference to geometry object
        dd4hep::Detector& m_lcdd;
    };
} // namespace DD4hepG4
} // namespace FW

#endif /* DETDESSERVICES_GEOCONSTRUCTION_H */
