// ----------------------------------------------------------------------------
// petalosim | TestGeomTilesThesis.h
//
// This class implements the geometry of a symmetric box of LXe.
//
// The PETALO Collaboration
// ----------------------------------------------------------------------------

#ifndef PET_BOX_H
#define PET_BOX_H

#include "nexus/GeometryBase.h"
#include "TileGeometryBase.h"

class G4GenericMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

class TileHamamatsuVUV;
class TileHamamatsuBlue;
class TileFBK;

namespace nexus
{
class SpherePointSampler;

class BoxPointSampler;
}

using namespace nexus;

class TestGeomTilesThesis : public GeometryBase
{

public:
  // Constructor
  TestGeomTilesThesis();
  //Destructor
  ~TestGeomTilesThesis();

  /// Generate a vertex within a given region of the geometry
  G4ThreeVector GenerateVertex(const G4String &region) const;

  TileGeometryBase *tile_;

private:
  void Construct();
  void BuildBox();
  void BuildSensors();

  G4LogicalVolume *lab_logic_;

  G4bool visibility_;
  G4double reflectivity_;
  G4bool tile_vis_;
  G4double tile_refl_;
  G4double sipm_pde_;

  G4ThreeVector source_pos_;

  G4String tile_type_;
  G4double time_binning_;

  G4double n_tile_rows_, n_tile_columns_;
  G4double tile_thickn_, full_row_size_, full_col_size_;

  G4double wls_depth_;

  G4double max_step_size_;

  /// Messenger for the definition of control commands
  G4GenericMessenger* msg_;

  SpherePointSampler* source_gen_;
};

#endif
