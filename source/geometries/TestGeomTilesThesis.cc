// ----------------------------------------------------------------------------
// petalosim | TestGeomTilesThesis.cc
//
// This class implements a test geometry to simulate the different tiles.
//
// The PETALO Collaboration
// ----------------------------------------------------------------------------

#include "TestGeomTilesThesis.h"
#include "TileHamamatsuVUV.h"
#include "TileHamamatsuBlue.h"
#include "TileFBK.h"
#include "PetMaterialsList.h"
#include "PetOpticalMaterialProperties.h"
#include "Na22Source.h"

#include "nexus/Visibilities.h"
#include "nexus/IonizationSD.h"
#include "nexus/FactoryBase.h"
#include "nexus/MaterialsList.h"
#include "nexus/OpticalMaterialProperties.h"
#include "nexus/SpherePointSampler.h"
#include "nexus/BoxPointSampler.h"

#include <G4GenericMessenger.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Material.hh>
#include <G4SDManager.hh>
#include <G4UserLimits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SubtractionSolid.hh>

using namespace nexus;

REGISTER_CLASS(TestGeomTilesThesis, GeometryBase)

TestGeomTilesThesis::TestGeomTilesThesis() : GeometryBase(),
                   visibility_(0),
                   reflectivity_(0),
                   tile_vis_(1),
                   tile_refl_(0.),
                   sipm_pde_(0.2),
                   source_pos_{},
                   tile_type_("HamamatsuVUV"),
                   n_tile_rows_(2),
                   n_tile_columns_(2),
                   wls_depth_(0.001 * mm),
                   max_step_size_(1. * mm)

{
  // Messenger
  msg_ = new G4GenericMessenger(this, "/Geometry/TestGeomTilesThesis/",
                                "Control commands of geometry TestGeomTilesThesis.");
  msg_->DeclareProperty("visibility", visibility_, "Visibility");
  msg_->DeclareProperty("surf_reflectivity", reflectivity_, "Reflectivity of the panels");
  msg_->DeclareProperty("tile_vis", tile_vis_, "Visibility of tiles");
  msg_->DeclareProperty("tile_refl", tile_refl_, "Reflectivity of SiPM boards");
  msg_->DeclareProperty("sipm_pde", sipm_pde_, "SiPM photodetection efficiency");

  msg_->DeclareProperty("n_tile_rows", n_tile_rows_, "Number of tiles per row");
  msg_->DeclareProperty("n_tile_columns", n_tile_columns_, "Number of tiles per column");

  msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  source_pos_,
                                "Set generation vertex.");

  msg_->DeclareProperty("tile_type", tile_type_,
                        "Type of the tile in the detection plane");

  G4GenericMessenger::Command &time_cmd =
      msg_->DeclareProperty("sipm_time_binning", time_binning_,
                            "Time binning for the sensor");
  time_cmd.SetUnitCategory("Time");
  time_cmd.SetParameterName("sipm_time_binning", false);
  time_cmd.SetRange("sipm_time_binning>0.");
}

TestGeomTilesThesis::~TestGeomTilesThesis()
{
}

void TestGeomTilesThesis::Construct()
{
  // LAB. Volume of air surrounding the detector ///////////////
  G4double lab_size = 1. * m;
  G4Box *lab_solid = new G4Box("LAB", lab_size/2., lab_size/2., lab_size/2.);

  G4Material *air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  lab_logic_ = new G4LogicalVolume(lab_solid, air, "LAB");
  lab_logic_->SetVisAttributes(G4VisAttributes::GetInvisible());
  this->SetLogicalVolume(lab_logic_);

  BuildSensors();

}


void TestGeomTilesThesis::BuildSensors()
{

  // TILES /////////////////////////////////////////////////////

  if (tile_type_ == "HamamatsuVUV") {
    tile_ = new TileHamamatsuVUV();
  } else if (tile_type_ == "HamamatsuBlue") {
    tile_ = new TileHamamatsuBlue();
  } else if (tile_type_ == "FBK") {
    tile_ = new TileFBK();
    tile_->SetPDE(sipm_pde_);
  } else {
    G4Exception("[TestGeomTilesThesis]", "BuildBox()", FatalException,
                "Unknown tile type for the detection plane!");
  }

  tile_->SetBoxGeom(1);
  // Construct the tile for the distances in the active
  tile_->SetTileVisibility(tile_vis_);
  tile_->SetTileReflectivity(tile_refl_);
  tile_->SetTimeBinning(time_binning_);

  tile_->Construct();
  tile_thickn_ = tile_->GetDimensions().z();



  // SiPMs /////////////////////////////////////////////////////

  G4double tile_size_x = tile_->GetDimensions().x();
  G4double tile_size_y = tile_->GetDimensions().y();
  full_row_size_ = n_tile_columns_ * tile_size_x;
  full_col_size_ = n_tile_rows_ * tile_size_y;

  G4LogicalVolume *tile_logic = tile_->GetLogicalVolume();

  G4String vol_name;
  G4int copy_no = 0;

  G4double z_pos = -10 *cm;
  for (G4int j = 0; j < n_tile_rows_; j++)
  {
    G4double y_pos = full_col_size_/2. - tile_size_y/2. - j*tile_size_y;
    for (G4int i = 0; i < n_tile_columns_; i++)
    {
      G4double x_pos = -full_row_size_/2. + tile_size_x/2. + i*tile_size_x;
      vol_name = "TILE_" + std::to_string(copy_no);

      new G4PVPlacement(0, G4ThreeVector(x_pos, y_pos, z_pos), tile_logic,
                        vol_name, lab_logic_, false, copy_no, false);
      copy_no += 1;
    }
  }

  G4RotationMatrix rot;
  rot.rotateY(pi);

}

G4ThreeVector TestGeomTilesThesis::GenerateVertex(const G4String &region) const
{
  G4ThreeVector vertex(0., 0., 0.);

  if (region == "CENTER")
    {
      return vertex;
    }
  else if (region == "AD_HOC")
    {
      vertex = source_pos_;
    }
  else if (region == "SOURCE")
    {
      vertex = source_gen_->GenerateVertex("VOLUME");
    }
  else
    {
      G4Exception("[TestGeomTilesThesis]", "GenerateVertex()", FatalException,
                  "Unknown vertex generation region!");
    }
  return vertex;
}
