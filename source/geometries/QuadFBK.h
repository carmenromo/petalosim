#ifndef QUAD_FBK_H
#define QUAD_FBK_H

#include "BaseGeometry.h"
#include <G4ThreeVector.hh>

class G4GenericMessenger;
namespace nexus {class SiPMpetFBK;}

namespace nexus {

  class QuadFBK: public BaseGeometry
  {
  public:
    /// Constructor
    QuadFBK();
    /// Destructor
    ~QuadFBK();
    
    /// Return dimensions of the SiPM
    //G4ThreeVector GetDimensions() const;
    
    /// Invoke this method to build the volumes of the geometry
    void Construct();
    
  private:
    //G4ThreeVector _dimensions; ///< external dimensions of the SiPMpet

    // Visibility of the tracking plane
    G4bool visibility_;

     // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    SiPMpetFBK* sipm_;

  };


} // end namespace nexus

#endif
