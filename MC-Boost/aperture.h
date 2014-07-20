//
//  aperture.h
//  Xcode
//
//  Created by Jacob Staley on 7/20/14.
//  Copyright 2011 BMPI, University of Twente. All rights reserved.
//

#include "vector3D.h"


typedef struct {
    double radius;
    double x_coord;
    double y_coord;
    double z_coord;
    
	bool xy_plane;
	bool xz_plane;
	bool yz_plane;
    
} Aperture_Properties;




class Aperture
{
public:
    
    Aperture(void);
    ~Aperture(void);
    
    
    virtual void setAperturePlaneXY(void)
    {
        // Set which plane the detector resides.
        xz_plane = false;
        yz_plane = false;
        xy_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(0.0f);
        normalVector.setDirY(0.0f);
        normalVector.setDirZ(1.0f); normalVector.location.z = 1.0f;
        
    }
    
    
    virtual void setAperturePlaneXZ(void)
    {
        // Set which plane the detector resides.
        yz_plane = false;
        xy_plane = false;
        xz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(0.0f);
        normalVector.setDirY(1.0f); normalVector.location.y = 1.0f;
        normalVector.setDirZ(0.0f);
    }
    
    
    virtual void setAperturePlaneYZ(void)
    {
        // Set which plane the detector resides.
        xz_plane = false;
        xy_plane = false;
        yz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        normalVector.setDirX(1.0f); normalVector.location.x = 1.0f;
        normalVector.setDirY(0.0f);
        normalVector.setDirZ(0.0f);
    }
    
    
    
    
protected:
    // Center coordinates of the detector in the medium. [meters]
    Vector3d aperture_center;
    
    // Vector that is normal to the plane.
    Vector3d normalVector;
    
    // possible planes that the detector can be placed in 3D space.
    bool xy_plane;
    bool xz_plane;
    bool yz_plane;
    
};

