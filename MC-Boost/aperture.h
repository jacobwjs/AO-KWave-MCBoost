//
//  aperture.h
//  Xcode
//
//  Created by Jacob Staley on 7/20/14.
//  Copyright 2011 BMPI, University of Twente. All rights reserved.
//

#ifndef APERTURE_H_
#define APERTURE_H_

#include "vector3D.h"
#include "coordinates.h"


typedef struct {
    double radius;

    coords coordinates;
    
	bool xy_plane;
	bool xz_plane;
	bool yz_plane;
    
} Aperture_Properties;




class Aperture
{
public:
    
    Aperture(void);
    ~Aperture(void);
    
    
    /// Return the radius of the aperture. [meters]
    double Get_radius()     const {return radius;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_x_coord()    const {return aperture_center.location.x;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_y_coord()    const {return aperture_center.location.y;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_z_coord()    const {return aperture_center.location.z;};
    
    
    /// Is the aperture on the x-y plane
    bool    Is_XY_plane()   const {return xy_plane;};
    
    /// Is the aperture on the y-z plane
    bool    Is_YZ_plane()   const {return yz_plane;};
    
    /// Is the aperture on the x-z plane
    bool    Is_XZ_plane()   const {return xz_plane;};
    
    
    /// Set the plan on which the aperture resides.
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
    
    /// Set the plan on which the aperture resides.
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
    
    /// Set the plan on which the aperture resides.    
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
    
    // Radius of the aperture. [meters]
    double radius;
    
    // Vector that is normal to the plane.
    Vector3d normalVector;
    
    // possible planes that the detector can be placed in 3D space.
    bool xy_plane;
    bool xz_plane;
    bool yz_plane;
    
};

#endif  /// APERTURE_H_

