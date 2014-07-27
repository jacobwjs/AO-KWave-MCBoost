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

#include <string>
using std::string;


typedef struct {
    double radius;

    Vector3d center_coords;
    Vector3d normalVector;
    
	bool xy_plane;
	bool xz_plane;
	bool yz_plane;
    
    std::string name;
    
} Aperture_Properties;




class Aperture
{
public:
    
    Aperture();
    Aperture(const Aperture_Properties &props);
    virtual ~Aperture();
    
    
    /// Return the radius of the aperture. [meters]
    double Get_radius()     const {return m_aperture_properties.radius;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_x_coord()    const {return m_aperture_properties.center_coords.location.x;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_y_coord()    const {return m_aperture_properties.center_coords.location.y;};
    
    /// Return the x-coordinate for the center of the aperture.
    double Get_z_coord()    const {return m_aperture_properties.center_coords.location.z;};
    
    
    /// Is the aperture on the x-y plane
    bool    Is_XY_plane()   const {return m_aperture_properties.xy_plane;};
    
    /// Is the aperture on the y-z plane
    bool    Is_YZ_plane()   const {return m_aperture_properties.yz_plane;};
    
    /// Is the aperture on the x-z plane
    bool    Is_XZ_plane()   const {return m_aperture_properties.xz_plane;};
    
    virtual std::string Get_name() const {return m_aperture_properties.name;};
    
    
    /// Set the plan on which the aperture resides.
    virtual void setAperturePlaneXY(void)
    {
        // Set which plane the detector resides.
        m_aperture_properties.xz_plane = false;
        m_aperture_properties.yz_plane = false;
        m_aperture_properties.xy_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        m_aperture_properties.normalVector.setDirX(0.0f);
        m_aperture_properties.normalVector.setDirY(0.0f);
        m_aperture_properties.normalVector.setDirZ(1.0f); m_aperture_properties.normalVector.location.z = 1.0f;
        
    }
    
    /// Set the plan on which the aperture resides.
    virtual void setAperturePlaneXZ(void)
    {
        // Set which plane the detector resides.
        m_aperture_properties.yz_plane = false;
        m_aperture_properties.xy_plane = false;
        m_aperture_properties.xz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        m_aperture_properties.normalVector.setDirX(0.0f);
        m_aperture_properties.normalVector.setDirY(1.0f); m_aperture_properties.normalVector.location.y = 1.0f;
        m_aperture_properties.normalVector.setDirZ(0.0f);
    }
    
    /// Set the plan on which the aperture resides.    
    virtual void setAperturePlaneYZ(void)
    {
        // Set which plane the detector resides.
        m_aperture_properties.xz_plane = false;
        m_aperture_properties.xy_plane = false;
        m_aperture_properties.yz_plane = true;
        
        // Set the direction that the vector that is normal to the plane.
        m_aperture_properties.normalVector.setDirX(1.0f); m_aperture_properties.normalVector.location.x = 1.0f;
        m_aperture_properties.normalVector.setDirY(0.0f);
        m_aperture_properties.normalVector.setDirZ(0.0f);
    }
    
    
    
    
protected:
    // Center coordinates of the detector in the medium. [meters]
    //Vector3d aperture_center;
    
    // Vector that is normal to the plane.
    //Vector3d normalVector;
    

    /// The properties of the aperture.
    Aperture_Properties m_aperture_properties;
    
};

#endif  /// APERTURE_H_

