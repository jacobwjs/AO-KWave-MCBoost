//
//  aperture.cpp
//  Xcode
//
//  Created by Jacob Staley on 7/20/14.
//  Copyright 2011 BMPI, University of Twente. All rights reserved.
//

#include "aperture.h"

Aperture::Aperture(void)
{
    m_aperture_properties.center_coords.location.x = 0.0f;
    m_aperture_properties.center_coords.location.y = 0.0f;
    m_aperture_properties.center_coords.location.z = 0.0f;
    
    m_aperture_properties.xy_plane = m_aperture_properties.yz_plane = m_aperture_properties.xz_plane = false;
    
    m_aperture_properties.name = "empty";
}


Aperture::Aperture(const Aperture_Properties &props)
{
    m_aperture_properties.center_coords.location.x = props.center_coords.location.x;
    m_aperture_properties.center_coords.location.y = props.center_coords.location.y;
    m_aperture_properties.center_coords.location.z = props.center_coords.location.z;
    
    m_aperture_properties.xy_plane = props.xy_plane;
    m_aperture_properties.yz_plane = props.yz_plane;
    m_aperture_properties.xz_plane = props.xz_plane;
    
    m_aperture_properties.radius = props.radius;
    
    m_aperture_properties.name = props.name;
}

Aperture::~Aperture(void)
{
    
}