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
    aperture_center.location.x = 0.0f;
    aperture_center.location.y = 0.0f;
    aperture_center.location.z = 0.0f;
    
    xz_plane = yz_plane = xz_plane = false;
}

Aperture::~Aperture(void)
{
    
}