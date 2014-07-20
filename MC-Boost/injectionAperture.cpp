//
//  circularDetector.h
//  Xcode
//
//  Created by jacob on 7/20/11.
//  Copyright 2011 BMPI. All rights reserved.
//
#include <iostream>
using std::cout;

#include "injectionAperture.h"



InjectionAperture::InjectionAperture(void)
{
    aperture_center.location.x = 0.0f;
    aperture_center.location.y = 0.0f;
    aperture_center.location.z = 0.0f;
}


InjectionAperture::InjectionAperture(const Aperture_Properties &props)
{
    aperture_center.location.x = props.coordinates.x;
    aperture_center.location.y = props.coordinates.y;
    aperture_center.location.z = props.coordinates.z;
    
    radius = props.radius;
    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();
    
    /// Set the plane on which the detector lays.
	if (props.xy_plane)
	{
		setAperturePlaneXY();
	}
	else if (props.xz_plane)
	{
		setAperturePlaneXZ();
	}
	else if (props.yz_plane)
	{
		setAperturePlaneYZ();
	}
	else
	{
		cout << "!!!ERROR: Detector plane has not been defined.\n";
		assert((xy_plane == true) || (xz_plane == true) || (yz_plane == true));  // One plane must be set.
	}
}


InjectionAperture::~InjectionAperture(void)
{
    
}