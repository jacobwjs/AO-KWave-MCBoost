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
    m_aperture_properties.center_coords.location.x = 0.0f;
    m_aperture_properties.center_coords.location.y = 0.0f;
    m_aperture_properties.center_coords.location.z = 0.0f;
}


InjectionAperture::InjectionAperture(const Aperture_Properties &props)
: Aperture(props)
{
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    m_aperture_properties.normalVector.withDirection();
    
	/// Set the plane on which the detector lays.
	if (props.xy_plane)
	{
		setAperturePlaneXY();
        cout << " plane: x-y\n";
	}
	else if (props.xz_plane)
	{
		setAperturePlaneXZ();
        cout << " plane: x-z\n";
	}
	else if (props.yz_plane)
	{
		setAperturePlaneYZ();
        cout << " plane: y-z\n";
	}
	else
	{
		cout << "!!!ERROR: Detector plane has not been defined.\n";
		assert((m_aperture_properties.xy_plane == true) ||
               (m_aperture_properties.xz_plane == true) ||
               (m_aperture_properties.yz_plane == true));  // One plane must be set.
	}
}


InjectionAperture::~InjectionAperture(void)
{
    
}