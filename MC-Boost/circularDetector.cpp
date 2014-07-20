//
//  circularDetector.cpp
//  Xcode
//
//  Created by jacob on 7/20/11.
//  Copyright 2011 BMPI. All rights reserved.
//
#include "debug.h"
#include "circularDetector.h"
#include <iostream>
using std::cout;

#undef DEBUG

// Allows some tolerance that might arise from
// rounding errors when working with doubles.
static const double THRESHOLD = 1e-10;



CircularDetector::CircularDetector(const Aperture_Properties &props)
:Detector(),
 Aperture()
{
    detector_center.location.x = aperture_center.location.x = props.coordinates.x;
    detector_center.location.y = aperture_center.location.y = props.coordinates.y;
    detector_center.location.z = aperture_center.location.z = props.coordinates.z;
    
    radius = props.radius;
    
	/// Notify where the detector is added to the medium.
    cout << "-----------------------------------------------------\n"
         << "Adding a circular detector to the medium /\n"
         << "-----------------------------------------\n"
         << " Location: [x=" << aperture_center.location.x
         << ", y=" << aperture_center.location.y
         << ", z=" << aperture_center.location.z << "] (meters)\n";
    cout << " Radius: " << props.radius << " (meters)\n";
    
    // initialize the vector normal to the plane to have
    // direction since it is only used to know the direction
    // and not location that the plane faces.
    normalVector.withDirection();
    
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
		assert((xy_plane == true) || (xz_plane == true) || (yz_plane == true));  // One plane must be set.
	}
    
    
    
}


CircularDetector::CircularDetector(const double radius, const Vector3d &centerPoint)
:Detector(centerPoint)
{
    cout << " radius: " << radius << " (meters)\n";
    this->radius = radius;
}



CircularDetector::CircularDetector(const double radius, const boost::shared_ptr<Vector3d> centerPoint)
:Detector(centerPoint)
{
    cout << " radius: " << radius << " (meters)\n";
    this->radius = radius;
}



CircularDetector::~CircularDetector()
{
    
}


bool CircularDetector::photonHitDetector(const boost::shared_ptr<Vector3d> p0)
{

	boost::mutex::scoped_lock lock(m_detector_mutex);
    
    // Create a line segment from the center of the detector to the location of the photon.
    boost::shared_ptr<Vector3d> centerToPoint = aperture_center - (*p0);
    
    // If the dot product of the normal vector from the detector surface and the line segment
    // found above is zero (i.e. they are perpendicular) we know that the photon has hit
    // the detector and lies on the plane formed by the detector.
    if (fabs((VectorMath::dotProduct((*centerToPoint), normalVector)) - 0.0) <= THRESHOLD)
    {
        
        // Now calculate the distance from the center of the detector plane to the intersection point, and if
        // that distance is larger than the radius the line segment missed the detector.
        if (fabs(VectorMath::Distance(aperture_center, (*p0))) <= radius)
        {
#ifdef DEBUG
            cout << "intersection point with detector = " << p0;
#endif
            // Write the photon intersection coordinates to file for use with post-processing routines.
            //savePhotonExitCoordinates(p0);
            
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
    
}

// To calculate if a line segment made from the previous location of the photon to it's current
// location passed through the detector plane we calculate the point the line segment intersects
// the plane.  If this point lies on the plane and within the bounds of the detector, then we
// calculate the exit angle and save it with the coordinate location of exit.
//
// This algorithm uses properties of parametric lines and vectors to find the intersection point
// and the coordinates at which this happens.
// Parametric line: P = p0 + t(p1-p0)
//      Note: 
//           - P is a point (in 3D space) on the line formed from the point vectors P1-P0.
//           - t is between 0 and 1 (inclusive) and gives the point along the line where the
//             intersection with the plane occured.
//
// Having a point on the plane ('aperture_center', which is a Vector3d) and the normal unit vector ('normalVector')
// it is possible to project the line segment from aperture_center - p0 onto the normal vector to calculate 't', where
// 't' gives the point on the line segment that intersects the plane.  With 't' calculated it's possible
// to find the coordinates of the intersection.
bool CircularDetector::photonPassedThroughDetector(const boost::shared_ptr<Vector3d> p0,
                                                   const boost::shared_ptr<Vector3d> p1)
{
    
    
    boost::mutex::scoped_lock lock(m_detector_mutex);


    // Create the new vector that is the line segment from the previous position of the photon (p0) to 
    // it's current position (p1).
    boost::shared_ptr<Vector3d> lineSegment = (*p1) - (*p0);
    
    // Create the vector that is from the previous location of the photon (p0) to the point on the plane (aperture_center);
    boost::shared_ptr<Vector3d> pointToPlane = aperture_center - (*p0);
    
    // Check if the photon has moved passed the plane.  If not there is no way it could have crossed the plane,
    // so return false.  Otherwise continue the calculation.
    // This is achieved by calculating 't' from the parametric form of a line.  't' gives the point on the segment
    // that connectos points p1 and p0 of the photon's location that passes through the plane.  If t < 0 the line
    // segment lies completely on one half of the space seperated by the plane (i.e. it has already moved beyond the
    // plane and not able to intersect it.  If t > 0 then it lies completely on the other half of the plane and has 
    // not intersected it.  If 0 <= t <= 1 then the line segment has crossed the plane.
    double t = VectorMath::dotProduct(normalVector, *pointToPlane) / VectorMath::dotProduct(normalVector, *lineSegment);
    
    
    
    if (t >= 0.0f && t <= 1.0f)
    {
        // If we have made it in here we know the line segment passed through the plane.
        // Knowing 't' we can calculate the intersection point of the line segment with the plane.
        // This is essentially P(x,y,z) = p0(x,y,z) + (p1(x,y,z) - p0(x,y,z))*t;
        // FIXME:  The dereferencing of returned pointers is confusing due to the returned object (boost::shared_ptr<>)
        //          from the operator overloading.
        boost::shared_ptr<Vector3d> intersectionPoint = (*p0) + (*((*lineSegment)*t));
        
        
        // Now calculate the distance from the center of the detector plane to the intersection point, and if
        // that distance is larger than the radius the line segment missed the detector.
        if (VectorMath::Distance(aperture_center, (*intersectionPoint)) > radius)
        {
            return false;
        }
        else
        {
            cout << "intersection point = " << intersectionPoint;
            
            // Write the photon intersection coordinates to file for use with post-processing routines.
            savePhotonExitCoordinates(intersectionPoint);
            
            return true;
        }
    }
    else
    {
        return false;
    }
    
}



void CircularDetector::savePhotonExitCoordinates(const boost::shared_ptr<Vector3d> exitCoords)
{
    // Use logger class to write out data.
    //Logger::getInstance()->writeExitData(exitCoords);
	cout << "CircularDetector::savePhotonExitCoordinates() ... STUB!\n";
}


void CircularDetector::savePhotonExitData(const boost::shared_ptr<Vector3d> photonVector,
                                          const double weight,
                                          const bool tagged)
{
    //Logger::getInstance()->writeExitData(photonVector, weight, tagged);
	cout << "CircularDetector::savePhotonExitData() ... STUB!\n";
}



void CircularDetector::savePhotonExitWeight(void)
{
    cout << "CircularDetector::savePhotonExitWeight() ... STUB!\n";
}



