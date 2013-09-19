//
//  sphere.cpp
//  Xcode
//
//  Created by jacob on 6/20/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "absorber.h"
#include "sphereAbsorber.h"
#include <cmath>
#include <iostream>
using std::cout;


SphereAbsorber::SphereAbsorber(const double radius, const double x, const double y, const double z)
:Absorber(x, y, z)
{
    this->radius = radius;
}

SphereAbsorber::SphereAbsorber(const double radius, const Vector3d &center)
:Absorber(center)
{
    this->radius = radius;
}

SphereAbsorber::SphereAbsorber(const double radius, const boost::shared_ptr<Vector3d> center)
:Absorber(center)
{
    this->radius = radius;
}


SphereAbsorber::~SphereAbsorber()
{
    // STUB
}


// FIXME: Need to check if the minimum distance from the line formed by the path
//        to the center of the absorber is within the radius bounds of the absorber.
//        If that is the case then the photon has moved through the absorber.
//        Based on the path length through the absorber a certain amount of 
//        energy needs to be dropped based on the local mu_a from the absorber.
bool SphereAbsorber::hitAbsorberBoundary(const boost::shared_ptr<Vector3d> photonLocation)
{
    // Convert the cartesian coordinates into a spherical coordinate structure.
    // sphereCoords scoords = cartesianToSpherical(center);
    
    // If the distance from the center of the absorber to the location of the photon is
    // larger than the radius of the absorber we know the photon hasn't crossed into the
    cout << "SphereAbsorber::hitAbsorberBoundary() Hit spherical absorber\n";
    
    return false;
    


}


// The idea is to use geometry and properties of vectors to find the shortest distance
// between a line extending from point A to point B and a central point C (i.e. the center
// of the absorber).  By taking the cross-product of A-B and A-C we have the area
// of the parallelogram between the vectors pointing from A to B and A to C.
// We know the area of a triangle is base*height/2.  By calculating the area of
// the parallelogram by taking the cross-product, as mentioned above, we have
// area = base*height*2 since the cross-product gives us twice the area of a triangle.
// Using this, we can calculate (A-B)X(A-C)/Length(AB) which
// leaves us the height, which we can compare to the radius of the absorber since
// this is the shortest distance to the line that MIGHT have passed through the absorber.
bool SphereAbsorber::crossedAbsorber(const boost::shared_ptr<Vector3d> A,
                             const boost::shared_ptr<Vector3d> B)
{
    // Subtract the previous location (B) of the photon from the current location (A)
    // to form a new vector.
    boost::shared_ptr<Vector3d> AB = (*A) - (*B);
    
    // Subtract the current location of the photon (A) from the center location of the
    // absorber (center) to yield a new vector.
    boost::shared_ptr<Vector3d> AC = (*A) - (*center);
    
    // Take the cross-product of AB and AC to get the area of the parallelogram formed
    // from A-B and A-C
    boost::shared_ptr<Vector3d> result = VectorMath::crossProduct(AB, AC);
    
    //
    cout << "SphereAbsorber::crossedAbsorber() stub\n";
    
    return false;
    
}


bool SphereAbsorber::inAbsorber(const boost::shared_ptr<Vector3d> photonLocation)
{
    if (inSphereVolume(photonLocation))
        return true;
    else
        return false;
}


bool SphereAbsorber::inSphereVolume(const boost::shared_ptr<Vector3d> photonLocation)
{
    // Calculate the distance from the photon to the center location
    // of the absorber.
    double dist_to_radius = VectorMath::Distance(photonLocation, center);
    
    // If the distance of the photon to the absorber's center is
    // larger than the absorber's radius, we return false to indicate
    // the photon is not within the absorber's bounds.  Otherwise return
    // true, since the photon has "hopped" to the absorber.
    if (dist_to_radius > radius)
        return false;
    else
        return true;
}


// FIXME: Verify this is correct.
sphereCoords SphereAbsorber::cartesianToSpherical(const coords &center)
{
    sphereCoords temp;
    temp.r      = radius;
    temp.theta  = acos(center.z / temp.r);
    temp.phi    = atan2(center.y, center.x);
}
