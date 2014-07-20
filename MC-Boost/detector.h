//
//  detector.h
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#ifndef DETECTOR_H
#define DETECTOR_H

#include <boost/thread/mutex.hpp>
#include "vector3D.h"
#include "logger.h"
#include "vectorMath.h"
using namespace VectorMath;
//#include <boost/math/complex/fabs.hpp>


class Detector
{
public:
    Detector(void);
    Detector(const double x, const double y, const double z);
    Detector(const Vector3d &centerPoint);
    Detector(const boost::shared_ptr<Vector3d> centerPoint);
    virtual ~Detector();
        
    virtual bool photonPassedThroughDetector(const boost::shared_ptr<Vector3d> p0,
                                             const boost::shared_ptr<Vector3d> p1) = 0;
    virtual bool photonHitDetector(const boost::shared_ptr<Vector3d> p0) = 0;
    virtual void savePhotonExitCoordinates(const boost::shared_ptr<Vector3d> exitCoords) = 0;
    virtual void savePhotonExitWeight(void) = 0;
    
    
    
protected:
    // Center coordinates of the detector in the medium. [cm]
    Vector3d detector_center;
    
	// Mutex to serialize access to the detector.
	boost::mutex m_detector_mutex;
};


#endif



