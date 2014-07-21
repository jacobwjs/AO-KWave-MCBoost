//
//  circularDetector.h
//  Xcode
//
//  Created by jacob on 7/20/11.
//  Copyright 2011 BMPI. All rights reserved.
//

#ifndef INJECTIONAPERTURE_H_
#define INJECTIONAPERTURE_H_

#include "aperture.h"

class InjectionAperture : public Aperture 
{
public:

    InjectionAperture(void);
    InjectionAperture(const Aperture_Properties &props);
    ~InjectionAperture(void);

};

#endif   /// INJECTIONAPERTURE_H_


