//
//  logger.h
//  Xcode
//
//  Created by jacob on 7/18/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


// Logger singleton.
// NOTE:  Construction is NOT thread-safe, must be initialized in main before any threads are spawned.
#ifndef LOGGER_H
#define LOGGER_H

#include <vector>
#include <map>
#include <ctime>
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <boost/thread/mutex.hpp>
#include <boost/lexical_cast.hpp>

#include "RNG.h"
#include "multikey.h"
#include "loggerBase.h"



class Logger : public LoggerBase
{
public:
    Logger();
    ~Logger();
    
    virtual void Destroy();
    

private:
    Logger(const Logger&){};             // copy constructor is private
    
};

#endif  /// LOGGER_H
