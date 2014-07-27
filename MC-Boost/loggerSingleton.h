//
//  loggerSingleton.h
//  Xcode
//
//  Created by Jacob Staleu on 7/18/11.
//  Copyright 2014 BMPI, University of Twente. All rights reserved.
//

// Logger singleton.
// NOTE:  Construction is NOT thread-safe, must be initialized in main before any threads are spawned.
#ifndef LOGGERSINGLETON_H
#define LOGGERSINGLETON_H

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







class LoggerSingleton : public LoggerBase
{
public:
    void Destroy();
    
    static LoggerSingleton * getInstance();
    
	
    
private:
    LoggerSingleton() {};                            // default constructor is private
    LoggerSingleton(const LoggerSingleton&){};             // copy constructor is private
    ~LoggerSingleton();
    
    LoggerSingleton& operator=(const LoggerSingleton& rhs) {};  // assignment operator is private
    
    static LoggerSingleton * pInstance;
    
};

#endif  /// LOGGERSINGLETON_H
