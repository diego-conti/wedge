/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unimib.it  *
 *                                                                         *
 *   This file is part of Wedge.                                           *
 *                                                                         *
 *   Wedge is free software; you can redistribute it and/or modify         *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   Wedge is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Wedge; if not, write to the                                *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "logging.h"
#include <time.h>
#include <cstdlib>

#define LOGGING_MAXIMUM_LENGTH 600

Wedge::LoggingLevel logging_level=Wedge::LOGGING_LEVEL_INFO;

namespace Wedge {
using namespace std;
using namespace GiNaC;

namespace internal { 
Log log;
LogBuffer log_buffer;

LogBuffer::LogBuffer() {
	maximum_length=LOGGING_MAXIMUM_LENGTH;
}

static terminate_handler terminate_function;

static void TerminateLoggingSession() {
	log<<"TERMINATED after " <<
		log.time_elapsed()  <<" milliseconds."<<
		 	endl;
	set_terminate(terminate_function);
	terminate();
}


Log::Log() : ofstream((string(getenv("HOME"))+"/wedge.log").c_str()) {
	start_clock=new clock_t(clock());
	time_t start_time;
	time(&start_time);
	(*this)<<"Logging session BEGINS at "<<
		ctime (&start_time)<<	 		
			endl;
	terminate_function=set_terminate(TerminateLoggingSession);
}		 
Log::~Log()
{		
	(*this)<<"Logging session ENDS after " << 
		time_elapsed()  <<" milliseconds."<<
		 	endl;	
	delete start_clock;
}


int Log::time_elapsed() 
{
	return (1000* (clock()-*start_clock))/CLOCKS_PER_SEC;
}

void Log::Message(const char* severity, const char* file, int line, const char* e, std::ostream& logBuffer)
{
	*this<<severity<<" at "<<time_elapsed()<< "ms in " << file << ", line "<<line<<", "<<e<<" = ";
	LogBuffer& f= static_cast<LogBuffer&>(logBuffer);
	*this<<f.text()<<endl;
}

string LogBuffer::text() {
	string result=str();
	if (result.size()>maximum_length)
		result="argument exceeding maximum length";
	str("");
	return result; 
}

}

std::ostream& operator<<(std::ostream& os, MaxLength x)
{
	internal::LogBuffer* log=dynamic_cast<internal::LogBuffer*>(&os);
	if (log!=NULL) log->maximum_length=x.maximum_length;
	else os<<"Error: manipulator MaxLength invoked on non-log stream"<<endl;
	return os;
}


}
