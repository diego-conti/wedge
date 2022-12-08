/***************************************************************************
 *   Copyright (C) 2007, 2008, 2009 by Diego Conti, diego.conti@unipi.it  *
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
#ifndef LOGGING_H_
#define LOGGING_H_
/** @ingroup Base */ 

/** @{ 
 * @file logging.h
 * @brief Macros for selective logging
 * 
 * This file defines several macros to send messages to a log file. A severity level is associated to each message;
 * a log message is ignored if its severity is below the value of the variable loggingLevel. 
 * 
 * @remark When debugging, most debuggers will enable one to change the value of loggingLevel at execution time.
 * This can be useful to tailor the behaviour of the logging mechanism according to one's need.  
 * 
 * @remark The argument of each LOG macro is output to a stream; therefore, it can have the form, say, latex<<expression.
 */
 

#define WEDGE_LOGGING_ENABLED
#ifdef WEDGE_LOGGING_ENABLED
#define LOG_DEBUG(e) if (logging_level<=LOGGING_LEVEL_DEBUG) Wedge::internal::log.Message("DEBUG", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_INFO(e) if (logging_level<=LOGGING_LEVEL_INFO) Wedge::internal::log.Message("INFO", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_WARN(e) if (logging_level<=LOGGING_LEVEL_WARN) Wedge::internal::log.Message("WARN", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_ERROR(e) if (logging_level<=LOGGING_LEVEL_ERROR) Wedge::internal::log.Message("ERROR", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_CRIT(e) if (logging_level<=LOGGING_LEVEL_CRIT) Wedge::internal::log.Message("CRIT", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_ALERT(e) if (logging_level<=LOGGING_LEVEL_ALERT) Wedge::internal::log.Message("ALERT", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_FATAL(e) if (logging_level<=LOGGING_LEVEL_FATAL) Wedge::internal::log.Message("FATAL", __FILE__,__LINE__,#e, Wedge::internal::log_buffer<<e);
#define LOG_MSG(e) if (logging_level<=LOGGING_LEVEL_WARN) Wedge::internal::log.Message("MSG", __FILE__,__LINE__, "message", Wedge::internal::log_buffer<<e);

#else
#define LOG_DEBUG(e) 
#define LOG_INFO(e) 
#define LOG_NOTICE(e)
#define LOG_WARN(e) 
#define LOG_ERROR(e)
#define LOG_CRIT(e) 
#define LOG_ALERT(e)
#define LOG_FATAL(e)
#define LOG_MSG(e)
#endif

namespace GiNaC {
	std::ostream & latex(std::ostream & os);
}

namespace Wedge {
/** @brief The eight possible levels of severity associated to a log message
 */
enum LoggingLevel {	
	LOGGING_LEVEL_DEBUG=1,
	LOGGING_LEVEL_INFO =2,
	LOGGING_LEVEL_NOTICE =3,
	LOGGING_LEVEL_WARN =4,
	LOGGING_LEVEL_ERROR =5,
	LOGGING_LEVEL_CRIT =6,
	LOGGING_LEVEL_ALERT =7,
	LOGGING_LEVEL_FATAL =8
};

/** @brief Base class for selective logging
 *
 * Classes derived from Logging have an associated logging level, which is the minimum severity level of messages
 * that are actually output to the log, relative to LOG macros invoked inside the class.  
 */
  
struct Logging {	
	LoggingLevel logging_level; ///< The minimum severity level for log messages that are actually output to the log  
	
	/** @brief Construct the object, initializing the logging level.
	 */
	Logging(LoggingLevel loggingLevel=LOGGING_LEVEL_INFO) {this->logging_level=loggingLevel;}
};

/** @brief Manipulator to change the log message maximum length 
 * 
 * Should only be used inside the LOG_* macros
 */
struct MaxLength {
	int maximum_length;
	MaxLength(int m) {maximum_length=m;}
};
std::ostream& operator<<(std::ostream& os, MaxLength x);

namespace internal {
class Log : public std::ofstream {
	clock_t* start_clock;
public:
	Log();
	~Log();
	void Message(const char* severity, const char* file, int line, const char* e, std::ostream& logBuffer);
	int time_elapsed();
};

class LogBuffer : public std::stringstream {
public:
	LogBuffer();
	int maximum_length;
	std::string text();
};

extern Log log;
extern LogBuffer log_buffer;

}
} 

/** @brief Global logging level
 * 
 * This is the minimum severity level of messages that are actually output to the log, relative to LOG
 * macros invoked either in the global namespace or in classes that do not derive from Logging.
 */
extern Wedge::LoggingLevel logging_level;  

 /** @} */
#endif /*LOGGING_H_*/
