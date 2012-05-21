// Created: olav, Mon Nov 17 16:56:02 CET 2008

#include "olav.h"

void status(unsigned long x, unsigned long R, unsigned long N)
{
	// this assumes that x will take any value between 0 and R, so no jumping ahead
  if(R>N+1) R=N;
  if(x==0) // init
  {
  	//cout <<"[";
  	for(unsigned long i=0; i<R; i++)
  	  std::cout <<".";
  	//cout <<"]\b";
  	for(unsigned long i=0; i<R; i++)
  	  std::cout <<"\b";
  }

  // if((R>0)&&(x%((long)floor((double)N/R))==0))
	// Neuer Code, sollte jetzt auch funktionieren, wenn R kein Teiler von N ist.
  if((R>0) && ((int)floor((((double)x)*R)/N) > (int)floor((((double)x)-1.)*R/N)))
    std::cout <<"#"<<std::flush;
}

std::string sec2string(double seconds) {
  if(seconds > double(std::numeric_limits<long>::max())) return "inf";
  if(seconds < double(std::numeric_limits<long>::min())) return "-inf";
  return sec2string((long)seconds);
}
std::string sec2string(long seconds)
{
	std::ostringstream text(std::ostringstream::out);
	
	if(seconds<0) {
		text <<"-";
		seconds *= -1;
	}

	if(seconds>3600*24) {
		text << seconds/(3600*24) << "d ";
		seconds = seconds % (3600*24);
	}
	
	if(seconds>3600) {
		text << seconds/3600 << "h ";
		seconds = seconds % 3600;
	}
	
	if(seconds>60) {
		text << seconds/60 << "m ";
		seconds = seconds % 60;
	}
	text << seconds << "s";

	return text.str();
}


std::string ETAstring(long completed, long total, double elapsedtime)
{
	std::string result = "unknown";

	if (completed < total) {
		if ((completed > 0)&&(elapsedtime > .75))
			result = sec2string(double(total-completed)*elapsedtime/double(completed));
	}
	else result = sec2string((long)0);

	return result;
}


std::string ETAstring(long completed, long total, time_t start)
{
  time_t now;
	time(&now);

	return ETAstring(completed,total,difftime(now,start));
}
