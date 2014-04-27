/**
 * Copyright 2010-2011 Sebastiaan J. van Schaik
 *
 * This file is part of PWAHStackTC.
 *
 * PWAHStackTC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PWAHStackTC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PWAHStackTC. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PERFORMANCETIMER_H_
#define PERFORMANCETIMER_H_
#include <sys/time.h>
#include <time.h>

class PerformanceTimer {
private:
	struct timeval* _startTime;
	double _storedRunTime;

public:
	PerformanceTimer();
	virtual ~PerformanceTimer();

	static PerformanceTimer start();
	static double diffTimeMilliSecs(const timeval* time1, const timeval* time2);
	static long diffTimeMicroSecs(const timeval* time1, const timeval* time2);
	double reset();
	double resetAndStop();
	double currRunTime();
	long currRunTimeMicro();
	void pause();
	void resume();
	bool running();
};

#endif /* PERFORMANCETIMER_H_ */
