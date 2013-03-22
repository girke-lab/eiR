#include "profiling.h"
#include <iostream>

Timer::Timer()
{
	reset();
}

Timer::~Timer()
{
}

void Timer::reset()
{
	s = UNINITIALIZED;
	accumulated = 0;
}
void Timer::start()
{
	if (s == RUNNING) 
		std::cout << "WARNING: timer restarted when it's running"
			<< std::endl;
	s = RUNNING;
	struct timeval now;
	gettimeofday(&now, NULL);
	last = now;
}

double Timer::pause()
{
	if (s != RUNNING) {
		std::cout << "ERROR: timer is not running when you try to pause."
			<< std::endl;
		return 0;
	}
	s = PAUSED;
	struct timeval now;
	gettimeofday(&now, NULL);
	double diff = now.tv_sec - last.tv_sec + 
		(now.tv_usec - last.tv_usec) * 0.000001;
	accumulated += diff;
	return accumulated;
}

double Timer::read()
{
	return accumulated;
}

/*
int main()
{
	Timer t;
	t.start();
	int k = 0;
	for (int i = 0; i < 1000000; i ++) {
		k += i * i;
		if (k > 10000) k -= i * i * i;
	}
	double tt = t.pause();
	std::cout << tt << std::endl;
	t.reset();
	t.start();
	k = 0;
	for (int i = 0; i < 1000000; i ++) {
		k += i * i;
		if (k > 10000) k -= i * i * i;
	}
	tt = t.pause();
	std::cout << tt << std::endl;
	k = 0;
	t.start();
	for (int i = 0; i < 1000000; i ++) {
		k += i * i;
		if (k > 10000) k -= i * i * i;
	}
	tt = t.pause();
	std::cout << tt << std::endl;

	return 0;
}
*/
