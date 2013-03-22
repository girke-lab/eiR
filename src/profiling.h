#include <time.h>
#include <sys/time.h>

class Timer
{
	private:
		enum TimerStatus {
			UNINITIALIZED,
			RUNNING,
			PAUSED
		};
		struct timeval last;
		TimerStatus s;
		double accumulated;
	public:
		Timer();
		~Timer();
		void start();
		void reset();
		double read();
		double pause();
};
