#include <sys/resource.h>
#include <sys/time.h>
#include <stdio.h>

static double realtime0;

double sys_cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double sys_realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return (tp.tv_sec + tp.tv_usec * 1e-6) - realtime0;
}

void sys_liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

void sys_init()
{
	sys_liftrlimit();
	realtime0 = sys_realtime();
}

const char *sys_timestamp()
{
	static char buf[256];
	double rt, ct;
	rt = sys_realtime();
	ct = sys_cputime();
	snprintf(buf, 255, "%.3f*%.2f", rt, ct/rt);
	return buf;
}
