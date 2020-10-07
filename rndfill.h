#ifndef RNDFILL_H__
#define RNDFILL_H__

#ifdef __GNUC__
#include <sched.h>
#include <unistd.h> 
#endif

#ifdef _MSC_VER
#define DUAL_CORE

#if defined( DUAL_CORE )
#  define WINDOWS_LEAN_AND_MEAN
#  include <windows.h>
#endif
#include <string.h>
#include <math.h>

#include <intrin.h>
#pragma intrinsic( __rdtsc )

__inline unsigned long long read_tsc(void)
{
	return __rdtsc();
}

#else
#ifdef __GNUC__
#include <sys/resource.h>
#include <x86intrin.h>
inline unsigned long long read_tsc(void)
{
#if defined(__i386__)
	unsigned long long cycles;
	__asm__ volatile (".byte 0x0f, 0x31" : "=A"(cycles));
	return cycles;
#else
#if defined(__x86_64__)
	unsigned int hi, lo;
	__asm__ volatile ("rdtsc" : "=a"(lo), "=d"(hi));
	return (((unsigned long long)lo) | ((unsigned long long)(hi)<<32));
#else
#error "Unsupported architecture for counting cycles"
#endif
#endif
}
#endif
#endif

#define RAND(a,b) (((a = 36969 * (a & 65535) + (a >> 16)) << 16) + \
	(b = 18000 * (b & 65535) + (b >> 16))  )

void block_rndfill(unsigned char *buf, const int len)
{
	static unsigned int a[2], mt = 1, count = 4;
	static unsigned char r[4];
	int                  i;

	if(mt) { mt = 0; *(unsigned long long*)a = read_tsc(); }

	for(i = 0; i < len; ++i)
	{
		if(count == 4)
		{
			*(unsigned int*)r = RAND(a[0], a[1]);
			count = 0;
		}

		buf[i] = r[count++];
	}
}
#endif
