#include <time.h>

#if !defined ADDON
#define ADDON
#endif /* ADDON */
#include "TECADDON.h"

#include "pctdone.h"

#define CLOCK_INTERVAL_SECONDS (0.1)

static int     Count;
static int     ClockInterval;
static clock_t PreviousClock;
static clock_t CurrentClock;

extern void ResetPercentDone(void)
{
    Count = 0;
    ClockInterval = 1;
    CurrentClock = clock();
}

extern Boolean_t CheckPercentDone(int PercentDone)
{
    Boolean_t IsOk = TRUE;

    Count++;
    if (Count % ClockInterval == 0)
    {
        double ElapsedTime;

        PreviousClock = CurrentClock;
        CurrentClock  = clock();
        ElapsedTime = (double)(CurrentClock - PreviousClock) / CLOCKS_PER_SEC;
        if (ElapsedTime < CLOCK_INTERVAL_SECONDS)
        {
            ClockInterval *= 2;
        }
        else
        {
            if (ClockInterval > 1 && ElapsedTime > 2.0 * CLOCK_INTERVAL_SECONDS)
            {
                ClockInterval = ClockInterval * 2 / 3;
            }

            if (!TecUtilDialogCheckPercentDone(PercentDone))
            {
                IsOk = FALSE;
            }
        }
    }

    return(IsOk);
}

