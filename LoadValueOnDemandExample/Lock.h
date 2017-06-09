#ifndef LOCK_H
#define LOCK_H

#include "ADDGLBL.h"

class Lock
{
public:
    Lock()
    {
        TecUtilLockStart(AddOnID);
    }
    ~Lock()
    {
        TecUtilLockFinish(AddOnID);
    }

private:
    Lock(const Lock&);
    Lock& operator=(const Lock&);
};

#endif // LOCK_H
