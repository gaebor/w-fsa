#include "Isatty.h"

#ifdef _MSC_VER
#include <io.h>
bool IsTty(FILE* f)
{
    return _isatty(_fileno(f));
}
#else
#include <unistd.h>
bool IsTty(FILE* f)
{
    return isatty(fileno(f));
}

#endif // _MSCVER
