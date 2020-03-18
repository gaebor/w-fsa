#include <Isatty.h>

#include <stdio.h>

#ifdef _MSC_VER
#include <io.h>

bool IsStderrTty()
{
    return _isatty(_fileno(stderr));
}
#else
#include <unistd.h>
bool IsStderrTty()
{
    return isatty(fileno(stderr));
}

#endif // _MSCVER
