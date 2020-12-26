//sign function

#include "functions.h"

double sign(double a, double b)
{
a=abs(a);

if(b!=0) a=a*b/abs(b);

return(a);
}
