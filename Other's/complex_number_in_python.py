import cmath
import math
  
# Initializing real numbers
x = 1.0
y = 1.0
  
# converting x and y into complex number
z = complex(x, y);
  
# printing exponent of complex number
print ("The exponent of complex number is : ", end="")
print (cmath.exp(z))
  
# printing log form of complex number
print ("The log(base 10) of complex number is : ", end="")
print (cmath.log(z,10))