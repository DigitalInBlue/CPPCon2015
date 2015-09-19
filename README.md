#Demystifying Floating Point
John Farrier, CPPCon 2015

Slides and Source Code from "Demystifying Floating Point Numbers" given at CPPCon 2015.

##Notes
* IEEE Single Format - 1 sign, 8 exponent, 23 fraction
* IEEE Double Format - 1 sign, 11 exponent, 52 fraction
* Has several rounding options (Java only permits "Round to Nearest")
* Do not test for equality
* SSE and SSE2 may be too underspecified to be deterministic
* Order of operations matters
* The relative magnitudes of numbers involved in a computataion matters
* Prefer to multiply, add, and subtract - divide only when necessary
* Avoid elementary functions i.e. logarithmic, exponential, trigonometric, hyperbolic

##General Online References
* [A Logarithm Too Clever by Half](https://www.cs.berkeley.edu/~wkahan/LOG10HAF.TXT)
* [Anatomy of a floating point number](http://www.johndcook.com/blog/2009/04/06/anatomy-of-a-floating-point-number/)
* [Comparing Floating Point Numbers](http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm)
* [Demystifying Floating Point](http://dlang.org/d-floating-point.html)
* [Don't Store That in a Float](https://randomascii.wordpress.com/2012/02/13/dont-store-that-in-a-float/)
* [Five Tips for Floating Point Programming](http://www.codeproject.com/Articles/29637/Five-Tips-for-Floating-Point-Programming)
* [Floating Point Comparison](http://floating-point-gui.de/errors/comparison/)
* [Floating Point Determinisim](http://www.gamedev.net/topic/499435-floating-point-determinism/)
* [Floating point number representation](http://www.cprogramming.com/tutorial/floating_point/understanding_floating_point_representation.html)
* [Floating Point](http://www.tfinley.net/notes/cps104/floating.html)
* [IEEE 754 Converter](http://www.h-schmidt.net/FloatConverter/IEEE754.html)
* [IEEE floating-point exceptions in C++](http://www.johndcook.com/blog/ieee_exceptions_in_cpp/)
* [Is floating  point math deterministic?](http://blogs.msdn.com/b/shawnhar/archive/2009/03/25/is-floating-point-math-deterministic.aspx)
* [Is Your Model Susceptible to Floating-Point Errors?](http://jasss.soc.surrey.ac.uk/9/4/4.html)
* [Mathematics in Video Games](http://www.gamasutra.com/view/feature/131605/mathematics_in_videogames.php)
* [Visualizing Floats](http://www.gamasutra.com/view/feature/1965/visualizing_floats.php?print=1)
* [We Really Don't Know How To Compute!](http://www.infoq.com/presentations/We-Really-Dont-Know-How-To-Compute)
* [What Every Computer Scientist Should Know About Floating-Point Arithmetic](http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html)
* [Lecture Notes on the Status of IEEE 754, Kahan, October 1, 1997](http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF)
* [https://www.cs.fsu.edu/~engelen/courses/HPC-adv/FP.pdf](https://www.cs.fsu.edu/~engelen/courses/HPC-adv/FP.pdf)

##General Offline References
* Hacker's Delight, (Warren)
* Handbook of Floating-Point Arithmetic (Muller, Brisebarre)
* Numerical Computing with IEEE Floating Point Arithmetic (Overton)

##Compiler Specific Information
###Visual Studio
[Microsoft Visual C++ Floating-Point Optimization](https://msdn.microsoft.com/library/aa289157.aspx)
https://msdn.microsoft.com/library/aa289157.aspx

####Compiler Flags
* /fp:[precise | except[-] | fast | strict ]

###GCC
[Math_Optimization_Flags](https://gcc.gnu.org/wiki/Math_Optimization_Flags)
[Semantics of Floating Point Math in GCC](https://gcc.gnu.org/wiki/FloatingPointMath)
https://gcc.gnu.org/wiki/FloatingPointMath
https://gcc.gnu.org/wiki/Math_Optimization_Flags

####Compiler Flags
* -ffinite-math-only
* -fno-rounding-math
* -fno-signaling-nans
* -fno-signed-zeros
* -fnotrapping-math
* -frounding-math
* -fsignaling-nans
* -funsafe-math-optimizations
* -funsafe-math-optimizations
