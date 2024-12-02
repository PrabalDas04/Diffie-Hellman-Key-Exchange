This is a project on Diffie-Hellman key Exchange in Multiplicative and Elliptic Curve group on field Z_p.
I have implemented the finite field arithmetics on finite field which is base of DHKEX.
During key exchange each party will do those computation locally prior to the key exchange.
The whole code is implemented in basic C language with minimum no of library used. Refer to documentation pdf for better understanding.
In the code, "Barrett_Reduction.c" contains all the function needed for Barrett reduction as well as some common function which is needed for both case operation.
"ECC_Scalar_Mult.c" contains all the function in elliptic curve arithmetic like point doubling, point addition, scalar multiple, inverse of a field element using Fermat etc.
"Mult_Expo.c" contains all the function needed for arithmetic in multiplicative group like left to right multiplication, right to left, Montgomery with branching and without branching etc.
