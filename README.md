# Python bindings to Succinct Data Structure Library 2.0

The Succinct Data Structure Library ([SDSL][SDSL]) is a powerful and flexible C++11 library implementing succinct data structures. In total, the library contains the highlights of 40 [research publications][SDSLLIT]. Succinct data structures can represent an object (such as a bitvector or a tree) in space close to the information-theoretic lower bound of the object while supporting operations of the original object efficiently. The theoretical time complexity of an operation performed on the classical data structure and the equivalent succinct data structure are (most of the time) identical.

## Bindings

Core classes:

 * `pysdsl.IntVector(size, default_value, bit_width=64)` — dynamic bit width
 * `pysdsl.BitVector(size, default_value)` — static bit width (1)
 * `pysdsl.Int8Vector(size, default_value)` — static bit width (8)
 * `pysdsl.Int16Vector(size, default_value)` — static bit width (16)
 * `pysdsl.Int32Vector(size, default_value)` — static bit width (32)
 * `pysdsl.Int64Vector(size, default_value)` — static bit width (64)


[SDSL]: https://github.com/simongog/sdsl-lite
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature "Succinct Data Structure Literature"
