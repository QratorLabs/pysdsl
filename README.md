# Python bindings to Succinct Data Structure Library 2.0

The Succinct Data Structure Library ([SDSL][SDSL]) is a powerful and flexible C++11 library implementing succinct data structures. In total, the library contains the highlights of 40 [research publications][SDSLLIT]. Succinct data structures can represent an object (such as a bitvector or a tree) in space close to the information-theoretic lower bound of the object while supporting operations of the original object efficiently. The theoretical time complexity of an operation performed on the classical data structure and the equivalent succinct data structure are (most of the time) identical.

## Bindings

Core classes:

 * `pysdsl.IntVector` --- dynamic bit size
 * `pysdsl.BitVector` --- static bit size
 * `pysdsl.Int8Vector` --- static bit size
 * `pysdsl.Int16Vector` --- static bit size
 * `pysdsl.Int32Vector` --- static bit size
 * `pysdsl.Int64Vector` --- static bit size


[SDSL]: https://github.com/simongog/sdsl-lite
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature "Succinct Data Structure Literature"
