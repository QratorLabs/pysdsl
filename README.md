# Python bindings to Succinct Data Structure Library 2.0

The Succinct Data Structure Library ([SDSL][SDSL]) is a powerful and flexible C++11 library implementing succinct data structures. In total, the library contains the highlights of 40 [research publications][SDSLLIT]. Succinct data structures can represent an object (such as a bitvector or a tree) in space close to the information-theoretic lower bound of the object while supporting operations of the original object efficiently. The theoretical time complexity of an operation performed on the classical data structure and the equivalent succinct data structure are (most of the time) identical.

## Bindings

Core classes:

 * `pysdsl.IntVector(size, default_value, bit_width=64)` — dynamic bit width
 * `pysdsl.BitVector(size, default_value)` — static bit width (1)
 * `pysdsl.Int4Vector(size, default_value)` — static bit width (4)
 * `pysdsl.Int8Vector(size, default_value)` — static bit width (8)
 * `pysdsl.Int16Vector(size, default_value)` — static bit width (16)
 * `pysdsl.Int24Vector(size, default_value)` — static bit width (24)
 * `pysdsl.Int32Vector(size, default_value)` — static bit width (32)
 * `pysdsl.Int64Vector(size, default_value)` — static bit width (64)

Construction from python sequences is also supported.

```python

In [1]: import pysdsl

In [2]: %time v = pysdsl.IntVector(1024 * 1024 * 256)
CPU times: user 914 ms, sys: 509 ms, total: 1.42 s
Wall time: 1.42 s

In [3]: v.size_in_mega_bytes
Out[3]: 2048.000008583069

In [4]: %time v.set_to_id()  # like *v = range(len(v))
CPU times: user 8.19 s, sys: 1.3 ms, total: 8.19 s
Wall time: 8.19 s

In [5]: v.width
Out[5]: 64

In [6]: %time v.bit_compress()
CPU times: user 23.3 s, sys: 155 ms, total: 23.5 s
Wall time: 23.5 s

In [7]: v.width
Out[7]: 28

In [8]: v.size_in_mega_bytes
Out[8]: 896.0000085830688

```

## Immutable compressed integer vectors

 * `EncVectorEliasDelta(IntVector)`
 * `EncVectorEliasGamma(IntVector)`
 * `EncVectorFibonacci(IntVector)`
 * `EncVectorComma2(IntVector)`
 * `EncVectorComma4(IntVector)`

```python
In [9]: %time ev = pysdsl.EncVectorEliasDelta(v)
CPU times: user 26.5 s, sys: 31.8 ms, total: 26.5 s
Wall time: 26.5 s

In [10]: ev.size_in_mega_bytes
Out[10]: 45.75003242492676
```

Encoding values with variable length codes:

 * `VlcVectorEliasDelta(IntVector)`
 * `VlcVectorEliasGamma(IntVector)`
 * `VlcVectorFibonacci(IntVector)`
 * `VlcVectorComma2(IntVector)`
 * `VlcVectorComma4(IntVector)`

Encoding values with "escaping" technique:

 * `DacVector(IntVector)`
 * `DacVectorDP(IntVector)` — number of layers is choosen
                              with dynamic programming

Construction from python sequences is also supported.

## Immutable compressed bit (boolean) vectors

 * `BitVectorIL512(BitVector)` — A bit vector which interleaves the
                                 original `BitVector` with rank information.
 * `SDVector(BitVector)` — A bit vector which compresses very sparse populated
                           bit vectors by representing the positions of 1 by the
                           Elias-Fano representation for
                           non-decreasing sequences
 * `RRRVector63(BitVector)` — An H₀-compressed bitvector representation.
 * `HybVector16(BitVector)` — A hybrid-encoded compressed bitvector
                              representation

## Objects memory structure

Any object has a `.structure` property with technical information about an
object. `.structure_json` also provided for web-view implementations.
`.write_structure_json()` method puts that information into a file.

`.size_in_bytes` and `.size_in_mega_bytes` properties show how much memory the
object is occupying.

## Saving/Loading objects

All objects provide `.store_to_checked_file()` method allowing one to save
object into a file.

All classes provide `.load_from_checkded_file()` static method allowing one to
load object stored  with `.store_to_checked_file()`


[SDSL]: https://github.com/simongog/sdsl-lite
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature
"Succinct Data Structure Literature"
