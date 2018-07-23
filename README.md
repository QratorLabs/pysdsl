# Python bindings to Succinct Data Structure Library 2.0

The Succinct Data Structure Library ([SDSL][SDSL]) is a powerful and flexible C++11 library implementing succinct data structures. In total, the library contains the highlights of 40 [research publications][SDSLLIT]. Succinct data structures can represent an object (such as a bitvector or a tree) in space close to the information-theoretic lower bound of the object while supporting operations of the original object efficiently. The theoretical time complexity of an operation performed on the classical data structure and the equivalent succinct data structure are (most of the time) identical.

Most of examples from [SDSL cheat sheet][SDSL-CHEAT-SHEET] and [SDSL tutorial][SDSL-TUTORIAL] are implemented.

## Mutable bit-compressed vectors

Core classes (see `pysdsl.int_vector` for dict of all of them):

 * `pysdsl.IntVector(size, default_value, bit_width=64)` — dynamic bit width
 * `pysdsl.BitVector(size, default_value)` — static (fixed) bit width (1)
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

Buffer interface:

```python
In [9]: import array

In [10]: v = pysdsl.Int64Vector([1, 2, 3])

In [11]: array.array('Q', v)
Out[11]: array('Q', [1, 2, 3])
```

## Immutable compressed integer vectors

(See `pysdsl.enc_vector`):

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

Encoding values with variable length codes (see `pysdsl.variable_length_codes_vector`):

 * `VariableLengthCodesVectorEliasDelta(IntVector)`
 * `VariableLengthCodesVectorEliasGamma(IntVector)`
 * `VariableLengthCodesVectorFibonacci(IntVector)`
 * `VariableLengthCodesVectorComma2(IntVector)`
 * `VariableLengthCodesVectorComma4(IntVector)`

Encoding values with "escaping" technique (see `pysdsl.direct_accessible_codes_vector`):

 * `DirectAccessibleCodesVector(IntVector)`
 * `DirectAccessibleCodesVector8(IntVector)`,
 * `DirectAccessibleCodesVector16(IntVector)`,
 * `DirectAccessibleCodesVector63(IntVector)`,
 * `DirectAccessibleCodesVectorDP(IntVector)` — number of layers is chosen
                                                with dynamic programming
 * `DirectAccessibleCodesVectorDPRRR(IntVector)` — same but built on top of
                                                   RamanRamanRaoVector (see later)

Construction from python sequences is also supported.

## Immutable compressed bit (boolean) vectors

(See pysdsl.`all_immutable_bitvectors`)

 * `BitVectorInterLeaved64(BitVector)`
 * `BitVectorInterLeaved128(BitVector)`
 * `BitVectorInterLeaved256(BitVector)`
 * `BitVectorInterLeaved512(BitVector)` — A bit vector which interleaves the
                                          original `BitVector` with rank information
                                          (see later)
 * `SDVector(BitVector)` — A bit vector which compresses very sparse populated
                           bit vectors by representing the positions of 1 by the
                           Elias-Fano representation for
                           non-decreasing sequences
 * `RamanRamanRaoVector15(BitVector)`
 * `RamanRamanRaoVector63(BitVector)`
 * `RamanRamanRaoVector256(BitVector)` — An H₀-compressed bitvector representation.
 * `HybVector8(BitVector)`
 * `HybVector16(BitVector)` — A hybrid-encoded compressed bitvector
                              representation

See also: `pysdsl.raman_raman_rao_vectors`, `pysdsl.sparse_bit_vectors`,
`pysdsl.hybrid_bit_vectors` and `pysdsl.bit_vector_interleaved`.

## Rank and select operations on bitvectors

For bitvector `v` `rank(i)` for pattern `P` (by default `P` is a bitstring of
len 1: `1`) is the number of patterns `P` in the prefix `[0..i)` in vector `v`.

For bitvector `v` `select(i)` for pattern `P` (by default `P`=`1`) is the
position of the `i`-th occurrence of pattern `P` in vector `v`.

Create support instances for rank and/or select for different patterns via:

 * `v.init_rank()` or `v.init_rank_1()` for ranks of pattern `1`
    (e.g. the number of set bits in `v`)
 * `v.init_rank_0()` for ranks of pattern `0`
 * `v.init_rank_00()` (if supported by vector class) for ranks of pattern `00`
 * `v.init_rank_01()` (if supported by vector class) for ranks of pattern `01`
 * `v.init_rank_10()` (if supported by vector class) for ranks of pattern `10`
 * `v.init_rank_11()` (if supported by vector class) for ranks of pattern `11`
 * `v.init_support()` or `v.init_support_1()` for support of pattern `1`
    (e.g. the positions of set bits)
 * `v.init_support_0()` for ranks of pattern `0`
 * `v.init_support_00()` (if supported by vector class) for ranks of pattern `00`
 * `v.init_support_01()` (if supported by vector class) for ranks of pattern `01`
 * `v.init_support_10()` (if supported by vector class) for ranks of pattern `10`
 * `v.init_support_11()` (if supported by vector class) for ranks of pattern `11`

Once support instance `s` is created call it (`s(idx)` or `s.__call__(idx)`)
or use corresponding methods `s.rank(idx)` or `s.select(idx)` to get
the results.

`s.rank(idx)` and `s.select(idx)` are undefined if original bitvector is
mutable and was modified.


## Wavelet trees

The wavelet tree is a data structure that provides three efficient methods:

* The `[]`-operator: `wt[i]` returns the `i`-th symbol of vector for which the wavelet tree was build for.
* The rank method: `wt.rank(i, c)` returns the number of occurrences of symbol `c` in the prefix `[0..i-1]` in the vector for which the wavelet tree was build for.
* The select method: `wt.select(j, c)` returns the index `i` from `[0..size()-1]` of the `j`-th occurrence of symbol `c`.

## Comressed suffix arrays

Suffix array is a sorted array of all suffixes of a string.

SDSL supports bitcompressed and compressed suffix arrays.

Byte representaion of original IntVector should have no zero symbols in order to construct SuffixArray.

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


## Building

Requirements: static libraries for sdsl and divsufsort.

Call `pip` with binaries disabled to fetch sources and build the package:

```bash
pip install --no-binaries :all: pysdsl
```


[SDSL]: https://github.com/simongog/sdsl-lite
[SDSLLIT]: https://github.com/simongog/sdsl-lite/wiki/Literature
"Succinct Data Structure Literature"
[SDSL-CHEAT-SHEET]: https://simongog.github.io/assets/data/sdsl-cheatsheet.pdf
[SDSL-TUTORIAL]: https://simongog.github.io/assets/data/sdsl-slides/tutorial
