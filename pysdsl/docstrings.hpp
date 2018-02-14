#pragma once


const char* doc_int_vector(
    "This generic vector class could be used to generate a vector that "
    "contains integers of fixed width `w` in [1..64]."
);

const char* doc_bit_compress(
    "Bit compress the int_vector. Determine the biggest value X "
    "and then set the int_width to the smallest possible so that "
    "we still can represent X."
);

const char* doc_dac_vector(
    "A generic immutable space-saving vector class for unsigned integers.\n"
    "The values of a dac_vector are immutable after the constructor call.\n"
    "The `escaping` technique is used to encode values.\n"
    "This is defined as follows (see [1]):\n"
    "A k-bit integer is split into `K=k/(b-1)` bits each and "
    "encoded into `K` blocks of `b` bits each. All but the last block "
    "are marked with by a 1 in the most significant bit. Escaping with "
    "`b=8` is also known as vbyte-coding (see [2]). A experimental study "
    "of using escaping for the LCP array is given in [3].\n"
    "Time complexity: Order{log n/b} worst case, where b is the number "
    "of bits in a block\nReferences:\n"
    "[1] F. Transier and P. Sanders: `Engineering Basic Search Algorithms "
    "of an In-Memory Text Search Engine`, ACM Transactions on "
    "Information Systems, Vol. 29, No.1, Article 2, 2010\n"
    "[2] H.E. Williams and J. Zobel: `Compressing integers for fast file "
    "access`, Computing Journal Vol 43, No.3, 1999\n"
    "[3] N. Brisboa, S. Ladra, G. Navarro: `Directly addressable "
    "variable-length codes'', Proceedings of SPIRE 2009."
);

const char* doc_dac_vector_dp(
    "A generic immutable space-saving vector class for unsigned integers.\n"
    "The values of a dac_vector are immutable after the constructor call.\n"
    "The \"escaping\" technique is used to encode values. Bit widths of "
    "each encoding level are chosen optimally via dynamic programming.\n"
    "References\n [1] N. Brisaboa and S. Ladra and G. Navarro: `DACs: "
    "Bringing Direct Access to Variable-Length Codes`, "
    "Information Processing and Management (IPM) 2013"
);

const char* doc_bit_vector_il(
    "A bit vector which interleaves the original bit_vector with rank "
    "information. \nThis class is a uncompressed bit vector "
    "representation. It copies the original bit_vector and interleaves "
    "the data every t_bs bits with a cumulative sum of set bits before "
    "the current position. Each cumulative sum is stored in a 64 bit "
    "word."
);

const char* doc_rrr_vector(
    "An H_0-compressed bitvector representation.\n"
    "References:\n"
    "- Rasmus Pagh, Low redundancy in dictionaries with O(1) worst "
    "case lookup time, Technical Report 1998. "
    "ftp://ftp.cs.au.dk/BRICS/Reports/RS/98/28/BRICS-RS-98-28.pdf, "
    "Section 2.\n"
    "- Rajeev Raman, V. Raman and S. Srinivasa Rao, Succinct Indexable "
    "Dictionaries with Applications to representations of k-ary trees "
    "and multi-sets. SODA 2002.\n"
    "- Francisco Claude, Gonzalo Navarro: Practical Rank/Select "
    "Queries over Arbitrary Sequences. SPIRE 2008: 176-187\n"
    "- On the fly-decoding and encoding was discovered in; Gonzalo "
    "Navarro, Eliana Providel: Fast, Small, Simple Rank/Select on "
    "Bitmaps. SEA 2012"
);

const char* doc_sd_vector(
    "A bit vector which compresses very sparse populated bit vectors "
    "by representing the positions of 1 by the Elias-Fano "
    "representation for non-decreasing sequences\n"
    "References:\n"
    "- P. Elias: ''Efficient storage and retrieval by content and "
    "address of static files'', Journal of the ACM, 1974\n"
    "- R. Fano: ''On the number of bits required to implement an "
    "associative memory'', Memorandum 61. Computer Structures Group, "
    "Project MAC, MIT, 1971\n"
    "- D. Okanohara, K. Sadakane: ''Practical Entropy-Compressed "
    "Rank/Select Dictionary'', Proceedings of ALENEX 2007."
);

const char* doc_hyb_vector(
    "A hybrid-encoded compressed bitvector representation\n"
    "References:\n- Juha Karkkainen, Dominik Kempa and "
    "Simon J. Puglisi. ''Hybrid Compression of Bitvectors for the "
    "FM-Index.'' DCC 2014."
);

const char* doc_rank_v(
    "A rank structure proposed by Sebastiano Vigna\nSpace complexity: "
    "0.25n for a bit vector of length n bits.\n\nThe superblock size is "
    "512. Each superblock is subdivided into 512/64 = 8 blocks. "
    "So absolute counts for the superblock add 64/512 bits on top of each "
    "supported bit. Since the first of the 8 relative count values is 0, "
    "we can fit the remaining 7 (each of width log(512)=9) in a 64bit "
    "word. The relative counts add another 64/512 bits on top of each "
    "supported bit.\nIn total this results in 128/512=25% overhead.\n"
    "Reference\nSebastiano Vigna: Broadword Implementation of Rank/Select "
    "Queries. WEA 2008: 154-168"
);

const char* doc_rank_v5(
    "A class supporting rank queries in constant time.\n"
    "Space complexity: 0.0625n bits for a bit vector of length n bits.\n\n"
    "The superblock size is 2048. Each superblock is subdivided into "
    "2048/(6*64) = 5 blocks (with some bit remaining). So absolute counts "
    "for the superblock add 64/2048 bits on top of each supported bit. "
    "Since the first of the 6 relative count values is 0, we can fit the "
    "remaining 5 (each of width log(2048)=11) in a 64 bit word. The "
    "relative counts add another 64/2048 bits bits on top of each "
    "supported bit. In total this results in 128/2048=6.25% overhead."
);

const char* doc_rank_scan(
    "A class supporting rank queries in linear time.\n"
    "Space complexity: Constant.\n"
    "Time complexity: Linear in the size of the supported vector."
);

const char* doc_select_mcl(
    "A class supporting constant time select queries.\n"
    "Space usage: The space usage of the data structure depends on the "
    "number `m` of ones in the original bitvector `b`. We store the "
    "position of every 4096th set bit (called L1-sampled bits) of `b`. "
    "This takes in the worst case (m/4096) log(n) <= (n/64) bits.\n"
    "Next,\n(1) if the distance of two adjacent L1-sampled bits "
    "b[i] and b[j] is greater or equal than log^4 (n), then we store "
    "each of the 4096 positions of the set `b` in [i..j-1] with "
    "log(n) bits. This results in at most "
    "`4096 log(n) / log^4(n)=4096 / log^3(n)` bits per bit.\nFor a "
    "bitvector of 4GB, i.e. log(n) = 35 we get about 0.01 bits per bit.\n"
    "If the j-i+1 < log^4(n) then\n(2) we store the relative position of "
    "every 64th set bit (called L2-sampled bits) in b[i..j-1] in at most "
    "4 log log (n) bits per L2-sampled bits.\nAn pessimistic upper bound "
    "for the space would be `4 log log (n) / 64 <= 24/64 = 0.375` bit per "
    "bit (since `log log (n) <= 6`. It is very pessimistic, since we store "
    "the relative position in `log log (j-i+1) <= log log (n)` bits.\n\n"
    "The implementation is a practical variant of the following reference:"
    "\nDavid Clark: PhD Thesis: Compact Pat Trees, University of Waterloo, "
    "1996 (Section 2.2.2). "
    "http://www.nlc-bnc.ca/obj/s4/f2/dsk3/ftp04/nq21335.pdf"
);

const char* doc_select_scan(
    "A class supporting linear time select queries.\n"
    "Space complexity: Constant\n"
    "Time complexity: Linear in the size of the supported vector."
);