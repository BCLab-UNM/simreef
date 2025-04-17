
// bytell_hash_map.hpp — High-Performance Hash Map
// -------------------------------------------------
// This header provides an implementation of a fast and memory-efficient hash map.
// It is likely derived from ankerl::unordered_dense or a similar library,
// which uses open addressing with linear probing and robin-hood hashing.
//
// Advantages:
// - Extremely fast lookup and insertion
// - Smaller memory footprint than std::unordered_map
// - Does not require a separate hash table allocation for each node
//
// Typical usage in SimForager:
// - Used in place of unordered_map for grid-point indexing or concentration updates
//   where hash collisions are minimal and performance is critical.
//
// API Overview (typical):
// - unordered_dense::map<Key, Value>
// - Supports custom hash functions
// - Compatible with most STL containers and algorithms

// Note: As this is a third-party utility header, it should not be modified
// unless extending or debugging the library itself. Performance-critical.
