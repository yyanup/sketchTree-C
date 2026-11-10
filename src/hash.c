#include "common.h"
/*64-bit mixing function (the "finalizer" from SplitMix64).
*
* why hash at all? a canonical k-mer is just a number, and numbers from DNA
*are not random: k-mers strating wth A are all small, for example. if we
*kept the smallest raw values our sketch would be biased towards A-rich
*k-mers. mixing the bits frist makes evry k-mer equallylikely to be among
* the smallest, which is what MinHahs needs.
*
*every step 9xor-shift, multiply by an od number) can be undone, so the
**function is a bijection; two different k-mers can never get the same hash.
*/
uint64_t hash64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x>>30))* 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x>>27))* 0x94D049BB133111EBULL;
    return x ^ (x >> 31);

}