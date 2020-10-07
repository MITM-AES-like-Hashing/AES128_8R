# MITM Attack on the 8-round AES-128 Hashing

This repository provides the implementation of the MITM attack on the 8-round AES-128 hashing.

## On precomputation:

* We emphasize that the precomputation is done once for all. In the main loop of the attack procedure, the neutral bytes are all obtained by looking up these precomputed tables.

* The degree of freedom for the forward is two bytes (2^16). However, we use one byte (2^8) only, because the degree of freedom for backward and the degree of matching are all one byte. 

* In the program, `DoFF` is the degree of freedom for the forward in bits. Setting `DoFF` to be 8, the computational complexity of the precomputation for the forward chunk is 2^32, and the complexity for backward is 2^24. In our implementation, the memory requirement is about 60 ~ 62.5 GB. 

* To test the program on a PC with less than 64 GB memory, one can set `DoFF` to be 4. Then, the memory requirement is reduced to 4 ~ 5 GB, while the computational complexity for the whole attack will be increased by a factor of 2^4. For example,  when setting `DoFF` to be 8, the expected computational complexity of finding matches on 32-bit in the first column of (#MC1, #AK1) is 2^24, while the complexity will be 2^4 x 2^24 = 2^28 if we set `DoFF` to be 4.


## On verification of the whole attack:

* To verify the computational complexity of the attack, we count the number of calls of the forward computations, backward computations, and full computations of 8-round AES. They are denoted with prefix `complexity_F`, `complexity_B`, and `complexity_M` in the program. The total complexity is computed as `(complexity_F + complexity_B)/2 + complexity_M`.

* To make the verification feasible in limited time, we find matches on `PARTIAL` bits instead of 128 bits. In the program, one can set `PARTIAL` = `32`, `40`, `48`, `56`, `64`. The concrete patterns of the partial matches are defined using `PARTIAL_MATCH_TARGET_MASK` in [defines.h](https://github.com/MITM-AES-like-Hashing/AES128_8R/blob/master/defines.h) (we give priority to the **first column** and the **first diagonal** of the state).

## On the results:

* The attacks were verified on a server with 16 cores and 128 GB memory. 
  * To find a match on 32 bits, it took 0.546833 mins (see [F32.log](https://github.com/MITM-AES-like-Hashing/AES128_8R/blob/master/F32.log));
  * To find a match on 40 bits, it took 71.9793 mins (see [F40.log](https://github.com/MITM-AES-like-Hashing/AES128_8R/blob/master/F40.log));
  * To find a match on 48 bits, it took 9038.75 mins (see [F48.log](https://github.com/MITM-AES-like-Hashing/AES128_8R/blob/master/F48.log)).
  
  It used 16 OpenMP threads in parallel. The computational complexity is the sum of that counted by each thread, and the time is the combined CPU time of all threads as returned by clock().
