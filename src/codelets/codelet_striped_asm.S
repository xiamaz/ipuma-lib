#ifdef __IPU__
// #include "poplar/AvailableVTypes.h"
// #include "poplar/TileConstants.hpp"
// #include "popops/EncodingConstants.hpp"
#include "poplar/StackSizeDefs.hpp"

#define SWASM __runCodelet_StripedSWAffineAsm

#define VOFF_C               0
#define VOFF_BG              1
#define VOFF_TS              2

#define VOFF_SIM_MATRIX      3
#define VOFF_SIM_WIDTH       4

#define VOFF_MAX_N_PER_TILE	 5

#define VOFF_GAP_INIT     	 6
#define VOFF_GAP_EXT      	 7

#define VOFF_BUF_SIZE     	 8
#define VOFF_MAX_AB       	 9

#define VOFF_A            	 10
#define VOFF_B            	 11
#define VOFF_Alen          	 12
#define VOFF_Blen          	 13

#define VOFF_SCORE        	 14
#define VOFF_MISMATCHES    	 15
#define VOFF_A_RANGE         16
#define VOFF_B_RANGE         17

#define STRIPING_WIDTH 2

#define STACK_SIZE (40)

#define SOFF_SIMWIDTH 1
#define SOFF_A  2 // a including offset_a
#define SOFF_B  3 // b including offset_b
#define SOFF_SSCORE 4 // score offset
#define SOFF_GI 5 // value of gap init
#define SOFF_BI1 6 // store 16bit packed b offset
#define SOFF_BI2 7 // store 16bit packed b offset
#define SOFF_APOS 8
#define SOFF_BPOS 9


// CAVEAT
// ABI defines m11 - sp, m10 - lr, m9 - fp
// but we will use these as normal registers!

// define variables
#define n m10
#define c m9
#define bg m8
#define i m7
#define j m6
#define a m5
#define ts m4

#define gige a7
#define ge a6
#define lastnogap a5
#define agap a4
#define prevnogap a3

.macro dbgBreak NVAL REG
	cmpeq \REG, $j, \NVAL
	brz \REG, 3f
	setzi \REG, 0xFFF
	st32 $mzero, $mzero, $mzero, 0
3:
.endm

.macro zeroCbG
  // zero C and bG
	ld32 $c, $mzero, $mvertex_base, VOFF_C
	ld32 $bg, $mzero, $mvertex_base, VOFF_BG
	// load max size of C and bG
	ld32 $m0, $mzero, $mvertex_base, VOFF_MAX_AB
	ld32 $m0, $mzero, $m0, 0
	add $m0, $m0, 1
	shrs $m0, $m0, 1
	.align 8
	{
		rpt $m0, (2f - 1f)/8 - 1;
		fnop
	}
1:
	{
		st64step $azeros, $mzero, $c+=, 1
		fnop
	}
	{
		st64step $azeros, $mzero, $bg+=, 1
		fnop
	}
2:
	.align 4
.endm

.macro loadLen VOFF MDEST
	ld32 \MDEST, $mzero, $mvertex_base, \VOFF
	ld32 \MDEST, $mzero, \MDEST, $n
.endm

.macro pushArrOffset VOFF_LEN VOFF_ARR SOFF
	ld32 $m0, $mzero, $mvertex_base, \VOFF_LEN
	add $m1, $n, 1  // take offset of array from 2*n + 1
	ld32 $m0, $mzero, $m0, $m1
	ld32 $m1, $mzero, $mvertex_base, \VOFF_ARR
	add $m0, $m1, $m0
	st32 $m0, $mzero, $sp, \SOFF
.endm

.global SWASM
.type SWASM, @function
DEF_STACK_SIZE_OWN STACK_SIZE SWASM

.section .text.SWASM
.align 4
SWASM:
	add $sp, $mworker_base, -STACK_SIZE

	setzi $n, 0

	ld32 $m0, $mzero, $mvertex_base, VOFF_SIM_WIDTH
	ld32 $m0, $mzero, $m0, 0
	st32 $m0, $mzero, $sp, SOFF_SIMWIDTH

	ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_EXT
	ldb16 $ge, $m0, $mzero, 0  // ge

	ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_INIT
	ld32 $m1, $mzero, $m0, VOFF_GAP_INIT
	st32 $m1, $mzero, $sp, SOFF_GI
	ldb16 $gige, $m0, $mzero, 0  // ge

	f16v2add $gige, $gige, $ge

compN_start:
	zeroCbG

	pushArrOffset VOFF_Blen, VOFF_B, SOFF_B
	pushArrOffset VOFF_Alen, VOFF_A, SOFF_A

	setzi $i, 0
loopB_start:
	setzi $lastnogap, 0
	ldb16 $agap, $mzero, $sp, SOFF_GI

	// load two B values into a register
	ld32 $m1, $mzero, $sp, SOFF_B
	ldz16 $m0, $m1, $i, 0
	shuf8x8lo $m0, $m0, $mzero
	ld32 $m1, $mzero, $sp, SOFF_SIMWIDTH
	mul $m0, $m0, $m1  // b[i] * simwidth
	shl $m0, $m0, 1 // access 16bit values

	ld32 $m2, $mzero, $mvertex_base, VOFF_SIM_MATRIX
	sort4x16lo $m1, $m0, $mzero
	add $m1, $m1, $m2
	st32 $m1, $mzero, $sp, SOFF_BI1
	roll16 $m1, $m0, $mzero
	add $m1, $m1, $m2
	st32 $m1, $mzero, $sp, SOFF_BI2
	// store offsets of bi into stack

//// START INNER LOOP
	ld32 $a, $mzero, $sp, SOFF_A // a into a
	// sub $j, $j, 1  // not needed as we will need to loop+1 for proper vectorization
	// PRELOAD first aval
	ld32 $m2, $mzero, $sp, SOFF_BI1
	ldz8 $m1, $mzero, $a, 0
	ldb16 $prevnogap, $mzero, $m2, $m1
	sort4x16lo $prevnogap, $prevnogap, $azero
	// PRELOAD End

	ld32 $ts, $mzero, $mvertex_base, VOFF_TS

	loadLen VOFF_Alen, $m2
	add $m2, $m2, 1
	setzi $j, 0
.align 8
// we are unable to interleave loads, we will need to find a solution:
// possible approach:
// 		interleaved load → uninterleaved store into a second array and use for next cycle
// 	{
//  		rpt $m2, (loopA_end - loopA_start)/8 - 1
//  		fnop
// 	}
loopA_start:
	dbgBreak 1, $m0
	{
		ld32 $a0, $mzero, $bg, 0
		f16v4add $a4:5, $a4:5, $a6:7  // agap, lastnogap + gi/ge
	}
	{
		ld32 $a1, $mzero, $c, 0
		f16v2max $agap, $agap, $lastnogap // agap = max()
	}
	{
		ldz8 $m1, $mzero, $a, 0 // LA: Load a[i]
		f16v4add $a0:1, $a0:1, $a6:7
	}
	{
		ldz8step $m0, $mzero, $a+=, 1 // LA: Load a[i+1]
		f16v2max $a0, $a0, $a1  // bg[j] = max()
	}
	{
		st32 $a0, $mzero, $bg, 0  // save bg[j]
		f16v2max $lastnogap, $a0, $azero  // lastnogap = max(bg[j], 0)
	}
	{
		ld32 $prevnogap, $mzero, $c, 0
		f16v2max $lastnogap, $lastnogap, $prevnogap
	}
	{
		ld32 $m2, $mzero, $sp, SOFF_BI1 // LA
		f16v2max $lastnogap, $lastnogap, $agap // lastnogap = max(lastnogap, agap)
	}
	{
		ld32 $a0, $mzero, $sp, SOFF_SSCORE
		f16v2max $lastnogap, $lastnogap, $a0 // lasstnogap = max(lastnogap, simmatrix[bi][aj])
	}
	{
		st32 $lastnogap, $mzero, $c, 1
		f16v2max $a0, $a0, $lastnogap
	}
	{
		add $c, $c, 2
		fnop
	}
	{
		add $bg, $bg, 2
		fnop
	}
	{
		ldb16 $a2, $mzero, $m2, $m1  // LA: Load simmatrix a[i]
		fnop
	}
	{
		ld32 $m2, $mzero, $sp, SOFF_BI2  // LA: Load offset2
		fnop
	}
	{
		ldb16 $a1, $mzero, $m2, $m0 // LA: Load simmatrix a[i+1]
		fnop
	}
	{
		st32step $a0, $mzero, $ts+=, 1
		sort4x16lo $a2, $a2, $a1 // LA: pack lower 16bits into this register for fp16 usage
	}
	{
		st32 $a0, $mzero, $sp, SOFF_SSCORE
		f16v2add $prevnogap, $prevnogap, $a2  // LA
	}
	add $j, $j, 1
	loadLen VOFF_Alen, $m0
	cmpult $m0, $j, $m0
	brnz $m0, loopA_start
loopA_end:
.align 4
//// END INNER LOOP

	loadLen VOFF_Blen, $m0
	add $i, $i, STRIPING_WIDTH
	cmpult $m0, $i, $m0
	brnz $m0, loopB_start
loopB_end:

compN_end:
	exitz $m15

.size SWASM, . - SWASM 

.section .data
.printfd.str: .asciz "%d: %d %d\n"

#endif