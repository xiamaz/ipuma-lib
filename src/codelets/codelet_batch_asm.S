#ifdef __IPU__
// #include "poplar/AvailableVTypes.h"
// #include "poplar/TileConstants.hpp"
// #include "popops/EncodingConstants.hpp"
#include "poplar/StackSizeDefs.hpp"

#define SWASM __runCodelet_SWAffineAsm

#define VOFF_C               0
#define VOFF_BG              1

#define VOFF_SIM_MATRIX      2

#define VOFF_MAX_N_PER_TILE	 3

#define VOFF_GAP_INIT     	 4
#define VOFF_GAP_EXT      	 5

#define VOFF_BUF_SIZE     	 6
#define VOFF_MAX_AB       	 7

#define VOFF_SEQS          	 8
#define VOFF_META            9

#define VOFF_SCORE        	 10
#define VOFF_A_RANGE         11
#define VOFF_B_RANGE         12

#define VOFF_FORWARD_ONLY    13

#define SOFF_A               1
#define SOFF_B               2
#define SOFF_ALEN            3
#define SOFF_BLEN            4
#define SOFF_SSCORE          5
#define STACK_SIZE (40)

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
#define apos m4
#define bpos m3
#define bi m2

#define sscore a7
#define prevnogap a6
#define lastnogap a5
#define agap a4
#define gige a3
#define ge a2


.macro dbgBreak NVAL REG
	cmpeq \REG, $n, \NVAL
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
	ld32 $c, $mzero, $mvertex_base, VOFF_C
	ld32 $bg, $mzero, $mvertex_base, VOFF_BG
.endm

.macro pushNOffsets
	shl $m1, $n, 4 // 4 * n and 32bit addr
	ld32 $m0, $mzero, $mvertex_base, VOFF_META
	ld32step $m2, $m0, $m1+=, 1  // load a len
	st32 $m2, $mzero, $sp, SOFF_ALEN
	ld32step $m2, $m0, $m1+=, 1 // load a offset
	ld32 $m3, $mzero, $mvertex_base, VOFF_SEQS
	add $m2, $m3, $m2
	st32 $m2, $mzero, $sp, SOFF_A
	ld32step $m2, $m0, $m1+=, 1 // load b len
	st32 $m2, $mzero, $sp, SOFF_BLEN
	ld32step $m2, $m0, $m1+=, 1 // laod b offset
	add $m2, $m2, $m3
	st32 $m2, $mzero, $sp, SOFF_B
.endm

.global SWASM
.type SWASM, @function
DEF_STACK_SIZE_OWN STACK_SIZE SWASM

.section .text.SWASM
.align 4
SWASM:
	add $sp, $mworker_base, -STACK_SIZE

	ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_EXT
	ld32 $ge, $m0, $mzero, 0  // ge

	setzi $n, 0
compN_start:
 	zeroCbG

	pushNOffsets

	setzi $i, 0
	setzi $apos, 1
	setzi $bpos, 0
	setzi $sscore, 0

	// check if next sequence is zero
	// ld32 $m0, $mzero, $mvertex_base, VOFF_Blen
	ld32 $m0, $mzero, $sp, SOFF_BLEN
	brz $m0, compN_end // check if length of B is zero

loopB:
	ld32 $m0, $mzero, $sp, SOFF_B
	ldz8 $bi, $mzero, $m0, $i

{
	ld32 $bg, $mzero, $mvertex_base, VOFF_BG
	setzi $lastnogap, 0
}
{
	setzi $j, 0
	setzi $prevnogap, 0
}

	// load simmatrix[b[i]]
	ld32 $m0, $mzero, $mvertex_base, VOFF_SIM_MATRIX
	ld32 $bi, $mzero, $m0, $bi

	ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_INIT
	ld32 $agap, $m0, $mzero, 0 // gi

{
	ld32 $c, $mzero, $mvertex_base, VOFF_C
	f32add $gige, $agap, $ge
}

	st32 $sscore, $mzero, $sp, SOFF_SSCORE

	ld32 $a, $mzero, $sp, SOFF_A
	ld32 $m0, $mzero, $sp, SOFF_ALEN

{
	rpt $m0, (loopA_end - loopA_start)/8 - 1
	fnop
}
loopA_start:
	{
		ld32 $a0, $mzero, $bg, 0
		f32v2add $a4:5, $a4:5, $a2:3
	}
	{
		ld32 $a1, $mzero, $c, 0
		f32max $agap, $agap, $lastnogap // agap calculation
	}
	{
		ldz8step $m0, $mzero, $a+=, 1 // m0 = a[j]
		f32v2add $a0:1, $a0:1, $a2:3 // calc bg + ge and c[j] + gi + ge
	}
	{
		nop
		f32max $a0, $a0, $a1  // a0 = bgj
	}
	{
		st32step $a0, $mzero, $bg+=, 1 // save bgj
		f32max $lastnogap, $a0, $azero
	}
	{
		ld32 $a0, $bi, $mzero, $m0
		f32max $lastnogap, $lastnogap, $agap
	}
	{
		add $j, $j, 1
		f32add $a0, $a0, $prevnogap
	}
	{
		ld32 $prevnogap, $mzero, $c, 0
		f32max $lastnogap, $lastnogap, $a0
	}
	{
		st32step $lastnogap, $mzero, $c+=, 1
		f32cmpgt $a0, $lastnogap, $sscore
	} 
	{
		mov $m0, $a0
		f32max $sscore, $sscore, $lastnogap
	}
	{
		movz $apos, $m0, $j
		fnop
	}
loopA_end:
	ld32 $a0, $mzero, $sp, SOFF_SSCORE
	f32cmpgt $a0, $sscore, $a0
	mov $m0, $a0
	movz $bpos, $m0, $i

	add $i, $i, 1
	ld32 $m0, $mzero, $sp, SOFF_BLEN
	cmpult $m0, $i, $m0
	brnz $m0, loopB // check if length of B is zero

loopB_end:
	ld32 $m0, $mzero, $mvertex_base, VOFF_SCORE
	f32toi32 $a0, $sscore
	st32 $a0, $mzero, $m0, $n

	// store end position
	sub $apos, $apos, 1
	shl $m1, $apos, 16
	ld32 $m0, $mzero, $mvertex_base, VOFF_A_RANGE
	stm32 $m1, $m0, $n

	shl $m1, $bpos, 16
	ld32 $m0, $mzero, $mvertex_base, VOFF_B_RANGE
	stm32 $m1, $m0, $n

	// check if forward only
	ld32 $m0, $mzero, $mvertex_base, VOFF_FORWARD_ONLY
	ldz8 $m0, $mzero, $m0, 0
	brnz $m0, loopReverse_end

// now run in reverse and find startpos
	zeroCbG

	// set i and j to apos and bpos
	mov $i, $bpos
	// mov $j, $apos
	// store j so it can be restored each cycle

	st32 $apos, $sp, $mzero, 0
 	setzi $sscore, 0

	// ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_EXT
	// ld32 $ge, $m0, $mzero, 0  // ge
loopBrev:
	// loading bi

	ld32 $m0, $mzero, $sp, SOFF_B
	ldz8 $bi, $mzero, $m0, $i

	// load simmatrix[b[i]]
	ld32 $m0, $mzero, $mvertex_base, VOFF_SIM_MATRIX
	ld32 $bi, $mzero, $m0, $bi

	ld32 $bg, $mzero, $mvertex_base, VOFF_BG

	{
		ld32 $m0, $mzero, $mvertex_base, VOFF_GAP_INIT
		setzi $lastnogap, 0
	}
	{
		ld32 $agap, $m0, $mzero, 0 // gi
		setzi $prevnogap, 0
	}
	{
		ld32 $c, $mzero, $mvertex_base, VOFF_C
		f32add $gige, $agap, $ge
	}

	// loading start position of a
	ld32 $a, $mzero, $sp, SOFF_A

	// repeats is endpos + 1
	ld32 $j, $mzero, $sp, 0
	add $m0, $j, 1
	st32 $sscore, $mzero, $sp, SOFF_SSCORE

{
	rpt $m0, (loopArev_end - loopArev_start)/8 - 1
	fnop
}
loopArev_start:
	{
		ld32 $a0, $mzero, $bg, $j
		f32v2add $a4:5, $a4:5, $a2:3
	}
	{
		ld32 $a1, $mzero, $c, $j
		f32max $agap, $agap, $lastnogap
	}
	{
		ldz8 $m0, $mzero, $a, $j // m0 = a[j]
		f32v2add $a0:1, $a0:1, $a2:3 // calc bg + ge and c[j] + gi + ge
	}
	{
		nop
		f32max $a0, $a0, $a1  // a0 = bgj
	}
	{
		st32 $a0, $mzero, $bg, $j // save bgj
		f32max $lastnogap, $a0, $azero
	}
	{
		ld32 $a0, $bi, $mzero, $m0
		f32max $lastnogap, $lastnogap, $agap
	}
	{
		ld32 $prevnogap, $mzero, $c, $j
		f32add $a0, $a0, $prevnogap
	}
	{
		nop
		f32max $lastnogap, $lastnogap, $a0
	}
	{
		st32 $lastnogap, $mzero, $c, $j
		f32cmpgt $a0, $lastnogap, $sscore
	} 
	{
		mov $m0, $a0
		f32max $sscore, $sscore, $lastnogap
	}
	{
		movz $apos, $m0, $j
		fnop
	}
	{
		sub $j, $j, 1
		fnop
	}
loopArev_end:
	ld32 $a0, $mzero, $sp, SOFF_SSCORE
	f32cmpgt $a0, $sscore, $a0
	mov $m0, $a0
	movz $bpos, $m0, $i

	sub $i, $i, 1
	cmpslt $m0, $i, $mzero
	brz $m0, loopBrev // check if length of B is zero

loopBrev_end:
	ld32 $m0, $mzero, $mvertex_base, VOFF_A_RANGE
	ld32 $m1, $mzero, $m0, $n
	or $m1, $m1, $apos
	stm32 $m1, $m0, $n

	ld32 $m0, $mzero, $mvertex_base, VOFF_B_RANGE
	ld32 $m1, $mzero, $m0, $n
	or $m1, $m1, $bpos
	stm32 $m1, $m0, $n

loopReverse_end:
	// increment n and check if we run again
	add $n, $n, 1
	ld32 $m0, $mzero, $mvertex_base, VOFF_MAX_N_PER_TILE
	ld32 $m0, $mzero, $m0, 0
	cmpult $m0, $n, $m0 // n < max_ab
	brnz $m0, compN_start
compN_end:
	exitz $m15

.size SWASM, . - SWASM 

.section .data
.printfd.str: .asciz "%d: %d %d\n"

#endif