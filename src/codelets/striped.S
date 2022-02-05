	.text
	.allow_optimizations
	.file	"3"
	.section	.text.__runCodelet_MultiSWStriped,"ax",@progbits
	.globl	__runCodelet_MultiSWStriped     # -- Begin function __runCodelet_MultiSWStriped
	.p2align	2
	.type	__runCodelet_MultiSWStriped,@function
__runCodelet_MultiSWStriped:            # @__runCodelet_MultiSWStriped
.Lfunc_begin0:
# %bb.0:
	mov	$m11, $m12
	get $m0, 1
	and $m1, $m0, 7
	mov	$m0, $m13
	call $m10, _ZN14MultiSWStriped7computeEj
	exitnz $m0
.Lfunc_end0:
	.size	__runCodelet_MultiSWStriped, .Lfunc_end0-__runCodelet_MultiSWStriped
	.section	.stack_sizes,"o",@progbits,.text.__runCodelet_MultiSWStriped
	.long	.Lfunc_begin0
	.byte	0
	.section	.text.__runCodelet_MultiSWStriped,"ax",@progbits
                                        # -- End function
	.section	.text._ZN14MultiSWStriped7computeEj,"axG",@progbits,_ZN14MultiSWStriped7computeEj,comdat
	.weak	_ZN14MultiSWStriped7computeEj   # -- Begin function _ZN14MultiSWStriped7computeEj
	.p2align	2
	.type	_ZN14MultiSWStriped7computeEj,@function
_ZN14MultiSWStriped7computeEj:          # @_ZN14MultiSWStriped7computeEj
.Lfunc_begin1:
# %bb.0:
	add $m11, $m11, -112
	sub $m6, $m12, $m11
	cmpult $m6, $m6, __worker_stack_max_growth_plus_one
	brz $m6, _stackoverflow
	st32 $m8, $m11, $m15, 4                 # 4-byte Folded Spill
	st32 $m9, $m11, $m15, 3                 # 4-byte Folded Spill
	st32 $m10, $m11, $m15, 2                # 4-byte Folded Spill
	st32 $m7, $m11, $m15, 1                 # 4-byte Folded Spill
	mov	$m4, $m0
	ld32 $m0, $m4, $m15, 5
	ld32 $m0, $m0, $m15, 0
	cmpslt $m0, $m0, 1
	brnz $m0, .LBB1_19
# %bb.1:
	mov	$m6, $m15
	ld32 $m0, $m4, $m15, 6
	ld32 $m2, $m0, $m15, 0
	ld32 $m0, $m4, $m15, 7
	ld32 $m0, $m0, $m15, 0
	st32 $m0, $m11, $m15, 27                # 4-byte Folded Spill
	add $m0, $m0, $m2
	st32 $m0, $m11, $m15, 20                # 4-byte Folded Spill
	shl $m0, $m1, 2
	ld32 $m2, $m4, $m15, 10
	st32 $m2, $m11, $m15, 5                 # 4-byte Folded Spill
	ld32 $m2, $m4, $m15, 11
	st32 $m2, $m11, $m15, 18                # 4-byte Folded Spill
	sub $m8, 5, $m1
	sub $m2, 0, $m1
	st32 $m2, $m11, $m15, 13                # 4-byte Folded Spill
	st32 $m0, $m11, $m15, 6                 # 4-byte Folded Spill
	sub $m0, 0, $m0
	st32 $m0, $m11, $m15, 11                # 4-byte Folded Spill
	add $m0, $m1, -1
	st32 $m0, $m11, $m15, 10                # 4-byte Folded Spill
	st32 $m8, $m11, $m15, 7                 # 4-byte Folded Spill
	st32 $m4, $m11, $m15, 23                # 4-byte Folded Spill
	st32 $m1, $m11, $m15, 8                 # 4-byte Folded Spill
	bri .LBB1_3
.LBB1_2:                                #   in Loop: Header=BB1_3 Depth=1
	ld32 $m0, $m4, $m15, 5
	ld32 $m6, $m11, $m15, 9                 # 4-byte Folded Reload
	add $m6, $m6, 1
	ld32 $m0, $m0, $m15, 0
	cmpslt $m0, $m6, $m0
	ld32 $m8, $m11, $m15, 7                 # 4-byte Folded Reload
	brz $m0, .LBB1_19
.LBB1_3:                                # =>This Loop Header: Depth=1
                                        #     Child Loop BB1_7 Depth 2
                                        #     Child Loop BB1_11 Depth 2
                                        #       Child Loop BB1_12 Depth 3
	ld32 $m2, $m4, $m15, 12
	shl $m3, $m6, 3
	ld32 $m0, $m2, $m3, 0
	st32 $m0, $m11, $m15, 21                # 4-byte Folded Spill
	brz $m0, .LBB1_19
# %bb.4:                                #   in Loop: Header=BB1_3 Depth=1
	ld32 $m5, $m4, $m15, 13
	ld32 $m0, $m5, $m3, 0
	st32 $m0, $m11, $m15, 16                # 4-byte Folded Spill
	brz $m0, .LBB1_19
# %bb.5:                                #   in Loop: Header=BB1_3 Depth=1
	st32 $m6, $m11, $m15, 9                 # 4-byte Folded Spill
	ld32 $m7, $m4, $m15, 9
	or $m3, $m3, 4
	ld32 $m2, $m2, $m3, 0
	ld32 $m0, $m5, $m3, 0
	st32 $m0, $m11, $m15, 15                # 4-byte Folded Spill
	ld32 $m3, $m7, $m15, 0
	cmpslt $m3, $m1, $m3
	brz $m3, .LBB1_8
# %bb.6:                                #   in Loop: Header=BB1_3 Depth=1
	ld32 $m3, $m4, $m15, 0
	ld32 $m5, $m4, $m15, 1
	ld32 $m0, $m11, $m15, 6                 # 4-byte Folded Reload
	add $m3, $m3, $m0
	add $m5, $m5, $m0
	mov	$m6, $m1
.LBB1_7:                                #   Parent Loop BB1_3 Depth=1
                                        # =>  This Inner Loop Header: Depth=2
	st32step $a15, $m15, $m3+=, 6
	st32step $a15, $m15, $m5+=, 6
	add $m6, $m6, 6
	ld32 $m0, $m7, $m15, 0
	cmpslt $m0, $m6, $m0
	brnz $m0, .LBB1_7
.LBB1_8:                                #   in Loop: Header=BB1_3 Depth=1
	st32 $m15, $m11, $m15, 24               # 4-byte Folded Spill
	ld32 $m0, $m11, $m15, 16                # 4-byte Folded Reload
	add $m0, $m8, $m0
	st32 $m0, $m11, $m15, 14                # 4-byte Folded Spill
	cmpslt $m3, $m1, $m0
	brz $m3, .LBB1_17
# %bb.9:                                #   in Loop: Header=BB1_3 Depth=1
	ld32 $m0, $m11, $m15, 21                # 4-byte Folded Reload
	add $m3, $m8, $m0
	mov	$m8, $m0
	ld32 $m0, $m11, $m15, 13                # 4-byte Folded Reload
	cmpslt $m3, $m0, $m3
	ld32 $m1, $m11, $m15, 20                # 4-byte Folded Reload
	brz $m3, .LBB1_17
# %bb.10:                               #   in Loop: Header=BB1_3 Depth=1
	st32 $m15, $m11, $m15, 24               # 4-byte Folded Spill
	ld32 $m0, $m11, $m15, 5                 # 4-byte Folded Reload
	add $m0, $m0, $m2
	st32 $m0, $m11, $m15, 19                # 4-byte Folded Spill
	add $m0, $m8, 5
	st32 $m0, $m11, $m15, 12                # 4-byte Folded Spill
	ld32 $m5, $m11, $m15, 8                 # 4-byte Folded Reload
.LBB1_11:                               #   Parent Loop BB1_3 Depth=1
                                        # =>  This Loop Header: Depth=2
                                        #       Child Loop BB1_12 Depth 3
	mov	$m3, $m15
	st32 $m15, $m11, $m15, 25               # 4-byte Folded Spill
	ld32 $m2, $m4, $m15, 6
	ld32 $m7, $m2, $m15, 0
	ld32 $m9, $m11, $m15, 12                # 4-byte Folded Reload
	add $m9, $m9, -1
	ld32 $m0, $m11, $m15, 16                # 4-byte Folded Reload
	cmpslt $m0, $m5, $m0
	st32 $m0, $m11, $m15, 26                # 4-byte Folded Spill
	ld32 $m0, $m11, $m15, 15                # 4-byte Folded Reload
	st32 $m5, $m11, $m15, 17                # 4-byte Folded Spill
	add $m0, $m5, $m0
	st32 $m0, $m11, $m15, 22                # 4-byte Folded Spill
	ld32 $m6, $m11, $m15, 10                # 4-byte Folded Reload
	ld32 $m2, $m11, $m15, 11                # 4-byte Folded Reload
	ld32 $m5, $m11, $m15, 13                # 4-byte Folded Reload
.LBB1_12:                               #   Parent Loop BB1_3 Depth=1
                                        #     Parent Loop BB1_11 Depth=2
                                        # =>    This Inner Loop Header: Depth=3
	shr $m0, $m6, 31
	cmpslt $m4, $m5, $m8
	ld32 $m10, $m11, $m15, 26               # 4-byte Folded Reload
	add $m0, $m0, $m10
	add $m0, $m0, $m4
	add $m4, $m1, $m3
	ld32 $m10, $m11, $m15, 27               # 4-byte Folded Reload
	add $m7, $m7, $m10
	cmpeq $m0, $m0, 3
	max $m7, $m4, $m7
	brz $m0, .LBB1_13
# %bb.14:                               #   in Loop: Header=BB1_12 Depth=3
	ld32 $m10, $m11, $m15, 23               # 4-byte Folded Reload
	ld32 $m0, $m10, $m15, 1
	ld32 $m3, $m10, $m15, 0
	ld32 $m4, $m0, $m2, 0
	add $m4, $m1, $m4
	ld32 $m1, $m3, $m2, 0
	ld32 $m8, $m11, $m15, 27                # 4-byte Folded Reload
	add $m1, $m1, $m8
	max $m1, $m4, $m1
	st32 $m1, $m3, $m2, 0
	ld32 $m3, $m11, $m15, 19                # 4-byte Folded Reload
	ldz8 $m3, $m3, $m15, $m5
	ld32 $m4, $m10, $m15, 4
	ld32 $m8, $m11, $m15, 21                # 4-byte Folded Reload
	ld32 $m3, $m4, $m15, $m3
	ld32 $m4, $m11, $m15, 18                # 4-byte Folded Reload
	ld32 $m10, $m11, $m15, 22               # 4-byte Folded Reload
	ldz8 $m4, $m4, $m15, $m10
	ld32 $m3, $m3, $m15, $m4
	ld32 $m4, $m11, $m15, 25                # 4-byte Folded Reload
	add $m3, $m3, $m4
	max $m3, $m3, $m7
	max $m1, $m3, $m1
	max $m3, $m1, 0
	ld32 $m1, $m11, $m15, 20                # 4-byte Folded Reload
	ld32 $m4, $m0, $m2, 0
	st32 $m4, $m11, $m15, 25                # 4-byte Folded Spill
	st32 $m3, $m0, $m2, 0
	ld32 $m0, $m11, $m15, 24                # 4-byte Folded Reload
	max $m0, $m0, $m3
	st32 $m0, $m11, $m15, 24                # 4-byte Folded Spill
	bri .LBB1_15
.LBB1_13:                               #   in Loop: Header=BB1_12 Depth=3
	#APP
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop

	#NO_APP
.LBB1_15:                               #   in Loop: Header=BB1_12 Depth=3
	add $m5, $m5, 1
	add $m2, $m2, 4
	add $m6, $m6, -1
	brnzdec $m9, .LBB1_12
# %bb.16:                               #   in Loop: Header=BB1_11 Depth=2
	ld32 $m5, $m11, $m15, 17                # 4-byte Folded Reload
	add $m5, $m5, 6
	ld32 $m0, $m11, $m15, 14                # 4-byte Folded Reload
	cmpslt $m0, $m5, $m0
	ld32 $m4, $m11, $m15, 23                # 4-byte Folded Reload
	brnz $m0, .LBB1_11
.LBB1_17:                               #   in Loop: Header=BB1_3 Depth=1
	ld32 $m0, $m4, $m15, 14
	setzi $m2, 4
	ld32 $m1, $m11, $m15, 8                 # 4-byte Folded Reload
	ld32 $m3, $m11, $m15, 24                # 4-byte Folded Reload
	stm32 $m3, $m0, $m1
	cmpult $m0, $m2, $m1
	ld32 $m2, $m4, $m15, 16
	stm32 $m3, $m2, $m1
	brnz $m0, .LBB1_2
# %bb.18:                               #   in Loop: Header=BB1_3 Depth=1
	#APP
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop
	nop

	#NO_APP
	bri .LBB1_2
.LBB1_19:
	setzi $m0, 1
	ld32 $m7, $m11, $m15, 1                 # 4-byte Folded Reload
	ld32 $m10, $m11, $m15, 2                # 4-byte Folded Reload
	ld32 $m9, $m11, $m15, 3                 # 4-byte Folded Reload
	ld32 $m8, $m11, $m15, 4                 # 4-byte Folded Reload
	add $m11, $m11, 112
	br $m10
.Lfunc_end1:
	.size	_ZN14MultiSWStriped7computeEj, .Lfunc_end1-_ZN14MultiSWStriped7computeEj
	.section	.stack_sizes,"Go",@progbits,_ZN14MultiSWStriped7computeEj,comdat,.text._ZN14MultiSWStriped7computeEj
	.long	.Lfunc_begin1
	.byte	112
	.section	.text._ZN14MultiSWStriped7computeEj,"axG",@progbits,_ZN14MultiSWStriped7computeEj,comdat
                                        # -- End function
	.ident	"clang version 14.0.0 (ssh://git@phabricator.sourcevertex.net/diffusion/LLVMPROJECT/llvm-project.git ad697584ed527786ca41847bf56b2019532cd738)"
	.section	".note.GNU-stack","",@progbits
	.addrsig
