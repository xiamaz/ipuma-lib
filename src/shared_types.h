#pragma once

#ifndef NSEEDS
#define NSEEDS 2
#endif

namespace ipu {

struct __attribute__((__packed__)) SeedPair {
	int32_t seedAStartPos;	
	int32_t seedBStartPos;	
};

struct __attribute__((__packed__)) SWMeta {
	int32_t sizeA;
	int32_t offsetA;
	int32_t sizeB;
	int32_t offsetB;
};

struct __attribute__((__packed__)) XDropMeta : public SWMeta {
	int32_t seedAStartPos;	
	int32_t seedBStartPos;	
};
}