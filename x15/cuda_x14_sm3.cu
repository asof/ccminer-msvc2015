#include <stdio.h>
#include <stdint.h>
#include <memory.h>

#include "cuda_helper.h"

#define SM3_DIGEST_LENGTH	32
#define SM3_BLOCK_SIZE		64
#define SM3_CBLOCK		(SM3_BLOCK_SIZE)
#define SM3_HMAC_SIZE		(SM3_DIGEST_LENGTH)


#define cpu_to_be16(v) (((v)<< 8) | ((v)>>8))
#define cpu_to_be32(v) (((v)>>24) | (((v)>>8)&0xff00) | (((v)<<8)&0xff0000) | ((v)<<24))
#define be16_to_cpu(v) cpu_to_be16(v)
#define be32_to_cpu(v) cpu_to_be32(v)

//ROTATELEFT(X,n) = (((X)<<(n)) | ((X)>>(32-(n))))
#define ROTATELEFT(x, bits) __funnelshift_l(x, x, bits)

#define P0(x) ((x) ^  ROTATELEFT((x),9)  ^ ROTATELEFT((x),17))
#define P1(x) ((x) ^  ROTATELEFT((x),15) ^ ROTATELEFT((x),23))

#define FF0(x,y,z) ( (x) ^ (y) ^ (z))
#define FF1(x,y,z) (((x) & (y)) | ( (x) & (z)) | ( (y) & (z)))

#define GG0(x,y,z) ( (x) ^ (y) ^ (z))
#define GG1(x,y,z) (((x) & (y)) | ( (~(x)) & (z)) )


__device__
void sm3_compress_gpu(uint32_t digest[8], const uint32_t pblock[64])
{
	int j;
	uint32_t W[68], W1[64];
	//const uint32_t *pblock = (const uint32_t *)block;

	uint32_t A = digest[0];
	uint32_t B = digest[1];
	uint32_t C = digest[2];
	uint32_t D = digest[3];
	uint32_t E = digest[4];
	uint32_t F = digest[5];
	uint32_t G = digest[6];
	uint32_t H = digest[7];
	uint32_t SS1, SS2, TT1, TT2, T[64];

	for (j = 0; j < 16; j++) {
		W[j] = cpu_to_be32(pblock[j]);
	}
	for (j = 16; j < 68; j++) {
		W[j] = P1(W[j - 16] ^ W[j - 9] ^ ROTATELEFT(W[j - 3], 15)) ^ ROTATELEFT(W[j - 13], 7) ^ W[j - 6];;
	}
	for (j = 0; j < 64; j++) {
		W1[j] = W[j] ^ W[j + 4];
	}

	for (j = 0; j < 16; j++) {

		T[j] = 0x79CC4519;
		SS1 = ROTATELEFT((ROTATELEFT(A, 12) + E + ROTATELEFT(T[j], j)), 7);
		SS2 = SS1 ^ ROTATELEFT(A, 12);
		TT1 = FF0(A, B, C) + D + SS2 + W1[j];
		TT2 = GG0(E, F, G) + H + SS1 + W[j];
		D = C;
		C = ROTATELEFT(B, 9);
		B = A;
		A = TT1;
		H = G;
		G = ROTATELEFT(F, 19);
		F = E;
		E = P0(TT2);
	}

	for (j = 16; j < 64; j++) {

		T[j] = 0x7A879D8A;
		SS1 = ROTATELEFT((ROTATELEFT(A, 12) + E + ROTATELEFT(T[j], j)), 7);
		SS2 = SS1 ^ ROTATELEFT(A, 12);
		TT1 = FF1(A, B, C) + D + SS2 + W1[j];
		TT2 = GG1(E, F, G) + H + SS1 + W[j];
		D = C;
		C = ROTATELEFT(B, 9);
		B = A;
		A = TT1;
		H = G;
		G = ROTATELEFT(F, 19);
		F = E;
		E = P0(TT2);
	}

	digest[0] ^= A;
	digest[1] ^= B;
	digest[2] ^= C;
	digest[3] ^= D;
	digest[4] ^= E;
	digest[5] ^= F;
	digest[6] ^= G;
	digest[7] ^= H;
}


__global__ void x14_sm3_gpu_hash_64(uint32_t threads, uint32_t startNounce, uint32_t *g_hash)
{
    uint32_t thread = (blockDim.x * blockIdx.x + threadIdx.x);
    if (thread < threads)
    {
		uint32_t digest[8];
		unsigned char block[64]; 
		//sm3_ctx_t ctxData;
		//sm3_ctx_t *ctx = &ctxData;
		uint32_t nounce = (startNounce + thread);

		uint32_t hashPosition = nounce - startNounce;
		uint32_t* data = &g_hash[hashPosition * 16];

		memset(digest, 0, 32);
		memset(block, 0, 64);

		digest[0] = 0x7380166F;
		digest[1] = 0x4914B2B9;
		digest[2] = 0x172442D7;
		digest[3] = 0xDA8A0600;
		digest[4] = 0xA96F30BC;
		digest[5] = 0x163138AA;
		digest[6] = 0xE38DEE4D;
		digest[7] = 0xB0FB0E4E;

		sm3_compress_gpu(digest, data);

		uint32_t *pdigest = (uint32_t *)data;
		uint32_t *count = (uint32_t *)(block + SM3_BLOCK_SIZE - 8);

		block[0] = 0x80;
		memset(block + 1, 0, SM3_BLOCK_SIZE - 9);

		count[0] = 0;
		count[1] = cpu_to_be32((uint32_t)1 << 9);

		sm3_compress_gpu(digest, (uint32_t *)block);
		for (int i = 0; i < 8; i++)
			pdigest[i] = cpu_to_be32(digest[i]);
		memset(&pdigest[8], 0, 32);
		//memset(digest, 0, 32);
		//memset(block, 0, 64);
		//if (thread == 0)
		//	printf("%X %X %X %X\n", pdigest[0], pdigest[1], pdigest[2], pdigest[3]);
    }
}


__host__ 
void x14_sm3_cpu_init(int thr_id, uint32_t threads)
{
}

__host__
void x14_sm3_cpu_hash_64(int thr_id, uint32_t threads, uint32_t startNounce, uint32_t *d_nonceVector, uint32_t *d_hash, int order)
{
    const uint32_t threadsperblock = 64;

    dim3 grid((threads + threadsperblock-1)/threadsperblock);
	dim3 block(threadsperblock);

	x14_sm3_gpu_hash_64 << <grid, block >> >(threads, startNounce, d_hash);
	MyStreamSynchronize(NULL, order, thr_id);
}
