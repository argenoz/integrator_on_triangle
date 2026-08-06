#include <stdint.h>
#include <stddef.h>

#define M_32_16 0x00ff00ff
#define M_32_8  0x0f0f0f0f
#define M_32_4  0x33333333
#define M_32_2  0x55555555


#define M_64_32 0x0000ffff0000ffff 
#define M_64_16 0x00ff00ff00ff00ff
#define M_64_8  0x0f0f0f0f0f0f0f0f
#define M_64_4  0x3333333333333333
#define M_64_2  0x5555555555555555


uint64_t e64 = -1;
uint32_t e32 = -1;
uint16_t e16 = -1;

struct coords32{
	uint32_t x,y;
};

struct coords16{
	uint16_t x,y;
};




uint32_t expa_16_32(uint16_t y)
{
	uint32_t x = y;
	x = (x|(x<<8))&M_32_16;//16 -> na dva
	x = (x|(x<<4))&M_32_8;//8 -> na dva
	x = (x|(x<<2))&M_32_4;//4->na dva
	x = (x|(x<<1))&M_32_2;//2->na dva
	return x;
}

uint32_t expa_32_64(uint32_t y)
{
	uint64_t x = y;
	x = (x<<16)&M_64_32;//32->na dva
	x = (x|(x<<8))&M_64_16;//16 -> na dva
	x = (x|(x<<4))&M_64_8;//8 -> na dva
	x = (x|(x<<2))&M_64_4;//4->na dva
	x = (x|(x<<1))&M_64_2;//2->na dva
	return x;
}

/*
#include <stdio.h>



int main()
{
uint32_t q=-1;
uint64_t p=expa_16_32(q);
int i = 32;
printf("%d %d\n",sizeof(uint32_t),sizeof(uint16_t));

do
{
	i--;
	if((q&(1<<i)))
	printf("1");
	else
	printf("0");
}
while(i);
printf("\n\n");
i = 64;
do
{
	i--;
	if((p&(1<<i)))
	printf("1");
	else
	printf("0");
}
while(i);

printf("\n\n");

}
*/

uint32_t morton16_32(uint16_t a, uint16_t b)
{
	return (expa_16_32(a) | (expa_16_32(b)<<1));
}

uint32_t morton32_64(uint32_t a, uint32_t b)
{
	return (expa_32_64(a) | (expa_32_64(b)<<1));
}
