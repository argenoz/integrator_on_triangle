



unsigned long morton(unsigned int a, unsigned int b)
{
unsigned long ans=0;
unsigned int i;
i = sizeof(unsigned int);
do	{
	i--;
	ans =((ans<<2) |(((a&1)<<1)|(b&1)));
	a = (a>>1);
	b = (b>>1);
	}
while(i);
return ans;
}


