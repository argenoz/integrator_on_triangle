



unsigned long morton(unsigned int a, unsigned int b)
{
unsigned long ans=0;
unsigned int i,l;
l = sizeof(unsigned int);
i=0;
while(i<l)
	{
	//ans =((ans<<2) |(((a&1)<<1)|(b&1)));
	ans = (ans |((((a&1)<<1)|(b&1))<<(i<<1))); 
	a = (a>>1);
	b = (b>>1);
	i++;
	}
return ans;
}


