
float ax,ay,bx,by,cx,cy,acx,acy,bcx,bcy,J,
	alpha_beta[3][2]={{0,1/2.0},{1/2.0,0},{1/2.0,1/2.0}};


extern float f(float x, float y);



float integrate_p3()
	{
		
		return (J*(f(ax,ay)+f(bx,by)+f(cx,cy)))/6;
	}



float integrate_p6()
	{
		
		float f4,f5,f6;
		f4 = f(acx*alpha_beta[0][0]+bcx*alpha_beta[0][1]+cx,
				acy*alpha_beta[0][0]+bcy*alpha_beta[0][1]+cy);
		f5 = f(acx*alpha_beta[1][0]+bcx*alpha_beta[1][1]+cx,
				acy*alpha_beta[1][0]+bcy*alpha_beta[1][1]+cy);		
		f6 = f(acx*alpha_beta[2][0]+bcx*alpha_beta[2][1]+cx,
				acy*alpha_beta[2][0]+bcy*alpha_beta[2][1]+cy);
		return (f4+f5+f6)*(J/6);
	}



void setV(float ax_,float ay_,float bx_,float by_,float cx_,float cy_)
{
	ax = ax_;
	ay = ay_;
	bx = bx_;
	by = by_;
	cx = cx_;
	cy = cy_;
	acx = (ax-cx);
	acy =(ay-cy);
	bcx = (bx-cx);
	bcy =(by-cy);
	J = acx*bcy-acy*bcx;
}



float optimuz[16]={729/6160.0,729/6160.0,-81/616.0,81/220.0,-859/55440.0,-859/55440.0,853/27720.0,64/3465.0,64/3465.0,61/1980.0,101/1260.0,61/1980.0,101/1260.0,-122/693.0,-122/693.0,416/3465.0,};




float ab_pairs[16][2]={{2/3.0,1/6.0},
{1/6.0,2/3.0},
{1/6.0,1/6.0},
{1/3.0,1/3.0},
{1/2.0,0.0},
{0.0,1/2.0},
{1/2.0,1/2.0},
{1/4.0,3/4.0},
{3/4.0,1/4.0},
{0.0,3/4.0},
{0.0,1/4.0},
{3/4.0,0.0},
{1/4.0,0.0},
{1/4.0,1/2.0},
{1/2.0,1/4.0},
{1/4.0,1/4.0}};

float integrate_p16()
{
	int i,j;
	float S;
	i=15;
	float alpha = ab_pairs[i][0],
				beta = ab_pairs[i][1];
	
	float x = acx*alpha+bcx*beta+cx,
				  y = acy*alpha+bcy*beta+cy;
	
	S = optimuz[i]*f(x,y);
	i=15;
	do
		{
			i--;
			alpha = ab_pairs[i][0];
			beta = ab_pairs[i][1];
			x = acx*alpha+bcx*beta+cx;
			y = acy*alpha+bcy*beta+cy;
			S+=(optimuz[i]*f(x,y));
		}
	while(i);
	
	return S*J;
	
}




