

float ab_pairs[16][2]={
{2/3.0,1/6.0},
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



float optimuz[16]=
{729/6160.0,
729/6160.0,
-81/616.0,
81/220.0,
-859/55440.0,
-859/55440.0,
853/27720.0,
64/3465.0,
64/3465.0,
61/1980.0,
101/1260.0,
61/1980.0,
101/1260.0,
-122/693.0,
-122/693.0,
416/3465.0,};


float S_p;
float S_r;

float ax=0,a=0,
	  bx=1,by=0,
	  cx=0,cy=1;
int n_alpha=100,n_beta=100;

void set_points(float,float,float,float,float,float);
void set_n(int,int);



extern float f(float,float);

float get_Sp(){return S_p;}
float get_Sr(){return S_r;}

void integrate_p();
void integrate_r();

extern void print_p();
extern void print_r();


void work()
{
	integrate_p();
	integrate_r();
	
}

void integrate_p()
{
	float acx = ax-cx,
		   bcx = bx-cx,
		   acy = ay-cy,
		   bcy = by-cy;
	float J = acx*bcy-bcx*acy;
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
	
	S_p= S*J;
	
}

void integrate_r()
	{
		
		float acx = ax-cx,
		   bcx = bx-cx,
		   acy = ay-cy,
		   bcy = by-cy;
		float J = acx*bcy-bcx*acy;
		float S = 0;
		float d_alpha = 1.0/(n_alpha-1),alpha;
		int i = 0;
		while(i<n_alpha)
		{
				alpha = i*d_alpha;
				float tS = 0,d_beta = (1-alpha)/(n_beta-1);
				int j=0;
				float xa = acx*alpha+cx,
						ya=acy*alpha+cy;
				while(j<n_beta)
					{
						float beta = d_beta*j;
						float x = xa + bcx*beta,
							  y = ya + bcy*beta;
					    tS += f(x,y);
						j++;
					}
				S += (tS*d_beta);
				
				i++;
		}

		S_r = S*J*d_alpha;
		
	}


void set_points(float tax,float tay, float tbx, float tby, float tcx, float tcy)
	{
		ax = tax;
		ay = tay;
		bx = tbx;
		by = tby;
		cx = tcx;
		cy = tcy;
	}

void set_n(int a, int b)
	{
		n_alpha = a;
		n_beta=b;
	}





