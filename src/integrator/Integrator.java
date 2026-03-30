package integrator;


import extJavaLib.extNumLib.Ariphmetical;
        


        
public class Integrator {
    
    private Ariphmetical J,ax,ay,bx,by,cx,cy,acx,acy,Kx,Ky;
    
    
    
    public Integrator(Ariphmetical[][] ABC)
        {
        ax = (new Ariphmetical(ABC[0][0])).cast(10);
        ay = new Ariphmetical(ABC[0][1]);
        ay = ay.cast(10);
        bx = new Ariphmetical(ABC[1][0]);
        bx = bx.cast(10);
        by = new Ariphmetical(ABC[1][1]);
        by = by.cast(10);
        cx = new Ariphmetical(ABC[2][0]);
        cx = cx.cast(10);
        cy = new Ariphmetical(ABC[2][1]);
        cy = cy.cast(10);
        acx = Ariphmetical.sub(ax,cx);
        acy = Ariphmetical.sub(ay,cy);
        Kx = Ariphmetical.sub(bx,cx);
        Ky = Ariphmetical.sub(by,cy);
        J = Ariphmetical.sub(
                Ariphmetical.prod(acx,Ky)
                , 
                Ariphmetical.prod(acy,Kx)
                );
        }
    
    public static Ariphmetical frac(Ariphmetical n)
        {
        Ariphmetical ans = Ariphmetical.E8;    
        while(Ariphmetical.cmp(n, Ariphmetical.E8)==2)
                {
                    ans = Ariphmetical.prod(n, ans);
                    n = Ariphmetical.sub(n,Ariphmetical.E8);
                }
        return ans;
        }
    
    public static Ariphmetical C_n_k(Ariphmetical n,Ariphmetical k)
        {
        return Ariphmetical.div(frac(n), Ariphmetical.prod(frac(k), frac(Ariphmetical.sub(n,k))));
        }
    
    public Ariphmetical integrate(Ariphmetical n, Ariphmetical m)
        {
        Ariphmetical S1,S2,S3,S4,S5,N = Ariphmetical.N8,E = Ariphmetical.E8;
        n = n.cast(8);
        m = m.cast(8);
        /*
        int _i;
        for(_i=0;_i<S.length;_i++)
            {
            S[_i] = Ariphmetical.N8;
            }
        */
        S1 = N;
        Ariphmetical La1,La,theta1,theta,n_p,m_q;
        Ariphmetical p,q,s,v,u,uv1,pq1,tmp;
        p = N;
        boolean flg = true;
        while(Ariphmetical.cmp(p,n)!=2)//1 
            {
            q = N;
            S2 = q;
            
            theta1 = Ariphmetical.prod(Ariphmetical.pow(Kx, p), C_n_k(n,p));
            n_p = Ariphmetical.sub(n,p);
            while(Ariphmetical.cmp(q,m)!=2)//2
                {
                v = N;
                S3= v;
                theta = Ariphmetical.prod(Ariphmetical.pow(Ky, q), C_n_k(m,q));
                theta = Ariphmetical.prod(theta,theta1);
                
                m_q = Ariphmetical.sub(m,q);
                pq1 = Ariphmetical.sum(E,Ariphmetical.sum(p,q));
                while(Ariphmetical.cmp(v,n_p)!=2)//3
                    {
                    u = N;
                    S4= u;
                    La1 = Ariphmetical.prod(C_n_k(n_p,v), 
                            Ariphmetical.prod(
                            Ariphmetical.pow(acx,v),
                            Ariphmetical.pow(cx, Ariphmetical.sub(n_p,v)))
                            
                            
                            );
                    while(Ariphmetical.cmp(u,m_q)!=2)
                        {
                        S5 = s = N;
                        uv1 = Ariphmetical.sum(E,Ariphmetical.sum(u,v));
                        La = Ariphmetical.prod(C_n_k(m_q,u), 
                            Ariphmetical.prod(
                            Ariphmetical.pow(acy,u),
                            Ariphmetical.pow(cy, Ariphmetical.sub(m_q,u)))
                            );
                        La = Ariphmetical.prod(La,La1);
                        La = Ariphmetical.div(La, uv1);
                        flg = true;
                        while(Ariphmetical.cmp(s,uv1)!=2)
                            {
                            tmp = Ariphmetical.div(C_n_k(uv1,s).cast(10),Ariphmetical.sum(pq1,s));
                            
                            if(flg)
                                S5 = Ariphmetical.sum(S5,tmp);
                            else
                                S5 = Ariphmetical.sub(S5,tmp);
                            
                            flg = ! flg;
                            s = Ariphmetical.sum(s, E);
                            }
                        
                        S5 = Ariphmetical.prod(La, S5);
                        
                        S4 = Ariphmetical.sum(S4, S5);
                        //System.out.print(S4.cast(6)+"\n");
                        S5=null;
                        u = Ariphmetical.sum(u, E);
                        }
                    S3 = Ariphmetical.sum(S3,S4);
                    S4 = null;
                    v = Ariphmetical.sum(v, E);    
                    }
                S2 = Ariphmetical.sum(S2,Ariphmetical.prod(S3, theta));
                S3 = null;
                q = Ariphmetical.sum(q, E);
                }
            S1 = Ariphmetical.sum(S2,S1);
            S2 = null;
            p = Ariphmetical.sum(p, E);
            }
        return Ariphmetical.prod(J,S1);
        }
    
    
    
    
    
    
    public static void main(String[] args) {
        // TODO code application logic here
        //System.out.print(Ariphmetical.pow((new Ariphmetical(3)).cast(8), (new Ariphmetical(123)).cast(8))+"\n");
         Ariphmetical inti;
         Ariphmetical ax,ay,bx,by,cx,cy;
         ax = (new Ariphmetical(-1)).cast(10);
         ay = (new Ariphmetical(2)).cast(10);
         bx = (new Ariphmetical(6)).cast(10);
         by = (new Ariphmetical(2)).cast(10);
         cx = (new Ariphmetical(4)).cast(10);
         cy = (new Ariphmetical(4)).cast(10);
         Ariphmetical[][] ABC = new  Ariphmetical[3][2];
         ABC[0][0]=ax;
         ABC[0][1]=ay;
         ABC[1][0]=bx;
         ABC[1][1]=by;
         ABC[2][0]=cx;
         ABC[2][1]=cy;
         Ariphmetical n = (new Ariphmetical(2)).cast(8),
                 m = (new Ariphmetical(3)).cast(8);
         Integrator i = new Integrator(ABC);
         inti = i.integrate(n, m);
         System.out.print(inti+"\t integrale\n");
         //System.out.print(Ariphmetical.pow(ax,n));
    
    }
    
}
