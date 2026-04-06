from copy import deepcopy as DPC
import sympy as SP


class Integrator:
    def fakto(self,n):
        if n < 2:
            return 1
        else:
            ans = 1
            while n !=1:
                ans = ans * n
                n = n - 1
            return ans
    #def C_n_k(self,n,k):
    #    return self.fakto(n)//(self.fakto(k)*self.fakto(n-k))
    def Cnk(self,n,k):
        while n >=len(self.Cnk_):
            qwe = self.Cnk_[-1]
            l = range(len(qwe)+1)
            ewq = [1 for i in l]
            self.Cnk_.append(ewq)
            l = l[1:-1:]
            for i in l:
                ewq[i]=qwe[i-1]+qwe[i]
                
            #print(len(self.Cnk_),n," len, n")
        return self.Cnk_[n][k]
    def C_n_k(self,n,k):
        #nuzhn=0
        #return self.fakto(n)//(self.fakto(k)*self.fakto(n-k))
        return self.Cnk(n,k)
    def __init__(self,abc):
        self.abc=DPC(abc)
        self.ax = abc[0][0]
        self.ay = abc[0][1]
        self.bx = abc[1][0]
        self.by = abc[1][1]
        self.cx = abc[2][0]
        self.cy = abc[2][1]
        self.Kx = self.bx-self.cx
        self.Ky = self.by-self.cy
        self.acx = self.ax-self.cx
        self.acy = self.ay-self.cy
        self.J = (self.ax-self.cx)*(self.by-self.cy)-(self.ay-self.cy)*(self.bx-self.cx)
        self.Cnk_=[[1],[1,1]]
    def integrate(self,n,m):
        S=[0 for i in range(5)]
        ln,lm = range(n+1),range(m+1)
        for p in ln:
            S[1]=0
            n_p = n - p
            theta1 = (self.Kx**p)*self.C_n_k(n,p)
            l_n_p = range(n_p+1)
            for q in lm:
                S[2]=0
                theta = theta1*(self.Ky**q)*self.C_n_k(m,q)
                m_q = m - q
                l_m_q =range(m_q+1)
                
                pq1 = p + q + 1.0                
                for v in l_n_p:
                    S[3]=0
                    La1 = self.C_n_k(n_p,v)*self.acx**v*self.cx**(n_p-v)
                    
                    for u in l_m_q:
                        uv1 = v+u+1
                        S[4]=0
                        La = La1*self.C_n_k(m_q,u)*self.acy**u*self.cy**(m_q-u)/(float(uv1))
                        
                        for s in range(uv1+1):
                            tmp = self.C_n_k(uv1,s)/(pq1+s)
                            
                            if s%2==0:
                                S[4]=S[4]+tmp
                            else:
                                S[4]=S[4]-tmp
                        S[3]=S[3]+S[4]*La
                        
                        S[4]=None
                    S[2]=S[2]+S[3]
                    
                    S[3]=None
                S[1]=S[1]+S[2]*theta
                
                S[2]=None
                thta=None
            S[0]=S[0]+S[1]
            #print(S[1],"\tS2")
            S[1]=None
        return S[0]*self.J

'''                        
ABC = [[-1,-31/2.0],[6, 2], [4,  4/19.0]]
ABC = [ [float(ABC[i][j]) for j in range(2)] for i in range(3)]       
qwe = Integrator(ABC)
n = 4
m = 3
x = SP.Symbol('x')
y = SP.Symbol('y')
i_ =  Integrator(ABC).integrate(n,m)
print(i_)


f = x**n * y**m

alpha = SP.Symbol('alpha')
beta = SP.Symbol('beta')
print(f)

f = f.subs('x',alpha*(ABC[0][0]-ABC[2][0])+beta*(ABC[1][0]-ABC[2][0])+ABC[2][0])
f = f.subs('y',alpha*(ABC[0][1]-ABC[2][1])+beta*(ABC[1][1]-ABC[2][1])+ABC[2][1])
i = SP.integrate(f,alpha)
i = i.subs(alpha,1-beta)-i.subs(alpha,0)
i = SP.integrate(i,beta)
i = i.subs(beta,1)-i.subs(beta,0)
J = ((ABC[0][0]-ABC[2][0])*(ABC[1][1]-ABC[2][1])-(ABC[0][1]-ABC[2][1])*(ABC[1][0]-ABC[2][0]))
i = i*J
print(i,"\tsymbolic drob")
print(float(i),"symbolic int")

print(abs(float(i)-i_))
'''
'''
ABC=[[2,1],[1,2],[0,0]]
n=12
m=12
x = SP.Symbol('x')
y = SP.Symbol('y')
f = x**n * y**m
alpha = SP.Symbol('alpha')
beta = SP.Symbol('beta')
f = f.subs('x',alpha*(ABC[0][0]-ABC[2][0])+beta*(ABC[1][0]-ABC[2][0])+ABC[2][0])
f = f.subs('y',alpha*(ABC[0][1]-ABC[2][1])+beta*(ABC[1][1]-ABC[2][1])+ABC[2][1])
i = SP.integrate(f,alpha)
i = i.subs(alpha,1-beta)-i.subs(alpha,0)
i = SP.integrate(i,beta)
i = i.subs(beta,1)-i.subs(beta,0)
J = ((ABC[0][0]-ABC[2][0])*(ABC[1][1]-ABC[2][1])-(ABC[0][1]-ABC[2][1])*(ABC[1][0]-ABC[2][0]))
i = i*J
print(i,"\tsymbolic drob")
print(float(i),"symbolic int")
#qwe = Integrator(ABC)
i_ =  Integrator(ABC).integrate(n,m)
print(i_)
#print(qwe.Cnk_)

'''

