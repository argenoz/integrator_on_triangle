import riman as RI
import sympy as SP
import integrator_ as INTI

'''
    dlya 


'''
#ABC [[ax,ay],[bx,by],[cx,cy]]
#Fi = [(dfi_i,fi_i)], fi_i lambda x,y dlya bazisnyh, a dfi_i eto cherez biharminicheskyi
#NM - ochevidno
#q - shtuchka sprava ot znaka ravenstyva
def Galerkin_(ABC,Fi,q=lambda x,y:0,NM=[1000,1000]):
    l = range(len(Fi))
    ans = [[0 for j in l] for i in l]
    #ans_ij eto proizvedenie dfi_j na fi_i
    def Riman_dlya_Galerkina(ABC,f,NM):
        summ = 0
        cx = ABC[2][0]
        cy = ABC[2][1]
        Ax = ABC[0][0]-cx
        Bx = ABC[1][0]-cx
        Ay = ABC[0][1]-cy
        By = ABC[1][1]-cy
        
        Jacob = Ax * By - Ay * Bx
        dbeta = 1.0/NM[0]
        ml = range(NM[1])
        for i in range(NM[0]):
            suma = 0
            beta = dbeta*(i+0.5)
            dalpha=(1-beta)/NM[1]
            for j in range(NM[1]):
                alpha = dalpha*(0.5+j)
                #x = (Ax*alpha+Bx*beta+cx)
                #y  = (Ay*alpha+By*beta+cy)
                #suma = suma+f(x,y)
                suma = suma+f(Ax*alpha+Bx*beta+cx,Ay*alpha+By*beta+cy)
            summ = summ+suma*dalpha
        return Jacob*dbeta*summ
    ml = None
    for i in l:
        for j in l:
            ml = lambda x,y:1
class Galerkin:
    '''
    Fiq - pervyi element eto spisok par: [(dfi,fi)], gde dfi eto lambda u biharmonimi i fi, a fi est' fi
            vtoroi element eto q(x,y)
    
    '''
    def __init__(self,ABC,Fiq,NM=[1000,1000]):
        self.ABC = ABC
        self.cx = ABC[2][0]
        self.cy = ABC[2][1]
        self.Ax = ABC[0][0]-self.cx
        self.Bx = ABC[1][0]-self.cx
        self.Ay = ABC[0][1]-self.cy
        self.By = ABC[1][1]-self.cy
        self.Fi=Fiq[0]#eto nabor par
        self.q=Fiq[1]#eto odna funkciya q(x,y)
        self.Jacob = self.Ax * self.By - self.Ay * self.Bx
        self.NM = NM
        self.dbeta = 1/NM[0]
        self.dalphas=[ (1-(0.5+i)*self.dbeta)/NM[1]   for i in range(NM[0])]
        self.lNM = [range(self.NM[0]),range(self.NM[1])]
    def do_analitic(self,f):
        return INTI.Integrator(self.ABC).integrate(f[0],f[1])
    def do_riman(self,f):
        summ=0
        for i in self.lNM[0]:
            suma = 0
            beta = self.dbeta*(i+0.5)
            dalpha = self.dalphas[i]
            
            for j in self.lNM[1]:
                alpha = dalpha*(0.5+j)
                x = self.Ax*alpha+self.Bx*beta+self.cx
                y = self.Ay*alpha+self.By*beta+self.cy
                suma = suma+f(x,y)
            summ = summ+suma*dalpha
        
        return summ*self.Jacob*self.dbeta
        
    def do(self):
        l = range(len(self.Fi))
        ans = [[0 for j in l] for i in l]
        qq=[]
        for i in l:#po kazdoi stroke.
            for j in l:#po kazhdomu stolbcy
                f = lambda x,y: self.Fi[i][1](x,y)*self.Fi[j][0](x,y)
                ans[i][j]=self.do_riman(f)
                #print([i,j])
            f = lambda x,y: self.Fi[i][1](x,y)*self.q(x,y)
            qq.append([self.do_riman(f)])
        return [ans,qq]

x = SP.Symbol('x')
y = SP.Symbol('y')

n0,m0=3,3
n = 4
m = 3

Fi = []

for i in range(n):
    for j in range(m):
        fi = x**(i+n0)*y**(j+m0)
        dfi = SP.diff(fi,x,x,x,x)+SP.diff(fi,x,x,y,y)*2+SP.diff(fi,y,y,y,y)
        print(fi)
        dfi = (lambda dfi,ex,ey:lambda xi,yi:dfi.subs(ex,xi).subs(ey,yi))(dfi,x,y)
        fi = (lambda fi,ex,ey:lambda xi,yi:fi.subs(ex,xi).subs(ey,yi))(fi,x,y)
        Fi.append((dfi,fi))
print(len(Fi))
def getXY(n,m,n0=0,m0=0):
    ans = []
    Z = lambda x,y:0
    for i in range(n+1):
        for j in range(m+1):
            fi = lambda x,y: x**(i+n0)*y**(j+m0)
            p = fi
            if i+n0>=4:
                if j+m0>=4:
                    fi =lambda x,y: (i+n0)*(i+n0-1)*(i+n0-2)*(i+n0-3)*x**(i+n0-4)+2*(i+n0)*(i+n0-1)*x**(i+n0-2)\
                    *(j+m0)*(j+m0-1)*y**(j+m0-2)+(j+m0)*(j+m0-1)*(j+m0-2)*(j+m0-3)*y**(j+m0-4)
                elif j+m0>=2:
                    fi =lambda x,y: (i+n0)*(i+n0-1)*(i+n0-2)*(i+n0-3)*x**(i+n0-4)+2*(i+n0)*(i+n0-1)*x**(i+n0-2)\
                    *(j+m0)*(j+m0-1)*y**(j+m0-2)
                else:
                    fi =lambda x,y: (i+n0)*(i+n0-1)*(i+n0-2)*(i+n0-3)*x**(i+n0-4)
            elif i+n0>=2:
                if j+m0>=4:
                    fi =lambda x,y: 2*(i+n0)*(i+n0-1)*x**(i+n0-2)\
                    *(j+m0)*(j+m0-1)*y**(j+m0-2)+(j+m0)*(j+m0-1)*(j+m0-2)*(j+m0-3)*y**(j+m0-4)
                elif j+m0>=2:
                    fi =lambda x,y: 2*(i+n0)*(i+n0-1)*x**(i+n0-2)\
                    *(j+m0)*(j+m0-1)*y**(j+m0-2)
                else:
                    fi =Z
            else:
                if j+m0>=4:
                    fi =lambda x,y: (j+m0)*(j+m0-1)*(j+m0-2)*(j+m0-3)*y**(j+m0-4)
                elif j+m0>=2:
                    fi =Z
                else:
                    fi =Z
            ans.append([fi,p])
    return ans            

ABC=[[2,1],[1,2],[0,0]]
E = lambda x,y:1
#Fi=[[E,lambda x,y:x**2],[E,lambda x,y:y*x**2],[E,lambda x,y:x**2*y**3]]
q = lambda x,y:1
Fiq=[getXY(5,5),q]
#print(Fiq)
galer = Galerkin(ABC,Fiq,NM=[50,50])
galer = galer.do()[0]
print(galer)

