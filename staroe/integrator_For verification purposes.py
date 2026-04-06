from copy import deepcopy as DPC
import sympy as SP


#       ax ay  bx by   cx  cy      
ABC = [[-1,2],[6, 2], [4,  4]]        
n = 2
m = 3
x = SP.Symbol('x')
y = SP.Symbol('y')




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


