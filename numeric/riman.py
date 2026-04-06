import math
def int_rim(a, b, c, n, m, n_alpha, m_beta, q = lambda x,y:1):
    ax = a[0]
    ay = a[1]
    bx = b[0]
    by = b[1]
    cx = c[0]
    cy = c[1]
    summ = 0
    dbeta = 1 / m_beta
    Ax = (ax - cx)
    Bx = (bx - cx)
    Ay = (ay - cy)
    By = (by - cy)
    #print(ay,by,cy)
    #print(Ay)
    Jacob = Ax * By - Ay * Bx
    for i in range(m_beta):
        beta = dbeta * (i + 0.5)
        dalpha = (1-beta)/n_alpha
        suma = 0
        #print(beta)
        for j in range(n_alpha):
            alpha = (j+0.5)*dalpha
            x = (Ax*alpha+Bx*beta+cx)
            y  = (Ay*alpha+By*beta+cy)
            suma = suma + x**n*y**m * q(x,y)
            
        '''
        while alpha<=beta:
            #print(alpha)
            suma = suma + (Ax*alpha+Bx*beta+cx)**n*(Ay*alpha+By*beta+cy)**m
            alpha = alpha+dalpha
        '''
        summ = summ+suma*dalpha
    #print(Jacob)
    return summ*dbeta*Jacob

ABC=[[-1,-31/2.0],[6,2],[4,4/19.0]]
n=4
m=24-14
an = 500*2
bm = 500*2
#ABC=[[0,1],[1,0],[0,0]]
#ABC = [[-3,2],[6, 2], [4,  4/3.0]]
#ABC = [[-3,3],[1,1],[-2,3]]
#ABC=[[0,1],[1,0],[0,0]]
ABC = [ [float(ABC[i][j]) for j in range(2)] for i in range(3)]   
i = int_rim(ABC[0],ABC[1],ABC[2],n,m,an,bm)
print(i)