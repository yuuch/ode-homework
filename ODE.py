import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
def explicit_euler_method(u_0=1 ,v_0=1,n=10000):
    # initial conditions
    t_interval=[0,100]
    mat_of_coe=[[98,198],[-99,-199]]# matrix of coefficients
    u_value_list=[u_0]
    v_value_list=[v_0]
    h = (t_interval[1]-t_interval[0])*1.0/n

    for i in range(n):
        #u_n+1 = (98*u_n+198*v_n) *h +u_n
        u_new = u_value_list[i] + (mat_of_coe[0][0]*u_value_list[i] + mat_of_coe[0][1] * \
                v_value_list[i] )*h
        u_value_list.append(u_new)
        v_new = v_value_list[i] + (mat_of_coe[1][0]*u_value_list[i] + mat_of_coe[1][1] * \
                v_value_list[i] )*h
        v_value_list.append(v_new)

    return u_value_list,v_value_list

def plot_lines(y1, y2, fig_name='myfig.svg',n=10000):
    #matplotlib.use('agg')
    x=[100.0/n*i for i in range(n+1)]
    assert len(y1)==len(x)
    plt.plot(x,y1,label='exact solution')
    plt.plot(x,y2,label='numerical solution')
    plt.title(fig_name[:-4])
    plt.legend()
    plt.savefig(fig_name,format='svg')
def exact_values(n=10000):
    t = np.array([100.0/n*i for i in range(n+1)])
    u_t = -3*np.exp(-100*t) + 4*np.exp(-1*t)
    v_t = 3*np.exp(-100*t) - 2*np.exp(-1*t)
    return u_t,v_t
def implicit_euler_method(u_0=1,v_0=1,n=10000):
    t_interval=[0,100]
    mat_of_coe=[[98,198],[-99,-199]]# matrix of coefficients
    u_value_list=[u_0]
    v_value_list=[v_0]
    h = (t_interval[1]-t_interval[0])*1.0/n

    for i in range(n):
        u_new = u_value_list[i]*(1+199*h)+198*h*v_value_list[i]
        u_new = u_new*(1.0/((1-98*h)*(1+199*h)+198*h*99*h))
        u_value_list.append(u_new)
        v_new = v_value_list[i]
        v_new = 99*h*u_value_list[i]-(1-98*h)*v_value_list[i]
        v_new = v_new*(1.0/((-198*h*99*h)-(1-98*h)*(1+199*h)))
        v_value_list.append(v_new)

    return u_value_list,v_value_list


    

if __name__=="__main__":
    n=10000
    v_list,u_list = explicit_euler_method(n=n)
    real_u ,real_v = exact_values(n=n)
    imp_eur_u,imp_eur_v =implicit_euler_method(n=n)
    #plot_lines(real_v,v_list,fig_name='u_lines.svg',n=n)
    #plot_lines(real_u,u_list,fig_name='v_lines.svg')
    plot_lines(real_u,imp_eur_u,fig_name='implicit_euler_u_lines.svg')

'''
class ODEs(object):
    def __init__(self,matrix,u_0,v_0,n):
        self.matrix = matrix
        self.exact_resolution = {'u':[u_0],'v':[v_0]}
        self.t_interval=[0,100]
        self.n = n 
        self.h = (self.t_interval[1]-self.t_interval[0])*1.0/n
        self.t = np.arange(n+1)*self.h
'''


        
            





