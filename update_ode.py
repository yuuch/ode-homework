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

def plot_lines(y1, y2,index, fig_name='myfig.svg',n=10000): # and computer the norm of error
    #matplotlib.use('agg')
    s =plt.figure(index)
    x=[100.0/n*i for i in range(n+1)]
    assert len(y1)==len(x)
    norm = max(abs(np.array(y1)-np.array(y2)))
    norm = float(norm)
    plt.subplot(1,3,1)
    plt.title('exact solution')
    plt.plot(x,y1,label='exact solution')
    plt.subplot(1,3,2)
    plt.title('numerical solution')
    plt.plot(x,y2,'r',label='numerical solution')
    plt.subplot(1,3,3)
    plt.plot(x,y1,label='exact')
    plt.plot(x,y2,'r--',label='numerical')
    plt.title(fig_name[:-4])
    plt.legend()
    plt.xlabel('the norm of the error is :%.5f'%norm)
    plt.savefig(fig_name,format='svg')
    return norm
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
def implicit_midpoint_method(u_0=1,v_0=1,n=10000):
    t_interval=[0,100]
    mat_of_coe=[[98,198],[-99,-199]]# matrix of coefficients
    u_value_list=[u_0]
    v_value_list=[v_0]
    h = (t_interval[1]-t_interval[0])*1.0/n
    #dace = 25*h**2+1+101*h/2
    #dace = (1-98*h/2)*(1+199*h/2)-(-99*h/2)*(198*h/2) # it is just a constant
    epsilon = 0.0001

    for i in range(n):
        u_s =u_value_list[i]
        v_s = v_value_list[i]
        while True:
            u_s_new =98*h/2*(u_value_list[i]+u_s)+198/2*h*(v_value_list[i]+v_s)
            v_s_new =-99/2*h*(u_value_list[i]+u_s)-199/2*h*(v_value_list[i]+v_s)
            if abs(u_s_new-u_s)<epsilon and abs(v_s_new-v_s)<epsilon:
                u_value_list.append(u_s_new)
                v_value_list.append(v_s_new)
                break
            else:
                u_s = u_s_new
                v_s = v_s_new
        #u_new = ((1+199.0*h/2)*(1+98.0*h/2)+198.0*h/2*(-99*h/2))*u_value_list[i]
        #u_new = u_new + 198*v_value_list[i]
        #u_new = (25*h**2+297*h/2+1)*u_value_list[i]+198*h*v_value_list[i]
        #u_new = u_new*1.0/dace
        #u_value_list.append(u_new)

        #v_new = (-99*h)*u_value_list[i]+((-99*h/2)*(99*h)+(1-49*h)*(1-199*h/2))*v_value_list[i]
        #v_new =-1*(99*h*u_value_list[i]+(25*h**2-1+197*h)*v_value_list[i])
        #v_new = v_new*1.0/dace
        #v_value_list.append(v_new)

    return u_value_list,v_value_list


    

if __name__=="__main__":
    n=10000
    exp_eul_u, exp_eul_v = explicit_euler_method(n=n)
    real_u, real_v = exact_values(n=n)
    imp_eul_u, imp_eul_v =implicit_euler_method(n=n)
    imp_mid_u, imp_mid_v = implicit_midpoint_method(n=n)
    #plot_lines(real_v,v_list,fig_name='u_lines.svg',n=n)
    #plot_lines(real_u,u_list,fig_name='v_lines.svg')
    exp_eul_u = plot_lines(real_u,exp_eul_u, index=1, fig_name='explicit_euler_u_lines.svg')
    exp_eul_v = plot_lines(real_v,exp_eul_v, index=2, fig_name='explicit_euler_v_lines.svg')
    imp_eul_u = plot_lines(real_u,imp_eul_u, index=3, fig_name='implicit_euler_u_lines.svg')
    imp_eul_v = plot_lines(real_v,imp_eul_v, index=4, fig_name='implicit_euler_v_lines.svg')
    imp_mid_u = plot_lines(real_u,imp_mid_u, index=5, fig_name='implicit_midpoint_u_lines.svg')
    imp_mid_v = plot_lines(real_v,imp_mid_v, index=6, fig_name='implicit_midpoint_v_lines.svg')

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


        
            





