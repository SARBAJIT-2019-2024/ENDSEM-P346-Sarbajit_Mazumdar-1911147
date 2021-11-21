#         #   #     #         #            #     # # # #        # # # #      #     #
# #     # #    #   #          #            #     #     #        #     #       #   #
#  #   #  #     # #           #            #     #     #        #     #        # #
#   # #   #      #            #            #     # # # #        # # # #         #
#    #    #      #            #            #     #     #        #               #
#         #      #            #            #     #     #        #               #  
#         #      #   #######  # # # # #    #     # # # #   #    #               # 

#####################################################
#Calling a matrix
''''''
#list_C=[]
#with open("txt file") as matC:
    #for k in matC:
        #list_C.append(list(map(float, k.split()))

       
                

## Custom library for importing functions
##################################################################################################
#Partial Pivot
def partial_pivoting(a,b,n):
    n=len(a)
    no_of_swaps=0
    for i in range(n-1):
        if abs(a[i][i]) == 0:
            for j in range(i+1,n):
                if abs(a[j][i]) > abs(a[i][i]):
                    a[j], a[i] = a[i], a[j]  # interchange ith and jth rows of matrix 'A'
                    b[j], b[i] = b[i], b[j]  # interchange ith and jth elements of vector 'b'
            no_of_swaps+=1
    return no_of_swaps        

###################################################################################################
#defining Gauss Jordan
def Gauss_jordan(list_C):
    a,b=Matrixmaker(list_C)
    n=len(a)
    
    partial_pivoting(a,b,n) #calling the partial pivot
    
    for i in range(n):
        pivot_element = a[i][i]
        
        
        if pivot_element==0:      #condition if pivot element turns out be zero in any step
            return (None,None)
        if type(b[i]) is list: 
            for l in range(n):
                if b[i][l] != 0:
                    b[i][l] = b[i][l]/pivot_element
            for r in range(i, n):
                a[i][r] = a[i][r]/pivot_element
            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j]  
                    for l in range(n):
                        if b[i][l] != 0:
                            b[k][l] = b[k][l] - balance_factor*b[i][l]

                        
        else: 
            b[i] = b[i]/pivot_element                  #Condition for pivot element
            for k in range(i, n):
                 a[i][k] = a[i][k]/pivot_element       #dividing the elements by pivot element    

            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    b[k] = b[k] - balance_factor*b[i]
                    for j in range(i, n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j] 
    return(a,b)
###############################################################################################
# Defining matrix maker which makes a matrix C in to A and b
    
def Matrixmaker(list_C):
    list_A=[[0 for x in range(len(list_C))] for y in range(len(list_C))]
    for i in range(len(list_C)):
            for j in range(len(list_C)):
                list_A[i][j]=list_C[i][j] 
    if len(list_C[0])==len(list_C)+1:          #Generating list of C  
        list_B=[0 for x in range(len(list_C))] #generating a list of B where we store the result
        for i in range(len(list_C)):
            list_B.append(0)
        
        for i in range(len(list_C)):
            list_B[i]=list_C[i][len(list_C)]
    else:
        list_B=[[0 for x in range(len(list_C))] for y in range(len(list_C))] 
        for i in range(len(list_C)):

            list_B[i][i]=1       
    return(list_A,list_B)

#####################################################################################################
 #Defining matrix multiplication

def matrix_mul(a,b):
    AB_=[[0 for x in range(len(a))] for y in range(len(a))] 
        
    for i in range(len(a)):
        for j in range(len(b[2])):#expanding along the coloumn 2 of B matrix
            for k in range(len(b)):
                AB_[i][j] += a[i][k]*b[k][j] 
    return AB_
#####################################################################################################    

# Defining function for determinant calculation

def determinant_calc(a):           
    if len(a) != len(a[1]):         
        print("The determinant is not defined, Provide  Square Matrix please !")#checking for entry of Square matrix.
    else:     
        no_of_swaps = 0
        for i in range(len(a)-1): 
            if abs(a[i][i])== 0:
                for j in range(i+1 , len(a)): 
                    if  abs(a[i][i]) < abs(a[j][i]): 
                        a[j], a[i] = a[i], a[j]        #swapping the rows of the matrix.
                no_of_swaps += 1  
        
        for i in range(len(a)):
            pivot_element = a[i][i]
            if pivot_element==0:
                return (0)
            for k in range(len(a)):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,len(a)):
                        a[k][j] = a[k][j] - balance_factor*(a[i][j]/pivot_element) #making the elements zero by balance factor
        det = 1
        for i in range(len(a)):
            det*=a[i][i]
            
        det*=(-1)**(no_of_swaps)    
        #print("The determinant for the given matrix is:")        
            
        return(det)     
                
#################################################################################################################################
#################################################################################################################################







#CROUT METHOD

def crout_method(a):
    n = len(a)
    for i in range(n):
        for j in range(i,n):
            sum_lower = 0
            for k in range(i):
                sum_lower += (a[i][k]*a[k][j])
            a[i][j] = a[i][j] - sum_lower
        for j in range(i,n):
            if (i == j):
                pass
            else:
                sum_upper = 0
                for k in range(i):
                    sum_upper += (a[j][k]*a[k][i])
                a[j][i] = (a[j][i]- sum_upper)/a[i][i]
    Lower_tri, Upper_tri = [[0 for x in range(n)] for y in range(n)], [[0 for x in range(n)] for y in range(n)] #Breaking the matrix into two parts upper and lower
    partial_pivoting(Upper_tri,Lower_tri,n)
    for i in range(n):
        for j in range(i+1):
            if (i == j):     # conditioning for crout U[i][i]=1
                Upper_tri[j][j] = 1
                Lower_tri[j][j] = a[j][j]
            else:
                Upper_tri[j][i], Lower_tri[i][j] = a[j][i], a[i][j]
    return  Upper_tri,Lower_tri
            
#####################################################################################################################################
#####################################################################################################################################


def combiner(Upper_tri,Lower_tri):  # combining upper and lower matrix in one matrix 
    n=len(Upper_tri)
    D=[[0 for x in range(n)] for y in range(n)]
    for i in range(len(Upper_tri)):
        for j in range(len(Upper_tri)):
            if i==j:
                D[j][j]=Lower_tri[i][j]   #replacing diagonal elements by Lower triangular's diagonal elements
            else:
                D[i][j]=Lower_tri[i][j]+Upper_tri[i][j]
    return D            


   
####################################################################################################################################
####################################################################################################################################

#FORWARD BACKWARD SUB FOR DO LITTLE

def forward_backward_sub_do_little(Upper_tri, Lower_tri, b):
    n=len(Upper_tri)
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += Lower_tri[i][j] * y[j]
        y[i] = b[i] - total

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += Upper_tri[i][j] * x[j]
        
        x[i] = (y[i] - total)/Upper_tri[i][i]

    return x
########################################################################################################################
########################################################################################################################

#FORWARD BACKWARD SUB FOR CROUT

def forward_backward_sub_crout(Upper_tri, Lower_tri, b):
    n=len(Upper_tri)
    partial_pivoting(Upper_tri,Lower_tri,n)
    D= combiner(Upper_tri,Lower_tri)
    
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += D[i][j] * y[j]
        y[i] = b[i] - total

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += D[i][j] * x[j]
        
        x[i] = (y[i] - total)/D[i][i]

    return x


  
#####################################################################################################################
#####################################################################################################################

#CROUT METHOD FOR solution
def linear_solver_crout(list_C):
    A,b=Matrixmaker(list_C)
    n= len(A)
    
    partial_pivoting(A, b, n)
    for i in range(len(A)):
        if A[i][i]==0:
            print("No solutions exist for this system !")
            return None
        
    Upper_tri, Lower_tri = crout_method(A)
    x = forward_backward_sub_crout(Upper_tri,Lower_tri, b)
    print("The solutions of the system of linear equations by Crout's method is")
    return x

################################################################################################################
################################################################################################################

# UPPER LOWER BY DO-LITTLE
def Do_little(a):
    
    n= len(a)
    Lower_tri = [[0 for x in range(n)]
             for y in range(n)]
    Upper_tri = [[0 for x in range(n)]
             for y in range(n)]
 
    for i in range(n):
        Lower_tri[i][i] = 1

    for j in range(n):
        for i in range(n):
            total = 0
            for k in range(i):
                total +=Lower_tri[i][k] * Upper_tri[k][j]

            if i == j:
                Upper_tri[i][j] = a[i][j] - total

            elif i > j:
                Lower_tri[i][j] = (a[i][j] - total)/Upper_tri[j][j]

            else :
                Upper_tri[i][j] = a[i][j] - total

    return Upper_tri, Lower_tri
################################################################################################################
################################################################################################################
#Dolittle METHOD FOR solution
def linear_solver_do_little(list_C):
    A,b=Matrixmaker(list_C)
    n= len(A)
    
    partial_pivoting(A, b, n)
    for i in range(len(A)):
        if A[i][i]==0:
            print("No solutions exist!")
            return None
    Upper_tri, Lower_tri = Do_little(A)
    x = forward_backward_sub_do_little(Upper_tri,Lower_tri, b)
    print("The solutions of the system of linear equations by Dolittle's method is")
    return x

##################################################################################################################
##################################################################################################################

#printing upper and lower triangular matrix

def print_matrix(Lower_tri,Upper_tri):
    
    for i in range(len(Lower_tri)):
        if Lower_tri[i][i]==1:
            for j in range(len(Lower_tri)):
                if i>=j:
                    print(Upper_tri[i][j],end=" ")
                else:
                    print(Lower_tri[i][j],end=" ")  
            print()          
        
        else :
            for j in range(len(Lower_tri)):
                if i>j:
                    print(Upper_tri[i][j],end=" ")
                else:
                    print(Lower_tri[i][j],end=" ")  
            print() 
    return  None

#####################################################################################################################
#####################################################################################################################

def Nullmatrix(m,n): # creating a null matrix
	p= [[0 for i in range(n)] for j in range(m)]
	return(p)
#Creating square matrix
def identity_mat(m):
	p=Nullmatrix(m,m)
	for i in range(m):
		p[i][i] = 1
	return(p)

###################################################################################################################
###################################################################################################################
def pivoting_matrix(a): #PA=LU
    
    m = len(a)

    # Creating an identity matrix                                                                                                                                                                                          
    identity_matrix = [[float(i ==j) for i in range(m)] for j in range(m)]

    #The largest element in each row will be placed into the diagonal                                                                                                                                                                                             
    for j in range(m):
        row = max(range(j, m), key=lambda i: abs(a[i][j]))
        if j != row:
            # Swapping the rows                                                                                                                                                                                                                            
            identity_matrix[j], identity_matrix[row] = identity_matrix[row], identity_matrix[j]

    return identity_matrix
################################################################################################################################
################################################################################################################################

#inverse using LU inverse
def LU_inverse(a):
    
    
    
    I=identity_mat(4)
    P = pivoting_matrix(a) #calling the pivoting matrix
    PA = matrix_mul(P, a) #PA=LU

    P_1=Nullmatrix(4,4) # as P_inverse= P_transpose 
    for i in range(4):
	    for j in range(4):
		    P_1[j][i]=round(P[i][j],3)


    Upper_tri,Lower_tri=crout_method(PA)
    
    det=1            # condition for inverse existing
    for i in range(len(Lower_tri)):
        det*=Lower_tri[i][i]
    if det==0:
        print("Inverse doesn't exist!")
        return None     
# storing matrix
    x_1=Nullmatrix(4,4)
    x_2=Nullmatrix(4,4)  


    for i in range(4):
	    x_1[i]=forward_backward_sub_crout(Upper_tri,Lower_tri,I[i])

    Px_i = matrix_mul(P_1,x_1)	


    for i in range(4):
	    for j in range(4):
		    x_2[j][i]=round(Px_i[i][j],5)

    return x_2



################################################################################################################
################################################################################################################    

#Forward backward method for cholesky's substitution 


def forward_backward_sub_cholesky(Lower_tri, Lower_tri_Transpose, b):
    n=len(Lower_tri)
    
    
    
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += Lower_tri[i][j] * y[j]
           
        y[i] = (b[i] - total)/Lower_tri[i][i]

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += Lower_tri_Transpose[i][j] * x[j]
        
        x[i] = (y[i] - total)/Lower_tri_Transpose[i][i]

    return x


#################################################################################################################
#################################################################################################################

#CHOLESKY's METHOD

def cholesky_decomp (C):
    a,b=Matrixmaker(C)
    n=len(a)
    Lower_tri = [[0 for x in range(n)]
             for y in range(n)]
    Lower_tri_Transpose = [[0 for x in range(n)]        #Generating two matrices L and L^T(Transpose) such that A=L*L^T
             for y in range(n)]         
                                                   

    for i in range(n):
        for j in range(i + 1):      #Fixing range for the iterations
            total = 0
 
            if j == i:                           #confition for diagonal entries
                for k in range(j):
                    total +=Lower_tri[i][k] **2
            
                Lower_tri[j][j] =(a[i][i] - total)**(0.5)   #taking the sqaure root
                
            else:
                for k in range(j):
                    total += (Lower_tri[i][k] *Lower_tri[j][k])
                if(Lower_tri[j][j] > 0):
                    Lower_tri[i][j] = ((a[i][j] - total) /
                                               Lower_tri[j][j])
    for i in range(len(Lower_tri)):
        for j in range(len(Lower_tri)):                             #making the L^T from L
                Lower_tri_Transpose[i][j]=Lower_tri[j][i]
    
    return Lower_tri,Lower_tri_Transpose

def Cholesky_Solver(A):    #making a function cholesky solver to solve equation by cholesky's method
    a,b=Matrixmaker(A)
    Lower_tri,Lower_tri_Transpose=cholesky_decomp(a)
    x = forward_backward_sub_cholesky(Lower_tri,Lower_tri_Transpose, b)
    print("The solutions of the system of linear equations by Cholesky's method is")
    return x

##################################################################################################################
##################################################################################################################
###################################################################################################################



#Function for Bracketing method
def bracketing(f, a_i, b_i): # defining bracketing method
    custom_iteration=0 # initialise an iteration
    beta = 0.05 # interval step of each bracket
    if a_i > b_i: # criteria for wrong ordering
        print("Please enter a value such that a < b")
        return None,None
    else:
        
        if f(a_i)*f(b_i) < 0: # correct initilisation
            return a_i, b_i
        else:
            while f(a_i)*f(b_i) > 0: 
                if abs(f(a_i)) < abs(f(b_i)):
                    a_i = a_i - beta*(b_i-a_i) # bracketing the interval
                    custom_iteration += 1 # adding up the iteration
                    prod = f(a_i) * f(b_i)

                elif abs(f(a_i)) > abs(f(b_i)):
                    b_i = b_i + beta*(b_i-a_i) # bracketing the interval
                    custom_iteration += 1
                    prod = f(a_i) * f(b_i)
                if custom_iteration > 30: # No of iterations after which it will stop searching for another bracket numerically
                    print("Try another bracket range")
                    return None,None
                
            
            return a_i, b_i

# Bisection method of finding roots of a given equation
def bisection_method(f, a_i, b_i, precision): # defining bisection function
    a_i, b_i = bracketing(f, a_i, b_i) #Calling the bracketing function for correcting the bracket
    if a_i==None and b_i==None:
        return None,None,None,None
    No_of_iterations = [] #blank list for iteration 
    root_val_i = [] #blank list for rootvalue c_i 
    x_i = [] # blank list for error in root
    f_x_i=[] # blank list for functional value after each iteration
    i=0  #No of iterations
    while abs(a_i-b_i)>precision: #precision condition
        
        c_i = (a_i+b_i)/2  # intial c
        
        if f(a_i) * f(c_i) < 0: #Intermediate value property
            b_i = c_i
        elif f(a_i) * f(c_i)  > 0:
            a_i = c_i

        No_of_iterations.append(i) # adding iteration in the blank list
        root_val_i.append(c_i) # adding value in root values blank list
        f_x_i.append(f(c_i)) # Adding f(x_i) value after each iteration
        x_error = abs(root_val_i[i]-root_val_i[i-1]) # defining error value
        
        x_i.append(x_error)# adding value in blank list of root error 
        i+=1
        
    else:
        return c_i, x_i,f_x_i, No_of_iterations
#########################################################################################################
#########################################################################################################
#Code for Regula Falsi method
def regula_falsi_method(f, a_i, b_i, precision): # defining regula falsi function
    a_i, b_i = bracketing(f, a_i, b_i) # calling the bracketing function
    No_of_iterations = []#blank list for iteration 
    root_val_i = []#blank list for rootvalue c_i 
    x_i = []# blank list for error in root
    f_x_i=[]# blank list for functional value after each iteration
    i=0
    c_i = a_i  # intial c
    while abs(a_i-b_i) >= precision:
        
        c_i = b_i - ((b_i-a_i)*f(b_i))/(f(b_i) - f(a_i)) #using general formula
        if f(c_i)==0:
            return c_i, x_i,f_x_i, No_of_iterations
        if f(a_i) * f(c_i) < 0:    # intermediate value property
            b_i = c_i
        elif f(a_i) * f(c_i) > 0:
            a_i = c_i
        
        
        root_val_i.append(c_i)# adding value in root values blank list
        x_error = abs(root_val_i[i] - root_val_i[i-1])# defining error value
        No_of_iterations.append(i)# adding iteration in the blank list
        x_i.append(x_error)# adding value in blank list of root error 
        f_x_i.append(f(c_i))# Adding f(x_i) value after each iteration
        i+=1
    else:
        return c_i, x_i,f_x_i, No_of_iterations

##############################################################################################
##############################################################################################
# Differentiation function
def differentiation(f,x_0):
    
    h=10**(-10)
    f1_x_0 = (f(x_0+h) - f(x_0-h))/(2*h) #using general formula
    return f1_x_0


#function for Newton Rhapson Method
def Newton_Raphson_Method(f, x_0, precision): #defining Newton Rhapson Method
    No_of_iterations = [] # Making a blank list for no of iterations
    f_x_i=[] # making a blank list for functional value after each iteration
    i=0
    if f(x_0)==0:
        return x_0,f_x_i,No_of_iterations

    while abs(f(x_0)/differentiation(f, x_0))>precision: #conition for precision value
        
        x_0 = x_0 - (f(x_0)/differentiation(f, x_0)) # using general formula for newton Rhapson
        No_of_iterations.append(i) # adding iteration to the iteration blank list
        f_x_i.append(f(x_0)) # adding function value into the functional value  blank list
        i+=1 # adding iterations

    else:
        return x_0,f_x_i,No_of_iterations


def polynomial_function(a_i, x): #defining a polynomial function
    
    P_x = 0
    for i in range(len(a_i)):
       P_x=P_x+ a_i[i]*(x**(len(a_i)-1-i))

    return P_x

def polynomial_df(a_i,x_0): # defing 1st derivative at x_0 
    h=10**(-4) #precision
    p1_x_0 = (polynomial_function(a_i,x_0+h) - polynomial_function(a_i,x_0-h))/(2*h)
    return p1_x_0
def polynomial_ddf(a_i,x_0):# defing 2nd derivative at x_0 
    h=10**(-5)
    p2_x_0 = (polynomial_df(a_i,x_0+h)-polynomial_df(a_i,x_0-h))/(2*h)
    
    return p2_x_0  



##################################################################################################
##################################################################################################    

# Function for lauguerre method

def laguerre(a_i, alpha, precision): #defining the lauguerre function
    n = len(a_i) - 1    #fixing the degree of the polynomial
    
   
    while abs(polynomial_function(a_i, alpha)) == 0: # if alpha=0 is a root it will give alpha as an output.
        return alpha

    else:  #else part
        
        alpha_next=10**(10)  #Fixing a value such that |alpha_next-alpha| never exceeds the precision value.
        while  abs(alpha_next-alpha)>precision:  
            
            if alpha_next!=10**10:
                alpha = alpha_next #changing alpha value after each iteration
            G = polynomial_df(a_i,alpha)/(polynomial_function(a_i, alpha)) # Defining G
            H = G**2 - (polynomial_ddf(a_i, alpha)/polynomial_function(a_i, alpha)) # Defining H
            
            if abs((G + ((n-1)*(n*H - G**2))**0.5))>abs((G - ((n-1)*(n*H - G**2))**0.5)): #Choosing the sign in the denominator to fix the larger absolute value. 
                a = n/(G + ((n-1)*(n*H - G**2))**0.5)
            else:
                a = n/(G - ((n-1)*(n*H - G**2))**0.5)
            
            alpha_next = alpha -a
            
        else:
            return alpha

def syn_div(a_i,alpha): # defining synthetic division
    for i in range(len(a_i)-1):
        a_i[i+1]=a_i[i+1]+alpha*a_i[i] # dividing original polynomial and changing the co efficients.
    a_i.pop()
    return a_i


def driver_function_lag(a_i,alpha,precision): #Defining driver function for lauguerre as the upper function will only give one root per time.
    rootlist=[] # creating an empty rootlist
    for i in range(len(a_i)-1):
        rootlist.append(laguerre(a_i, alpha, precision))
        a_i=syn_div(a_i,rootlist[i])  #adding roots to the blank list

    return rootlist    

###########################################################################################################################################3333
#############################################################################################################################################33

# Integration using the midpoint method
def integration_midpoint(f, upper_limit,lower_limit, N):
    
    sum = (upper_limit - lower_limit)/N #defining h
    total = 0
    for p in range(1,N+1):
        x = (p*sum +(p-1)*sum + 2*lower_limit)/2
        total += sum*f(x)

    return total


# Integration using the trapezoidal rule
def integration_trap(f, upper_limit, lower_limit, N):
   

    sum = (upper_limit - lower_limit)/N #defining h
    total = 0
    for i in range(1, N+1):
        x_n = lower_limit + i*sum # x ~ x+ih
        x_0 = x_n - sum 
        T_n = (f(x_0) + f(x_n))*(sum/2) # defining T_n as h/2*(f(x_n-1)+f(x_n))
        total += T_n

    return total    

# Integration via Simpson method
# Direct Functional approach

#def integration_simp(f, upper_limit, lower_limit, N):
    

    #h = (upper_limit - lower_limit)/N
    #total = 0
    #for i in range(1, N+1, 2):   
        #x_n_2 = lower_limit + i*h
        #x_n = x_n_2 - h
        #x_n_1 = (x_n + x_n_2)/2
        #S_n = (f(x_n) + 4*f(x_n_1) + f(x_n_2))*(h/3)
        #total += S_n

    #return total    


def integral_simp(f,b,a,N):                                         
    h, sum_val = (b-a)/N, 0              #defining h                               
    for i in range(N+1):
        x_i=a+i*h          #x ~ x+ih
        if i == 0 or i == N:
            sum_val += f(x_i)                                       # w = 1 
        elif i % 2 != 0:
            sum_val += 4*f(x_i)                                   # w = 4 for i=2k+1
        else:
            sum_val += 2*f(x_i)                                     # w = 2 for i=2k
    return sum_val*(h/3)





# monte carlo integration method
import random
def monte_carlo(f, lower_limit,upper_limit, N):
    # Generate list N set of randomly generated numbers
    Random_list = [] # creating a random list
    for i in range(N):
        Random_list.append(random.uniform(lower_limit, upper_limit)) #generating system build random data points

    sum = 0
    for i in range(N):
        sum += f(Random_list[i])

    total = (upper_limit-lower_limit)/N * sum

    return total, Random_list








    ##############################################################################################################################
    ##################################################################################################################################
    #############################################ODE SOLVE#########################################################################################

def euler_b(func,h,x_0,y_0,x_n,guess = 1):
        X=[x_0]
        Y=[y_0]
        N=int(abs(x_0-x_n)/h)                              
        for i in range(N):
            x=X[-1]+h
            f_raph = lambda y_raphson: Y[-1] + h*func(x,y_raphson) - y_raphson
            y_nr,b,c=Newton_Raphson_Method(f_raph,guess,10**-5)
            y=Y[-1]+h*func(x,y_nr)                       
            X.append(x)
            Y.append(y)
        return X,Y





def f_euler(df, y_0, x_0, x_n, h):
    
    x = [x_0] #listing
    y = [y_0]
         
    for i in range(int((x_n-x_0)/h)): #iteration upto  no. of steps
        x.append(x[i] + h)

    for i in range(int((x_n-x_0)/h)):
        y.append(y[i] + h* df(y[i], x[i]))

    return x, y

def predictor_corrector(df,y_0,x_0,x_n,h):
    X = [x_0]
    Y = [y_0]
    for i in range(int((x_n-x_0)/h)):
        k_1 = h*df(Y[i],X[i])
        yp = Y[i]+k_1
        x = X[i]+h
        k_2 = h*df(yp,x)
        y_c = Y[i] + (k_1+k_2)/2
        X.append(x)
        Y.append(y_c)
    return(X,Y)
   
#RK4

def rk_4(df_2, df, x0, y0, z0, x_N, h): # for NEUMANN BOUNDARY CONDITION
    
    x = []
    x.append(x0)
    y = []
    y.append(y0)
    z = []
    z.append(z0)
    
    
    

    
    for i in range(int((x_N-x0)/h)): #range up to the iterations will go
        x.append(x[i] + h)
        k1 = h * df(x[i], y[i], z[i])
        k1v = h * df_2(x[i], y[i], z[i])
        k2 = h * df(x[i] + h/2, y[i] + k1/2, z[i] + k1v/2)
        k2v = h * df_2(x[i] + h/2, y[i] + k1/2, z[i] + k1v/2)
        k3 = h * df(x[i] + h/2, y[i] + k2/2, z[i] + k2v/2)
        k3v = h * df_2(x[i] + h/2, y[i] + k2/2, z[i] + k2v/2)
        k4 = h * df(x[i] + h, y[i] + k3, z[i] + k3v)
        k4v = h * df_2(x[i] + h, y[i] + k3, z[i] + k3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (k1v + 2*k2v + 2*k3v + k4v)/6)

    return x, y, z


#LAGRANGE INTERPOLATION
def lagrange_interpolation(z_h, z_l, y_z_h, y_z_l, y):
    zeta = z_l + (z_h - z_l) * (y - y_z_l)/(y_z_h - y_z_l)
    return zeta

#Shooting Method
# Solves 2nd order ODE given Dirichlet boundary conditions


#def lagrange_interpolator(c_h,c_l,y_ich,y_icl,y_i):
  #  return c_l + ((c_h - c_l)/(y_ich - y_icl))*(y_i - y_icl)

def shooting_method_2(df,df_2,x_a,x_b,y_a,y_b,h,v_l,v_h):
   
    for i in range(int((x_b-x_a)/h)): 
        x_l,y_l =  rk_4(df,df_2,x_a,y_a,v_l,h,x_b)
        x_h,y_h =  rk_4(df,df_2,x_a,y_a,v_h,h,x_b)
        c_l=v_l
        c_h=v_h
        
        if abs(y_l[len(y_l)-1] - y_b) < 0.001:
            return x_l,y_l
        if abs(y_h[len(y_h)-1] - y_b) < 0.001:
            return x_h,y_h
        if y_l[len(y_l)-1] < y_b and y_h[len(y_h)-1] > y_b:
            c = lagrange_interpolation(c_h,c_l,y_h[len(y_h)-1],y_l[len(y_l)-1],y_b)
            
            x_c,y_c = rk_4(df,df_2,x_a,y_a,c,h,y_b)
            if abs(y_c[len(y_c)-1] - y_b) < 0.001:
                return x_c,y_c
            elif c < y_b:
                c_l = c
            else:
                c_h = c
        #if y_l[len(y_l)-1] > y_b and y_h[len(y_h)-1] < y_b:
            c = lagrange_interpolation(c_l,c_h,y_h[len(y_h)-1],y_l[len(y_l)-1],y_b)
            x_c,y_c = rk_4(df,df_2,x_a,y_a,c,h,y_b)
            if abs(y_c[len(y_c)-1] - y_b) < 0.001:
                return x_c,y_c
            elif c < y_b:
                c_h = c
            else:
                c_l = c
                                




def shooting_method(df_2, df, x0, y0, x_n, y_n, z_1, z_2,h, resilience=10**(-5)):
   
    x, y, z = rk_4(df_2, df, x0, y0, z_1, x_n, h)
    yn = y[-1]

    if abs(yn - y_n) > resilience:
        if yn < y_n:
            zeta_l = z_1
            yl = yn

            x, y, z = rk_4(df_2, df, x0, y0, z_2, x_n,h)
            yn = y[-1]

            if yn > y_n:
                zeta_h = z_2
                yh = yn

                
                zeta = lagrange_interpolation(zeta_h, zeta_l, yh, yl, y_n) #interpolation

               
                x, y, z = rk_4(df_2, df, x0, y0, zeta, x_n, h)
                return x, y, z

            else:
                print("CHOOSE ANOTHER VALUE")


        elif yn > y_n:
            zeta_h = z_1
            yh = yn

            x, y, z = rk_4(df_2, df, x0, y0, z_2, x_n, h)
            yn = y[-1]

            if yn < y_n:
                zeta_l = z_2
                yl = yn

                
                zeta = lagrange_interpolation(zeta_h, zeta_l, yh, yl, y_n) #interpolation
                x, y, z = rk_4(df_2, df, x0, y0, zeta, x_n, h)
                return x, y, z

            else:
                print("CHOOSE ANOTHER VALUE")


    else:
        return x, y, z         




################################################################################################################################
#################################################################################################################################
#CURVE FITTING

#Least square fitting
#Least Square Fit # straight line
def LeastSquare(X,Y):
    n = len(X)
    x_1=0;y_1=0;x_2=0;xy=0;y_2=0
    for i in range(n):
        x_1+=X[i]
        y_1+=Y[i]
        x_2+=pow(X[i],2)
        y_2 += pow(Y[i], 2)
        xy+=X[i]*Y[i]
    x1=x_1/n
    y1=y_1/n
    a = ((y1*x_2)-(x1*xy))/(x_2-(n*pow(x1,2)))
    b = (xy - n*x1*y1)/(x_2-n*pow(x1,2))
    Sxx = x_2-n*pow(x1,2)
    
    Syy = y_2-n*pow(y1,2)
    
    Sxy = xy-n*x1*y1
    sigma_x = Sxx/n
    sigma_y = Syy/n
    cov_xy = Sxy/n
    r_xy = pow(Sxy,2)/(Sxx*Syy)

    if b<0:
        r_xy=-r_xy**(0.5)
        
    else:
        r_xy=r_xy**(0.5)
    return a,b,sigma_x,sigma_y,cov_xy,r_xy




      
#Polynomial Fitting
def polynomial_curve_fitting(X,Y,n):
    if len(X) == len(Y):
        N, n = len(X), n + 1
        Mat_X, Mat_Y = [[0 for i in range(n)] for j in range(n)], [0 for i in range(n)] #list of list 
        for i in range(n):
            for j in range(n):
                valX = 0
                for k in range(N):
                    valX += X[k]**(i+j)
                Mat_X[i][j] = valX
        for i in range(n):
            valY = 0
            for j in range(N):
                valY += (X[j]**(i))*(Y[j])
            Mat_Y[i] = valY
        L,U=crout_method(Mat_X)
        print("The matrix required for polynomial fit ",Mat_X)
    p=determinant_calc(Mat_X)
    if p!=0:
        print("Inverse Exists and determinant is ",p)
        Mat_result = forward_backward_sub_crout(L,U,Mat_Y) #calling crout for forward and backward substitutition
        return Mat_result
    else:
        print("provide same number data values for x and y !")
