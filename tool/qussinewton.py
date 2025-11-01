import numpy as np
import numpy.linalg as linalg
## Check is element of array is NaN
def check_nan(A):
    return np.sum(np.isnan(A))

## Calculate gradient with central finite difference
def my_grad_approx(x,f,eps1,verbose=False):
    '''
    Calculate gradient of function my_f using central difference formula
    
    Inputs:
        x - point for which to evaluate gradient
        f - function to consider
        eps1 - perturbation size
        
    Outputs:
        grad - gradient (vector)
    '''
    
    n = len(x)
    grad = np.zeros(n)
    
    if(verbose):
        print("***** my_grad_approx at x = ",x,"*****")
        
    for i in range(0,n):
        
        # Create vector of zeros except eps in position i
        e = np.zeros(n)
        e[i] = eps1
        
        # Finite difference formula
        my_f_plus = f(x + e)
        my_f_minus = f(x - e)
        
        # Diagnostics
        if(verbose):
            print("e[",i,"] = ",e)
            print("f(x + e[",i,"]) = ",my_f_plus)
            print("f(x - e[",i,"]) = ",my_f_minus)
        
        
        grad[i] = (my_f_plus - my_f_minus)/(2*eps1)
    
    if(verbose):
        print("***** Done. ***** \n")
    
    return grad

## Calculate Hessian using cental finite difference
def my_hes_approx(x,grad,eps2):
    '''
    Calculate gradient of function my_f using central difference formula and my_grad
    
    Inputs:
        x - point for which to evaluate gradient
        grad - function to calculate the gradient
        eps2 - perturbation size (for Hessian NOT gradient approximation)
        
    Outputs:
        H - Hessian (matrix)
    '''
    
    n = len(x)
    H = np.zeros([n,n])
    
    for i in range(0,n):
        # Create vector of zeros except eps in position i
        e = np.zeros(n)
        e[i] = eps2
        
        # Evaluate gradient twice
        grad_plus = grad(x + e)
        grad_minus = grad(x - e)
        
        # Notice we are building the Hessian by column (or row)
        H[:,i] = (grad_plus - grad_minus)/(2*eps2)

    return H

## Linear algebra calculation
def xxT(u):
    '''
    Calculates u*u.T to circumvent limitation with SciPy
    
    Arguments:
    u - numpy 1D array
    
    Returns:
    u*u.T
    
    Assume u is a nx1 vector.
    Recall: NumPy does not distinguish between row or column vectors
    
    u.dot(u) returns a scalar. This functon returns an nxn matrix.
    '''
    
    n = len(u)
    A = np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,n):
            A[i,j] = u[i]*u[j]
    
    return A

## Analyze Hessian
def analyze_hes(B):
    print(B,"\n")
    
    l = linalg.eigvals(B)
    print("Eigenvalues: ",l,"\n")


def alg1_sr1(x0,calc_f,calc_grad,eps1,eps2,verbose=False,max_iter=1000):
    '''
    Arguments:
        x0 - starting point
        calc_f - funcation that calculates f(x)
        calc_grad - function that calculates gradient(x)
    
    Outputs:
        x - iteration history of x (decision variables)
        f - iteration history of f(x) (objective value)
        p - iteration history of p (steps)
        B - Hessian approximation
    '''
    
    # Allocate outputs as empty lists
    x = []
    f = []
    p = []
    grad = []
    B = []
    
    # Store starting point
    x.append(x0)
    k = 0
    
    flag = True
    
    print("Iter. \tf(x) \t\t||grad(x)|| \t||p|| \t\tmin(lambda)")
    
    while flag and k < max_iter:
        # Evaluate f(x) at current iteration
        f.append(calc_f(x[k]))
        
        # Evaluate gradient
        grad.append(calc_grad(x[k]))
        
        if(check_nan(grad[k])):
            print("WARNING: gradiant calculation returned NaN")
            break
        
        if verbose:
            print("\n")
            print("k = ",k)
            print("x = ",x[k])
            print("grad = ",grad[k])

        
        # Update Hessian approximation
        if k == 0:
            # Initialize with identity
            B.append(np.eye(len(x0)))

        else:
            # Change in x
            s = x[k] - x[k-1]

            # Change in gradient
            y = grad[k] - grad[k-1]

            # SR1 formulation
            u = y - B[k-1].dot(s)
            denom = (u).dot(s)
            
            # Formula: dB = u * u.T / (u.T * s) if u is a column vector.
            dB = xxT(u)/denom
            
            if(verbose):
                print("s = ",s)
                print("y = ",y)
                print("SR1 update denominator, (y-B[k-1]*s).T*s = ",denom)
                print("SR1 update u = ",u)
                print("SR1 update u.T*u/(u.T*s) = \n",dB)
            
            B.append(B[k-1] + dB)

        if verbose:
            print("B = \n",B[k])
            
        if(check_nan(B[k])):
            print("WARNING: Hessian update returned NaN")
            break
            
        c = linalg.cond(B[k])
        if c > 1E12:
            flag = False
            print("Warning: Hessian approximation is near singular.")
            print("B[k] = \n",B[k])
        
        else:
            # Calculate step
            p.append(linalg.solve(B[k],-grad[k]))

            if verbose:
                print("p = ",p[k])

            # Take step
            x.append(x[k] + p[k])

            # Calculate norms
            norm_p = linalg.norm(p[k])
            norm_g = linalg.norm(grad[k])

            # Calculate eigenvalues of Hessian (for display only)
            ev = np.real(linalg.eigvals(B[k]))

            # print("k = ",k,"\t"f[k],"\t",norm_g,"\t",norm_p)
            print(k,'  \t{0: 1.4e} \t{1:1.4e} \t{2:1.4e} \t{3: 1.4e}'.format(f[k],norm_g,norm_p,np.min(ev)))

            # Check convergence criteria
            flag = (norm_p > eps1) and (norm_g > eps2)

            # Update iteration counter
            k = k + 1
            
    print("Done.")
    print("x* = ",x[-1])
    
    return x,f,p,B