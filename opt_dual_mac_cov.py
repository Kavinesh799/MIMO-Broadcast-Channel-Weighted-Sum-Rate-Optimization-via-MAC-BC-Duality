import numpy as np
from scipy.optimize import minimize_scalar


def compute_gradient(H_list, Q_list, q_values, sigma2):
    """
    Compute the gradient of the weighted sum rate function with respect to Q_i.
    
    Parameters:
    H_list  : List of channel matrices H_i (each H_i is MxN)
    Q_list  : List of covariance matrices Q_i (each Q_i is NxN)
    q_values: List of weights q_k
    sigma2  : Noise variance
    
    Returns:
    grad_list: List of gradients ∇f_i for each Q_i
    """
    
    K = len(H_list)  # Number of users
    M = H_list[0].shape[0]  # Number of antennas transmitter
    N = H_list[0].shape[1]  # Number of antennas receiver
    
    # Compute sum term for first part
    interference_sum = sum(H_j @ Q_list[j] @ H_j.conj().T for j, H_j in enumerate(H_list))
    first_inv_term = np.linalg.inv(np.eye(M) + (1 / sigma2) * interference_sum)

    grad_list = []
    
    for i in range(K):
        H_i = H_list[i]
        
        # First term
        Q_i_grad = q_values[0] * (H_list[i].conj().T @ first_inv_term @ H_list[i])  
        
        # Second summation term
        sum_term = np.zeros_like(Q_i_grad)
        for k in range(1, i+1):
            sub_interference = sum(H_list[j] @ Q_list[j] @ H_list[j].conj().T for j in range(k, K))
            sub_inv_term = np.linalg.inv(np.eye(M) + (1 / sigma2) * sub_interference)
            sum_term += (q_values[k] - q_values[k - 1]) * (H_i.conj().T @ sub_inv_term @ H_i)
        
        Q_i_grad += sum_term
        grad_list.append(Q_i_grad)

    return grad_list


def find_dominant_eigen(gradients):
    """ Find dominant eigenvalue and eigenvector for each user """
    eigenvalues = []
    eigenvectors = []
    for grad in gradients:
        eigvals, eigvecs = np.linalg.eigh(grad)
        idx = np.argmax(eigvals)
        eigenvalues.append(eigvals[idx])
        eigenvectors.append(eigvecs[:, idx])
    return eigenvalues, eigenvectors

#sum rate or objective function 
#pg 98 eq 3.35
def sum_rate(Q_list, H_list, q_values,sigma2):
    """
    Compute the sum rate for a MIMO(BC) with MIMO(MAC) formula using duality
    
    Parameters:
    H_list  : List of channel matrices H_k (each H_k is MxN).
    Q_list  : List of dual MAC covariance matrices Q_k (each Q_k is NxN).
    q_values: List of weights q_k (QoS).
    sigma2  : Noise variance.

    Returns:
    sum_rate: Weighted sum rate for the MIMO BC.
    """

    K = len(H_list)  # Number of users
    M = H_list[0].shape[0]  # Number of antennas at the transmitter
    sum_rate = 0

    # Compute first term: q_1log2(det(I_m+(sum{1,K}(HQHconjT))*1/(sigma2)))
    interference = sum(H_list[j] @ Q_list[j] @ H_list[j].conj().T for j in range(K))
    first_term = np.eye(M) + (1 / sigma2) * interference
    sign, logdet = np.linalg.slogdet(first_term)
    if sign <= 0:
        # warning or skip iteration
        print("Warning: num_Matrix is not positive definite. Skipping this iteration.")
        return np.inf  # or return np.inf or a large penalty
    else:
        first_term = logdet / np.log(2)  # convert from ln to log2
        first_term = q_values[0]*first_term
           
    #print("first_term = ", first_term)
    
    # Compute second term: sum((q_k-q_k-1)log2(det(I_m+(sum{k,K}(HQHconjT))*1/(sigma2))))
    #qk is to be in ordered fashion (q0<=q1<=q2......)
    for k in range(1,K):
        interference_next = sum(H_list[j] @ Q_list[j] @ H_list[j].conj().T for j in range(k, K))
        second_term = np.eye(M) + (1 / sigma2) * interference_next
        #print("denterm = ", denom_term)
        # Compute log-det ratio
        sign, logdet = np.linalg.slogdet(second_term)
        if sign <= 0:
            # warning or skip iteration
            print("Warning: num_Matrix is not positive definite. Skipping this iteration.")
            continue  # or return np.inf or a large penalty
        else:
            second_term = logdet / np.log(2)  # convert from ln to log2
        # Weighted sum rate contribution
        sum_rate += (q_values[k] - q_values[k-1]) * second_term
    
    sum_rate += first_term
    return sum_rate


"""convex optimizing function"""
def f(t, Q_list, H_list, sigma2, i_star, P, v_i_star, q):
    """Function to evaluate sum rate for given Q matrices and t value"""
    new_Q = Q_list.copy()
    for i in range(len(new_Q)):
        new_Q[i]=t*new_Q[i]
    new_Q[i_star] = t * Q_list[i_star] + (1 - t) * P * np.outer(v_i_star, v_i_star.conj().T)
    return -sum_rate(new_Q,H_list, q,sigma2)  # Negative because we are maximizing

def optimize_t(Q_list, H_list, v_i, i_star, P, sigma2,qs):
    " Solve t* = arg max_t∈[0,1] f(...) "
    res = minimize_scalar(f,bounds=(0,1), args=(Q_list,H_list,sigma2,i_star, P, v_i,qs), method='bounded')
    t_star = res.x
    return t_star

def optimal_dual_mac_covariances(H_list, P, q_values, sigma2, max_iter=1000, tol=1e-6):
    """
    Implements the algorithm for computing optimal dual MAC covariances.

    Args:
        H_list: List of channel matrices H_k (each MxN for user k).
        P: Total power constraint.
        max_iter: Maximum number of iterations.
        tol: Convergence tolerance.

    Returns:
        Optimal covariance matrices Q_list.
    """
    K = len(H_list)
    N = H_list[0].shape[1]
    M = H_list[0].shape[0]

    # Step 1: Initialize Q_k^(n) = 0
    Q_list = [np.zeros((N, N), dtype=complex) for _ in range(K)]

    for n in range(max_iter):
        Q_prev = [Q.copy() for Q in Q_list]

        # Compute gradients
        gradients = compute_gradient(H_list, Q_list, q_values, sigma2)

        # dominant eigenvalues and eigenvectors
        eigenvalues, eigenvectors = find_dominant_eigen(gradients)

        # i* with max eigenvalue
        i_star = np.argmax(eigenvalues)

        # Solve for optimal t*
        t_star = optimize_t(Q_list, H_list, eigenvectors[i_star], i_star, P, sigma2, q_values)

        # Update Q
        v_i = eigenvectors[i_star]
        Q_list[i_star] = t_star * Q_list[i_star] + (1 - t_star) * P * np.outer(v_i, v_i.conj())
        for j in range(K):
            if j != i_star:
                Q_list[j] *= t_star  # Scale other Q matrices

        # Step 3: Convergence check
        if all(np.linalg.norm(Q_list[k] - Q_prev[k]) < tol for k in range(K)):
            #print("stopped iterations-diff below tol")
            break
        # if n==max_iter-1:
        #     print("stopped iterations-crossed max iter","avg_delta_Q=",sum(np.linalg.norm(Q_list[k] - Q_prev[k]) for k in range(K))/K)    
    return Q_list, t_star, i_star, eigenvectors, eigenvalues, gradients

