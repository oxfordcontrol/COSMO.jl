using PyCall
@pyimport osqp as osqp
@pyimport scipy.sparse as spar
@pyimport numpy as np

# Define problem data
P = spar.csc_matrix([[4, 1], [1, 2]])
q = np.array([1, 1])
A = spar.csc_matrix([[1, 1], [1, 0], [0, 1]])
l = np.array([1, 0, 0])
u = np.array([1, 0.7, 0.7])


prob = osqp.OSQP()
prob[:setup](P, q, A, l, u, alpha=1.0)
res = prob[:solve]()