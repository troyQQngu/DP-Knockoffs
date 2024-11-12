# the algorithm:
# 1. input: X,y
# 2. generate knockoffs X_k
# 3. compute the gram matrix t([X X_k]) %*% [X X_k] and the inner product 
#    t([X X_k]) %*% y. 