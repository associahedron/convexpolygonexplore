def make_stack_rec(w, d, stack_index, wsum): 
    n = len(w)
    h = n-d-wsum+1 # where wsum = sum(w[d+1:])
    vals = list(range(h)) # Entry in w[d] for each substack
    if stack_index[d]%2 == 1: # Reverse every other d-dim substack
        vals = reversed(vals)
    stack_index[d] += 1
    for val in vals:
        w[d] = val
        if d == 1: ## Base case
            w[0] = n-1-(w[d]+wsum) # Equation 1 re-arranged
            print(w) # This will print the codewords in Hamiltonian order
        else:
            make_stack_rec(w, d-1, stack_index, wsum+val)
def print_hamiltonian(n): # Print a Hamiltonian traversal of Kn
    make_stack_rec([0]*n, n-1, [0]*n, 0)
print_hamiltonian(6)
