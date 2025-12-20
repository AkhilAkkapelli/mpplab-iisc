import random
import sys

# Constants for this case
BASE= 2**32
NCOEFFS= 4
NTESTS= 50000

def to_coeffs_clean(val):
    '''
    Convert positive value to coeffs
    '''
    assert val >= 0, "Please provide a non-negative value for `to_clean_coeffs`"
    coeffs= []
    for _ in range(NCOEFFS):
        coeffs.append(val % BASE)
        val //= BASE
    return coeffs

def generate_file():
    print(f"Generating {NTESTS} test cases for normalization into 'mpi_from_coeffs_tests.txt'...")

    with open('./bin/mpi_scalar_mult.txt', "w") as f:
        # write header: num tests
        f.write(f"{NTESTS}\n")

        for _ in range(NTESTS):
            ## generate a random 64-bit number
            val1= random.randint(-2**96+1, 2**96-1)
            
            # convert to coeffs
            c1=to_coeffs_clean(abs(val1))
            sign= 1 if val1 >= 0 else -1
            c1= [sign*c for c in c1]

            # scalar
            val2= random.randint(-2**31+1, 2**31 - 1)

            # output
            val_out= val1*val2
            c2= to_coeffs_clean(abs(val_out))
            sign= -1 if val_out < 0 else 1
            c2= [sign*c for c in c2]

            # line up to be written
            # format: In0 In1 In2 In3 scalar Ex0 Ex1 Ex2 Ex3
            line_items= c1 + [val2] + c2
            line_str=   " ".join(map(str, line_items))
            f.write(line_str + "\n")

if __name__ == "__main__":
    generate_file()