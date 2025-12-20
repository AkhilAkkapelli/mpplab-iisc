import random
import sys

# Constants for this case
BASE= 2**32
NCOEFFS= 4
NTESTS= 500000

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

def to_coeffs_mixed(val):
    '''
    Generates mixed-sign coefficients for a value
    '''
    # 1. Cannonical Positive coefficients
    absval= abs(val)
    coeffs= to_coeffs_clean(absval)

    # 2. Apply sign
    if val < 0:
        coeffs= [-c for c in coeffs]

    # 3. Add noise to generate some mixed coeffs
    # add (k*BASE) to c_i, subtract k from c_{i+1}
    for i in range(NCOEFFS-1):
        noise= random.randint(-2**28,2**28)
        coeffs[i] += noise*BASE
        coeffs[i+1] -= noise
    
    return coeffs

def generate_file():
    print(f"Generating {NTESTS} test cases for normalization into 'test_cases.txt'...")

    with open('./bin/mpi_normalization_tests.txt', "w") as f:
        # write header: num tests
        f.write(f"{NTESTS}\n")

        for _ in range(NTESTS):
            ## generate a random 124-bit number
            val= random.randint(0, 2**124)
            if random.choice([True, False]):
                val= -val
            
            # convert to coeffs
            input_coeffs= to_coeffs_mixed(val)
            output_coeffs=to_coeffs_clean(abs(val))
            sign= 1 if val >= 0 else -1

            # line up to be written
            # format: In0 In1 In2 In3 sign Out0 Out1 Out2 Out3
            line_items= input_coeffs + [sign] + output_coeffs
            line_str=   " ".join(map(str, line_items))
            f.write(line_str + "\n")

if __name__ == "__main__":
    generate_file()