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
    print(f"[Python] Generating {NTESTS} test cases for normalization into 'mpi_from_string_tests.txt'...")

    with open('./bin/mpi_from_string_tests.txt', "w") as f:
        # write header: num tests
        f.write(f"{NTESTS}\n")

        for _ in range(NTESTS):
            ## generate a random 64-bit number
            val= random.randint(-2**127 + 1, 2**127 - 1)
            
            # convert to coeffs
            output_coeffs=to_coeffs_clean(abs(val))
            sign= 1 if val >= 0 else -1
            output_coeffs= [sign*c for c in output_coeffs]

            # line up to be written
            # format: [val] inp_coeffs
            line_items= [val] + output_coeffs
            line_str=   " ".join(map(str, line_items))
            f.write(line_str + "\n")

if __name__ == "__main__":
    generate_file()