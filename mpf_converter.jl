using Printf

struct MPF
    coeffs::Vector{Int64}  
    exponent::Int32
    is_negative::Bool
end

const BASE = BigInt(2)^32

function parse_mpf_string(s::String, num_coeffs::Int)
    s = strip(s)
    
    # Count significant digits in the input to determine required precision
    # Remove sign, decimal point, and scientific notation
    digits_only = replace(s, r"[+\-\.eE]" => "")
    num_digits = length(digits_only)
    
    # Binary bits needed: decimal_digits / log10(2) ≈ decimal_digits * 3.32
    # Add safety margin
    input_precision = ceil(Int, num_digits * 3.4) + 128
    
    # Use the larger of: what we need to store, or what's in the input
    required_precision = max(input_precision, num_coeffs * 32 + 128)
    
    # Set precision BEFORE parsing
    old_prec = Base.precision(BigFloat)
    setprecision(BigFloat, required_precision)
    
    val = parse(BigFloat, s)
    
    is_negative = val < 0
    val = abs(val)
    
    if val == 0
        setprecision(BigFloat, old_prec)
        return MPF(zeros(Int64, num_coeffs), Int32(0), false), Int32(0)
    end
    
    # Use frexp to get the binary exponent
    # frexp returns (mantissa, exp) where mantissa ∈ [0.5, 1.0) and val = mantissa * 2^exp
    mantissa_norm, exp_binary = frexp(val)
    
    # We want to store the number as: integer_mantissa * 2^stored_exponent
    # where integer_mantissa fits in num_coeffs * 32 bits
    
    # Convert mantissa to an integer by scaling up
    # mantissa_norm is in [0.5, 1.0), so mantissa_norm * 2^(num_coeffs*32) gives us the integer
    total_bits = num_coeffs * 32
    scaled_mantissa = mantissa_norm * BigFloat(2.0)^total_bits
    
    mantissa_int = BigInt(round(scaled_mantissa))
    
    # The stored exponent should account for the scaling we did
    stored_exponent = exp_binary - total_bits
    
    # Extract coefficients in little-endian order (LSB first)
    coeffs = zeros(Int64, num_coeffs)
    temp = mantissa_int
    for i in 1:num_coeffs
        coeffs[i] = Int64(temp & (BASE - 1))
        temp >>= 32  # Right shift by 32 bits
    end
    
    setprecision(BigFloat, old_prec)
    
    return MPF(coeffs, Int32(stored_exponent), is_negative), Int32(stored_exponent)
end

function to_string(mpf::MPF)
    # Check if zero
    if all(c -> c == 0, mpf.coeffs)
        return "0.0"
    end
    
    # Reconstruct mantissa as BigInt from coefficients (little-endian)
    mantissa = BigInt(0)
    num_coeffs = length(mpf.coeffs)
    
    # Build from most significant to least significant
    for i in num_coeffs:-1:1
        mantissa = (mantissa << 32) | BigInt(mpf.coeffs[i])
    end
    
    # Set precision high enough for accurate reconstruction
    old_prec = Base.precision(BigFloat)
    required_precision = max(512, num_coeffs * 32 + 128)
    setprecision(BigFloat, required_precision)
    
    # Calculate the actual value: mantissa * 2^exponent
    result = BigFloat(mantissa) * BigFloat(2.0)^mpf.exponent
    
    if mpf.is_negative
        result = -result
    end
    
    # Calculate decimal precision based on binary precision
    # Each coefficient stores 32 bits, so we can represent ~9.6 decimal digits per coefficient
    decimal_precision = ceil(Int, num_coeffs * 32 * 0.30103)
    
    # Format with maximum precision
    format_str = Printf.Format("%.$(decimal_precision)e")
    output = Printf.format(format_str, result)
    
    # Restore original precision
    setprecision(BigFloat, old_prec)
    
    return output
end

function test_conversion(s::String, num_coeffs::Int=8)
    println("Testing: $s with num_coeffs=$num_coeffs")
    m, e = parse_mpf_string(s, num_coeffs)
    println("  Coeffs: $(m.coeffs)")
    println("  Exponent: $e")
    println("  Is negative: $(m.is_negative)")
    result = to_string(m)
    println("  Reconstructed: $result")
    
    # Check accuracy with high precision
    old_prec = Base.precision(BigFloat)
    
    # Use very high precision for comparison
    digits_only = replace(s, r"[+\-\.eE]" => "")
    num_digits = length(digits_only)
    comparison_precision = ceil(Int, num_digits * 3.4) + 256
    
    setprecision(BigFloat, comparison_precision)
    original = parse(BigFloat, s)
    reconstructed = parse(BigFloat, result)
    rel_error = abs(original) > 0 ? abs(original - reconstructed) / abs(original) : abs(reconstructed)
    println("  Relative error: $(rel_error)")
    setprecision(BigFloat, old_prec)
    println()
end

# Example usage
# test_conversion("123456789012345678901234567890.12345678901234567890123456789012345678901234567890", 4)
# test_conversion("123456789012345678901234567890.12345678901234567890123456789012345678901234567890", 5)
# test_conversion("123456789012345678901234567890.12345678901234567890123456789012345678901234567890", 6)
# test_conversion("123456789012345678901234567890.12345678901234567890123456789012345678901234567890", 7)
# test_conversion("123456789012345678901234567890.12345678901234567890123456789012345678901234567890", 8)
# test_conversion("098765432109876543210987654321.09876543210987654321098765432109876543210987654321", 4)
# test_conversion("098765432109876543210987654321.09876543210987654321098765432109876543210987654321", 5)
# test_conversion("098765432109876543210987654321.09876543210987654321098765432109876543210987654321", 6)
# test_conversion("098765432109876543210987654321.09876543210987654321098765432109876543210987654321", 7)
# test_conversion("098765432109876543210987654321.09876543210987654321098765432109876543210987654321", 8)
