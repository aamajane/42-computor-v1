#!/usr/bin/python3

# Computor V1 - Polynomial Equation Solver
# This program solves polynomial equations of degree 0, 1, or 2
# It handles both standard and free-form input formats
# Features:
# - Reduced form display
# - Degree calculation
# - Solutions with discriminant analysis for degree 2
# - Support for complex solutions
# - Irreducible fraction display

import sys
import re
from fractions import Fraction

def normalize_expression(expression):
    """
    Convert free-form expressions to standard form (bonus part)
    
    Transforms expressions like "5 + 4X + X^2" into "5 * X^0 + 4 * X^1 + 1 * X^2"
    
    Handles:
    - Coefficients adjacent to X (5X → 5 * X)
    - Missing powers (X → X^1)
    - Standalone numbers (5 → 5 * X^0)
    - Mixed case (x → X)
    """
    # Replace patterns like "5X" with "5 * X"
    expression = re.sub(r'(\d+)([Xx])', r'\1 * \2', expression)
    
    # Replace "X" with "X^1" and add " * X^0" to standalone numbers
    parts = []
    for term in re.split(r'([+-])', expression):
        term = term.strip()
        if not term:
            continue
            
        if term in "+-":
            parts.append(term)
            continue
            
        if 'X' in term or 'x' in term:
            # Standardize to uppercase X
            term = term.replace('x', 'X')
            
            # Add power if missing
            if not re.search(r'X\^', term):
                term = term.replace('X', 'X^1')
                
            parts.append(term)
        else:
            # It's a constant term, add * X^0
            try:
                float(term)  # Check if it's a number
                parts.append(f"{term} * X^0")
            except ValueError:
                parts.append(term)  # Keep as is if not a number
    
    return ''.join(parts)

def parse_expression(expression):
    """
    Parse an expression like '5 * X^0 + 4 * X^1 - 9.3 * X^2'
    
    This method:
    1. Normalizes the expression using _normalize_expression()
    2. Splits the expression into individual terms
    3. Extracts the coefficient and power from each term
    4. Builds a dictionary mapping powers to their coefficients
    
    Returns a dictionary where keys are powers and values are coefficients
    """
    coefficients = {}
    
    # First normalize the expression to handle free-form entries (bonus part)
    expression = normalize_expression(expression)
    
    # Split by + or - signs, but keep the sign with the term
    # This clever trick transforms "a - b + c" into ["a", "-b", "c"]
    expression = expression.replace("-", "+-").replace("++", "+")
    if expression.startswith("+"):
        expression = expression[1:]
        
    terms = expression.split('+')
    terms = [term.strip() for term in terms if term.strip()]
    
    for term in terms:
        term = term.strip()
        if not term:
            continue
            
        # Parse term like "-9.3 * X^2"
        if '*' in term:
            parts = term.split('*', 1)
            try:
                # Standard order: coefficient * variable
                coefficient = float(parts[0].strip())
                power_part = parts[1].strip()
            except ValueError:
                # Handle reversed order like "X^2 * 5"
                try:
                    coefficient = float(parts[1].strip())
                    power_part = parts[0].strip()
                except ValueError:
                    raise ValueError(f"Invalid term format: {term}")
        else:
            # Handle terms like "X^2" (coefficient is implied 1)
            if term.startswith("-"):
                coefficient = -1.0
                power_part = term[1:].strip()
            else:
                coefficient = 1.0
                power_part = term.strip()
        
        # Extract power from the term
        if 'X^' in power_part:
            try:
                power = int(power_part.split('^')[1])
            except ValueError:
                raise ValueError(f"Invalid power in term: {term}")
        elif 'X' in power_part:
            power = 1  # X is X^1
        else:
            power = 0  # Constant term
        
        # Add or update coefficient for this power (combine like terms)
        if power in coefficients:
            coefficients[power] += coefficient
        else:
            coefficients[power] = coefficient
    
    return coefficients

def parse_equation(equation):
    """
    Parses a polynomial equation and reduces it to standard form.
    
    Steps:
    1. Split into left and right sides
    2. Parse each side to extract coefficients
    3. Move all terms to the left side (reduced form)
    4. Remove terms with coefficients close to zero
    
    Example: "5 * X^2 + 4 * X - 3 = X^2 + 2" becomes {0: -5, 1: 4, 2: 4}
    """
    # Split by equals sign
    sides = equation.split('=')
    if len(sides) != 2:
        raise ValueError("Invalid equation format: must contain exactly one equals sign")
    
    left_side, right_side = sides[0].strip(), sides[1].strip()
    
    # Parse both sides
    left_coefficients = parse_expression(left_side)
    right_coefficients = parse_expression(right_side)
    
    # Subtract right from left to get the reduced form
    # This effectively moves all terms to the left side of the equation
    for power, coefficient in right_coefficients.items():
        if power in left_coefficients:
            left_coefficients[power] -= coefficient
        else:
            left_coefficients[power] = -coefficient
    
    # Clean up coefficients close to zero (handle floating-point precision issues)
    for power in list(left_coefficients.keys()):
        if abs(left_coefficients[power]) < 1e-10:
            del left_coefficients[power]
    
    return left_coefficients

def reduced_form(coefficients):
    """
    Format the reduced form of the equation
    
    This takes the internal coefficients dictionary and formats it into a readable 
    polynomial equation in the standard form:
    ax^n + bx^(n-1) + ... + cx^1 + d = 0
    
    Returns the formatted equation as a string
    """
    if not coefficients:
        return "0 = 0"
    
    terms = []
    # Sort powers to display polynomial in standard form (ascending powers)
    for power in sorted(coefficients.keys()):
        coef = coefficients.get(power, 0)
        if abs(coef) < 1e-10:  # Skip terms with coefficient close to zero
            continue
            
        # Format the coefficient with the appropriate sign
        if terms and coef > 0:
            term = f"+ {coef}"  # Add plus sign for positive terms (except first)
        elif coef < 0:
            term = f"- {abs(coef)}"  # Format negative terms with minus sign
        else:
            term = f"{coef}"  # First term doesn't need a plus sign
        
        # Add the X^power part
        term = f"{term} * X^{power}"
        
        terms.append(term)
    
    if not terms:
        return "0 = 0"  # Special case: all coefficients are zero
    
    # Handle the case where the first term is positive (don't need +)
    if terms[0].startswith("+ "):
        terms[0] = terms[0][2:]
    
    return f"Reduced form: {' '.join(terms)} = 0"

def get_degree(coefficients):
    """
    Determine the polynomial degree
    
    The degree of a polynomial is the highest power with a non-zero coefficient.
    For example:
    - 3x^2 + 2x + 1 has degree 2
    - 5x + 3 has degree 1
    - 7 has degree 0
    
    Returns:
        int: The polynomial degree
    """
    # Find the highest power with non-zero coefficient
    degree = 0
    for power, coef in coefficients.items():
        if abs(coef) > 1e-10 and power > degree:
            degree = power
    return degree

def as_fraction(value):
    """
    Convert a float to a fraction if possible (bonus feature)
    
    This method attempts to convert a floating-point number to an
    irreducible fraction if the fraction closely approximates the
    original value.
    
    Examples:
    - 0.5 becomes "1/2"
    - 0.25 becomes "1/4"
    - 2.0 becomes "2" (integer)
    - 0.3333333333 becomes "1/3" (approximate)
    
    Returns:
        str: The number as a fraction string or the original number as string
    """
    try:
        # Use Python's Fraction class to convert to a rational number
        frac = Fraction(value).limit_denominator(1000)
        
        # Check if the fraction is close enough to the original value
        if abs(float(frac) - value) < 1e-10:
            # If denominator is 1, just return the numerator (an integer)
            if frac.denominator == 1:
                return str(frac.numerator)
            return str(frac)
    except:
        pass  # If conversion fails, fall back to the string representation
    
    return str(value)

def solve(coefficients):
    """
    Solve the polynomial equation based on its degree
    
    This method handles:
    - Degree 0: Either all real numbers are solutions or no solution exists
    - Degree 1: Linear equation with one solution: x = -c/b
    - Degree 2: Quadratic equation using the discriminant method:
        * Discriminant > 0: Two real solutions
        * Discriminant = 0: One real solution (double root)
        * Discriminant < 0: Two complex solutions
    
    Returns a formatted string with the solution(s)
    """
    degree = get_degree(coefficients)
    
    # Format output
    result = [f"Polynomial degree: {degree}"]
    
    # Per project requirement, we can only solve up to degree 2
    if degree > 2:
        result.append("The polynomial degree is strictly greater than 2, I can't solve.")
        return "\n".join(result)
    
    # Get coefficients for easier reading: ax² + bx + c = 0
    a = coefficients.get(2, 0)
    b = coefficients.get(1, 0)
    c = coefficients.get(0, 0)
    
    # Handle different cases based on the polynomial degree
    if degree == 0:
        # Constant polynomial: c = 0
        if c == 0:
            result.append("All real numbers are solutions.")
        else:
            result.append("This equation has no solution.")
    elif degree == 1:
        # Linear equation: bx + c = 0, solution: x = -c/b
        solution = -c / b
        result.append("The solution is:")
        result.append(as_fraction(solution))
    elif degree == 2:
        # Quadratic equation: ax² + bx + c = 0
        # Calculate discriminant: Δ = b² - 4ac
        discriminant = b*b - 4*a*c
        
        if discriminant == 0:
            # One real solution (double root): x = -b/(2a)
            solution = -b / (2*a)
            result.append("Discriminant is zero, the solution is:")
            result.append(as_fraction(solution))
        elif discriminant > 0:
            # Two real solutions: x = (-b ± √Δ)/(2a)
            sqrt_discriminant = discriminant ** 0.5
            solution1 = (-b + sqrt_discriminant) / (2*a)
            solution2 = (-b - sqrt_discriminant) / (2*a)
            result.append("Discriminant is strictly positive, the two solutions are:")
            result.append(as_fraction(solution1))
            result.append(as_fraction(solution2))
        elif discriminant < 0:
            # Two complex solutions: x = -b/(2a) ± i√(-Δ)/(2a)
            real_part = -b / (2*a)
            imag_part = ((-discriminant) ** 0.5) / (2*a)
            result.append("Discriminant is strictly negative, the two complex solutions are:")
            result.append(f"{as_fraction(real_part)} + {as_fraction(imag_part)}i")
            result.append(f"{as_fraction(real_part)} - {as_fraction(imag_part)}i")
    
    return "\n".join(result)

def main():
    if len(sys.argv) != 2:
        print("Usage: ./computor.py \"<equation>\"")
        sys.exit(1)

    try:
        coefficients = parse_equation(sys.argv[1])
        print(reduced_form(coefficients))
        print(solve(coefficients))
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
