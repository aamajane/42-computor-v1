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

def validate_equation(equation):
    """
    Validate if the equation follows the required format.
    - Must have exactly one equals sign
    - Must only contain allowed characters: digits, X/x, +, -, *, ^, spaces
    - Powers must be non-negative integers
    """
    # Check for exactly one equals sign
    if equation.count('=') != 1:
        raise ValueError("Equation must contain exactly one equals sign")
    
    # Validate characters
    valid_chars = set("0123456789.+-*/^xX= \t\n")
    for char in equation:
        if char not in valid_chars:
            raise ValueError(f"Invalid character in equation: '{char}'")
            
    return True

def parse_term(term):
    """
    Parse a single term like "5*X^2", "-X^3", or "42" into coefficient and power.
    Returns (coefficient, power) tuple.
    """
    # Handle empty terms
    if not term.strip():
        return None
    
    # Determine sign
    sign = 1
    if term.startswith('-'):
        sign = -1
        term = term[1:].strip()
    elif term.startswith('+'):
        term = term[1:].strip()
    
    # Handle constants (no variable)
    if 'x' not in term.lower():
        try:
            return (sign * float(term), 0)
        except ValueError:
            raise ValueError(f"Invalid constant term: {term}")
    
    # Handle terms with variables
    coefficient = sign
    power = 1  # Default power
    
    # Extract coefficient
    if '*' in term:
        parts = term.split('*', 1)
        if parts[0].strip():
            try:
                coefficient = sign * float(parts[0].strip())
            except ValueError:
                raise ValueError(f"Invalid coefficient: {parts[0]}")
        term = parts[1].strip()
    elif term.lower() != 'x' and term.lower()[0] not in 'x':
        # Handle format like "5x" without *
        match = re.match(r'^(\d+\.?\d*)([xX].*)', term)
        if match:
            coefficient = sign * float(match.group(1))
            term = match.group(2)
    
    # Validate that term starts with x/X now
    if not term.lower().startswith('x'):
        raise ValueError(f"Invalid term format: {term}")
    
    # Extract power
    if '^' in term:
        parts = term.split('^', 1)
        try:
            power = int(parts[1].strip())
            if power < 0:
                raise ValueError("Negative powers are not supported")
        except ValueError:
            raise ValueError(f"Invalid power: {parts[1]}")
    
    return (coefficient, power)

def parse_expression(expression):
    """
    Parse a polynomial expression into a dictionary mapping powers to coefficients.
    """
    # Remove all spaces
    expression = expression.replace(" ", "")
    
    # Add + sign at beginning if needed
    if expression and expression[0] not in ['+', '-']:
        expression = '+' + expression
    
    # Split into terms
    terms = []
    i = 0
    while i < len(expression):
        if expression[i] in ['+', '-']:
            # Find the next sign or end of expression
            j = i + 1
            while j < len(expression) and expression[j] not in ['+', '-']:
                j += 1
            terms.append(expression[i:j])
            i = j
        else:
            i += 1
    
    # Parse each term
    coefficients = {}
    for term in terms:
        if not term.strip():
            continue
            
        try:
            coef, power = parse_term(term)
            if power in coefficients:
                coefficients[power] += coef
            else:
                coefficients[power] = coef
        except ValueError as e:
            raise ValueError(f"Error parsing term '{term}': {str(e)}")
    
    return coefficients

def parse_equation(equation):
    """
    Parse a polynomial equation and reduce it to standard form.
    """
    # Validate equation format
    validate_equation(equation)
    
    # Split into left and right sides
    sides = equation.split('=')
    left_side, right_side = sides[0].strip(), sides[1].strip()
    
    # Parse both sides
    left_coefficients = parse_expression(left_side)
    right_coefficients = parse_expression(right_side)
    
    # Subtract right from left to get the reduced form
    result = {}
    for power, coef in left_coefficients.items():
        result[power] = coef
    
    for power, coef in right_coefficients.items():
        if power in result:
            result[power] -= coef
        else:
            result[power] = -coef
    
    # Clean up zero coefficients (accounting for floating point precision)
    result = {power: coef for power, coef in result.items() if abs(coef) > 1e-10}
    
    return result

def reduced_form(coefficients):
    """
    Format the reduced form of the equation
    
    This takes the internal coefficients dictionary and formats it into a readable 
    polynomial equation in the standard form:
    ax^n + bx^(n-1) + ... + cx^1 + d = 0
    
    Returns the formatted equation as a string
    """
    if not coefficients:
        return "Reduced form: 0 = 0"
    
    terms = []
    # Sort powers to display polynomial in standard form (ascending powers)
    for power in sorted(coefficients.keys()):
        coef = coefficients.get(power, 0)
        
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
        return "Reduced form: 0 = 0"  # Special case: all coefficients are zero
    
    # Handle the case where the first term is positive (don't need +)
    if terms[0].startswith("+ "):
        terms[0] = terms[0][2:]
    
    return f"Reduced form: {' '.join(terms)} = 0"

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
    degree = max(coefficients.keys()) if coefficients else 0
    
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
        print("Error: Invalid equation format.")
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
