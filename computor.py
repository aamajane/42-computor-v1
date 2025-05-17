#!/usr/bin/python3

import sys
import re
from fractions import Fraction

def normalize_term(term):
    """
    Convert free-form term to standard form using regex pattern matching
    
    Handles all variations including:
    - Standalone numbers (5.1 → 5.1*X^0)
    - Standalone variables (X → 1*X^1)
    - Missing coefficients (X^2 → 1*X^2)
    - Missing powers (5*X → 5*X^1)
    - No multiplication symbol (5X^2 → 5*X^2)
    """
    # If it's just a number (no X)
    if 'X' not in term:
        return f"{term}*X^0"
    
    # Match different parts of the term using regex
    # Group 1: coefficient
    # Group 2: X
    # Group 3: power
    match = re.match(r'^([-+]?\d*\.?\d*)?(\*)?X(\^[-+]?\d+)?$', term)
    
    if match:
        coef, multiply, power = match.groups()
        
        # Handle coefficient
        if not coef or coef in ['+', '-']:
            coef = f"{coef}1"
        elif coef == '':
            coef = '1'
            
        # Handle power
        if not power:
            power = '^1'
            
        # Ensure multiplication symbol
        return f"{coef}*X{power}"
    
    # Handle terms where X is in the middle (like 5X^2)
    # Group 1: coefficient
    # Group 2: power
    match = re.match(r'^([-+]?\d*\.?\d+)X(\^[-+]?\d+)?$', term)
    
    if match:
        coef, power = match.groups()
        if not power:
            power = '^1'
        return f"{coef}*X{power}"
    
    # If it doesn't match our patterns, return as is (will be handled by error checking later)
    return term

def parse_term(term):
    """
    Parse a single term like "5*X^2" into coefficient and power.
    """
    # Use regex to directly match the expected format (coefficient*X^power)
    # Group 1: coefficient
    # Group 2: power
    match = re.match(r'^([-+]?\d*\.?\d+)\*X\^([-+]?\d+)$', term)
    
    if not match:
        raise ValueError(f"Invalid term format: {term}")
    
    try:
        coefficient = float(match.group(1))
        power = int(match.group(2))
        return (coefficient, power)
    except ValueError:
        raise ValueError(f"Invalid numeric values in term: {term}")

def parse_expression(expression):
    """
    Parse a polynomial expression into a dictionary mapping powers to coefficients.
    """
    # Add + sign at beginning if needed
    if expression and expression[0] not in ['+', '-']:
        expression = '+' + expression
    
    # Split into terms using regex
    terms = re.findall(r'[+\-][^+\-]*', expression)
    
    # Parse each term
    coefficients = {}
    for term in terms:
        if not term:
            continue
        
        sign = -1 if term[0] == '-' else 1
        norm_term = normalize_term(term[1:])
        coef, power = parse_term(norm_term)
        if power in coefficients:
            coefficients[power] += coef * sign
        else:
            coefficients[power] = coef * sign
    
    return coefficients

def parse_equation(equation):
    """
    Parse a polynomial equation and reduce it to standard form.
    """
    # Remove all whitespace first
    equation = re.sub(r'\s+', '', equation)
    
    # Standardize variable case
    equation = equation.replace('x', 'X')
    
    # Validate characters
    valid_chars = set("0123456789.=-+*^X")
    for char in equation:
        if char not in valid_chars:
            raise ValueError(f"Invalid character in equation: '{char}'")
    
    # Split into left and right sides by equal sign
    sides = equation.split('=')
    if len(sides) != 2 or not sides[0] or not sides[1]:
        raise ValueError("Equation must be in the form 'left = right'")
    
    # Parse both sides
    left_coefficients = parse_expression(sides[0])
    right_coefficients = parse_expression(sides[1])
    
    # Subtract right from left to get the reduced form
    for power, coef in right_coefficients.items():
        if power in left_coefficients:
            left_coefficients[power] -= coef
        else:
            left_coefficients[power] = -coef
    
    # Clean up zero coefficients (accounting for floating point precision)
    # Create a copy of keys with list() to safely modify the dictionary during iteration
    for power in list(left_coefficients.keys()):
        if abs(left_coefficients[power]) < 1e-10:
            del left_coefficients[power]
    
    return left_coefficients

def reduced_form(coefficients):
    """
    Format the reduced form of the equation as a readable polynomial.
    """
    if not coefficients:
        return "Reduced form: 0 = 0"
    
    terms = []
    # Sort powers in descending order
    for i, power in enumerate(sorted(coefficients.keys(), reverse=True)):
        coef = coefficients[power]
        
        # Format with appropriate sign
        prefix = "" if i == 0 and coef > 0 else "+ " if coef > 0 else "- "
        value = abs(coef)
        
        # Create the term with its power
        term = f"{prefix}{value} * X^{power}"
        terms.append(term)
    
    return f"Reduced form: {' '.join(terms)} = 0"

def as_fraction(value):
    """
    Convert a float to a fraction if possible
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
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
