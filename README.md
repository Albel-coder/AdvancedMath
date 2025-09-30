# AdvancedMath

AdvancedMath is a modern C++ library for advanced mathematical computations, providing convenient and efficient tools for working with complex numbers, linear algebra, mathematical analysis, and computational geometry.

## 🚀 Features

### ✅ Implemented

· Complex Numbers - full support for arithmetic operations, mathematical functions, and transformations
  · Algebraic Form (Complex class) - traditional real/imaginary representation
  · Trigonometric Form (ComplexTrigonometric class) - magnitude/angle representation
  · All basic arithmetic operations (+, -, *, /)
  · Comparison and assignment operators
  · Magnitude, phase, and conjugate calculations
  · Stream I/O in convenient format
  · Easy conversion between algebraic and trigonometric forms

### 🔄 In Development

· Matrices and Linear Algebra - matrix operations, determinants, inverse matrices
· Mathematical Analysis - numerical integration, differentiation
· Computational Geometry - barycentric coordinates, geometric transformations
· Additional Mathematical Structures - quaternions, tensors

### 🆕 Trigonometric Form Features

The library now includes a dedicated ComplexTrigonometric class for working with complex numbers in trigonometric form:

### Key Features:

· Create from polar coordinates: magnitude and angle (radians or degrees)
· Efficient operations: multiplication, division, and exponentiation optimized for trigonometric form
· Root extraction: calculate all n-th roots of a complex number
· Easy conversion: seamless conversion between algebraic and trigonometric forms
· Angle normalization: automatic handling of angle periodicity

· Multiplication: O(1) in trigonometric form vs O(4) in algebraic form
· Division: O(1) in trigonometric form vs O(9) in algebraic form
· Exponentiation: O(1) in trigonometric form vs O(n) in algebraic form
· Root extraction: Natural and efficient in trigonometric form

### 🛠 Technical Features

· Pure C++17 - modern language features
· Template-based design - support for various numeric types
· High performance - minimal overhead, optimized operations
· Dual representation - choose the best form for your use case
· Header-only (planned) - easy integration into projects
· Fully documented code - clear API
· Mathematical constants - built-in support for π and other constants

### 🎯 Project Goals

· Create a universal mathematical library for scientific computations
· Provide intuitive API for complex mathematical operations
· Ensure high performance for resource-intensive calculations
· Support multiple mathematical representations
· Become a useful tool for students, researchers, and developers
· Optimize operations by choosing the most efficient representation

### 🔬 Use Cases

The trigonometric form is particularly useful for:

· Signal processing: working with phases and amplitudes
· Electrical engineering: AC circuit analysis
· Physics: wave functions and oscillations
· Computer graphics: rotations and transformations
· Control systems: frequency domain analysis