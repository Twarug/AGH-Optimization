package solution

import (
	"fmt"
	"math"

	"optimization/matrix"
)

// Solution struct equivalent to the C++ class 'solution'
type Solution struct {
	X, Y, G, H, Ud *matrix.Matrix
	Flag           int
}

// Static call counters (similar to the static variables in C++)
var FCalls, GCalls, HCalls int

// ClearCalls resets the static call counters
func ClearCalls() {
	FCalls = 0
	GCalls = 0
	HCalls = 0
}

// Constructors

// NewSolutionScalar creates a solution initialized with a scalar
func NewSolutionScalar(L float64) *Solution {
	return &Solution{
		X:    matrix.NewMatrixScalar(L),
		Y:    matrix.NewMatrixScalar(math.NaN()),
		G:    matrix.NewMatrixScalar(math.NaN()),
		H:    matrix.NewMatrixScalar(math.NaN()),
		Ud:   matrix.NewMatrixScalar(math.NaN()),
		Flag: -1,
	}
}

// NewSolutionMatrix creates a solution from a matrix
func NewSolutionMatrix(A *matrix.Matrix) *Solution {
	return &Solution{
		X:    A,
		Y:    matrix.NewMatrixScalar(math.NaN()),
		G:    matrix.NewMatrixScalar(math.NaN()),
		H:    matrix.NewMatrixScalar(math.NaN()),
		Ud:   matrix.NewMatrixScalar(math.NaN()),
		Flag: -1,
	}
}

// NewSolutionArray creates a solution from an array of values
func NewSolutionArray(n int, A []float64) (*Solution, error) {
	mat, err := matrix.NewMatrixFromArray(A)
	if err != nil {
		return nil, fmt.Errorf("solution.NewSolutionArray: %v", err)
	}
	return &Solution{
		X:    mat,
		Y:    matrix.NewMatrixScalar(math.NaN()),
		G:    matrix.NewMatrixScalar(math.NaN()),
		H:    matrix.NewMatrixScalar(math.NaN()),
		Ud:   matrix.NewMatrixScalar(math.NaN()),
		Flag: -1,
	}, nil
}

// CopySolution creates a copy of an existing solution
func CopySolution(A *Solution) *Solution {
	ud := matrix.NewMatrixScalar(math.NaN())
	if !math.IsNaN(A.Ud.M[0][0]) {
		ud = A.Ud
	}
	return &Solution{
		X:    A.X,
		Y:    A.Y,
		G:    A.G,
		H:    A.H,
		Ud:   ud,
		Flag: A.Flag,
	}
}

// FitFun computes the fitness function
func (s *Solution) FitFun(ff func(*matrix.Matrix, *matrix.Matrix, *matrix.Matrix) (*matrix.Matrix, error), ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	FCalls++
	y, err := ff(s.X, ud1, ud2)
	if err != nil {
		return nil, fmt.Errorf("solution.FitFun: %v", err)
	}
	s.Y = y
	return y, nil
}

// Grad computes the gradient of the function
func (s *Solution) Grad(gf func(*matrix.Matrix, *matrix.Matrix, *matrix.Matrix) (*matrix.Matrix, error), ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	GCalls++
	g, err := gf(s.X, ud1, ud2)
	if err != nil {
		return nil, fmt.Errorf("solution.Grad: %v", err)
	}
	s.G = g
	return g, nil
}

// Hess computes the Hessian of the function
func (s *Solution) Hess(hf func(*matrix.Matrix, *matrix.Matrix, *matrix.Matrix) (*matrix.Matrix, error), ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	HCalls++
	h, err := hf(s.X, ud1, ud2)
	if err != nil {
		return nil, fmt.Errorf("solution.Hess: %v", err)
	}
	s.H = h
	return h, nil
}

// GetDim returns the dimension of the solution's X matrix
func GetDim(A *Solution) (int, error) {
	length, err := matrix.GetLen(A.X)
	if err != nil {
		return 0, fmt.Errorf("solution.GetDim: %v", err)
	}
	return length, nil
}

// String implements the stringer interface for Solution
func (s *Solution) String() string {
	output := fmt.Sprintf("x = %v\n", s.X)
	output += fmt.Sprintf("y = %v\n", s.Y)
	output += fmt.Sprintf("f_calls = %d\n", FCalls)
	if GCalls > 0 {
		output += fmt.Sprintf("g_calls = %d\n", GCalls)
	}
	if HCalls > 0 {
		output += fmt.Sprintf("h_calls = %d\n", HCalls)
	}
	output += fmt.Sprintf("Exit flag: %d\n", s.Flag)
	return output
}
