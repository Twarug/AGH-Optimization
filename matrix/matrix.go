package matrix

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"strings"
	"time"
)

type Matrix struct {
	n, m int
	M    [][]float64
}

func (mat *Matrix) String() string {
	var sb strings.Builder
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			sb.WriteString(fmt.Sprintf("%f", mat.M[i][j]))
			if j < mat.m-1 {
				sb.WriteString(",")
			}
		}
		sb.WriteString("\n")
	}
	return sb.String()
}

// Constructors

// NewMatrixScalar creates a 1x1 matrix with a scalar value
func NewMatrixScalar(value float64) *Matrix {
	mat := &Matrix{
		n: 1,
		m: 1,
		M: [][]float64{{value}},
	}
	return mat
}

// NewMatrix creates an n x m matrix initialized with a given value.
func NewMatrix(n, m int, defaultValue float64) (*Matrix, error) {
	if n <= 0 || m <= 0 {
		return nil, errors.New("matrix dimensions must be positive")
	}
	mat := &Matrix{
		n: n,
		m: m,
		M: make([][]float64, n),
	}
	for i := 0; i < n; i++ {
		mat.M[i] = make([]float64, m)
		for j := 0; j < m; j++ {
			mat.M[i][j] = defaultValue
		}
	}
	return mat, nil
}

// NewMatrixFromArray creates a matrix from a 1D array.
// The resulting matrix will have dimensions len(arr) x 1.
func NewMatrixFromArray(arr []float64) (*Matrix, error) {
	if len(arr) == 0 {
		return nil, errors.New("array length must be positive")
	}
	mat := &Matrix{
		n: len(arr),
		m: 1,
		M: make([][]float64, len(arr)),
	}
	for i, val := range arr {
		mat.M[i] = []float64{val}
	}
	return mat, nil
}

// NewMatrixFrom2DArray creates an n x m matrix from a 2D array.
func NewMatrixFrom2DArray(arr [][]float64) (*Matrix, error) {
	if len(arr) == 0 || len(arr[0]) == 0 {
		return nil, errors.New("matrix dimensions must be positive")
	}
	mat := &Matrix{
		n: len(arr),
		m: len(arr[0]),
		M: arr,
	}
	return mat, nil
}

// CopyMatrix creates a copy of a matrix
func CopyMatrix(mat *Matrix) *Matrix {
	newMat := &Matrix{
		n: mat.n,
		m: mat.m,
		M: make([][]float64, mat.n),
	}
	for i := 0; i < mat.n; i++ {
		newMat.M[i] = make([]float64, mat.m)
		copy(newMat.M[i], mat.M[i])
	}
	return newMat
}

// Matrix Methods

// GetElement gets the value at position (i, j).
func (mat *Matrix) GetElement(i, j int) (float64, error) {
	if i >= mat.n || j >= mat.m || i < 0 || j < 0 {
		return 0, errors.New("index out of range")
	}
	return mat.M[i][j], nil
}

// SetElement sets the value at position (i, j)
func (mat *Matrix) SetElement(i, j int, value float64) error {
	if i >= mat.n || j >= mat.m || i < 0 || j < 0 {
		return errors.New("index out of range")
	}
	mat.M[i][j] = value
	return nil
}

// SetCol sets the values of a column
func (mat *Matrix) SetCol(col *Matrix, j int) error {
	if j >= mat.m || j < 0 {
		return errors.New("column index out of range")
	}
	if col.n != mat.n || col.m != 1 {
		return errors.New("column dimensions must match")
	}
	for i := 0; i < mat.n; i++ {
		mat.M[i][j] = col.M[i][0]
	}
	return nil
}

// SetRow sets the values of a row
func (mat *Matrix) SetRow(row *Matrix, i int) error {
	if i >= mat.n || i < 0 {
		return errors.New("row index out of range")
	}
	if row.m != mat.m || row.n != 1 {
		return errors.New("row dimensions must match")
	}
	for j := 0; j < mat.m; j++ {
		mat.M[i][j] = row.M[0][j]
	}
	return nil
}

// AddCol adds a column to the matrix
func (mat *Matrix) AddCol(value float64) {
	for i := 0; i < mat.n; i++ {
		mat.M[i] = append(mat.M[i], value)
	}
	mat.m++
}

// AddRow adds a row to the matrix
func (mat *Matrix) AddRow(value float64) {
	newRow := make([]float64, mat.m)
	for i := range newRow {
		newRow[i] = value
	}
	mat.M = append(mat.M, newRow)
	mat.n++
}

// Arithmetic operations (Add, Subtract, Multiply, Divide)

// AddMatrix adds two matrices
func (mat *Matrix) AddMatrix(other *Matrix) (*Matrix, error) {
	if mat.n != other.n || mat.m != other.m {
		return nil, errors.New("matrix dimensions must match")
	}
	result, _ := NewMatrix(mat.n, mat.m, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			result.M[i][j] = mat.M[i][j] + other.M[i][j]
		}
	}
	return result, nil
}

// SubtractMatrix subtracts one matrix from another
func (mat *Matrix) SubtractMatrix(other *Matrix) (*Matrix, error) {
	if mat.n != other.n || mat.m != other.m {
		return nil, errors.New("matrix dimensions must match")
	}
	result, _ := NewMatrix(mat.n, mat.m, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			result.M[i][j] = mat.M[i][j] - other.M[i][j]
		}
	}
	return result, nil
}

// MultiplyMatrix multiplies two matrices
func (mat *Matrix) MultiplyMatrix(other *Matrix) (*Matrix, error) {
	if mat.m != other.n {
		return nil, errors.New("matrix dimensions are not compatible for multiplication")
	}
	result, _ := NewMatrix(mat.n, other.m, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < other.m; j++ {
			for k := 0; k < mat.m; k++ {
				result.M[i][j] += mat.M[i][k] * other.M[k][j]
			}
		}
	}
	return result, nil
}

// TransposeMatrix returns the transpose of the matrix
func (mat *Matrix) TransposeMatrix() *Matrix {
	result, _ := NewMatrix(mat.m, mat.n, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			result.M[j][i] = mat.M[i][j]
		}
	}
	return result
}

// Utility functions

// M2D converts a 1x1 matrix to a scalar
func M2D(mat *Matrix) float64 {
	if mat.n != 1 || mat.m != 1 {
		panic(errors.New("M2D: conversion is only valid for 1x1 matrices"))
	}
	return mat.M[0][0]
}

// Determinant calculates the determinant of a square matrix
func (mat *Matrix) Determinant() (float64, error) {
	if mat.n != mat.m {
		return 0, errors.New("matrix must be square")
	}
	if mat.n == 1 {
		return mat.M[0][0], nil
	}
	if mat.n == 2 {
		return mat.M[0][0]*mat.M[1][1] - mat.M[0][1]*mat.M[1][0], nil
	}
	det := 0.0
	for p := 0; p < mat.m; p++ {
		subMatrix, _ := mat.Minor(0, p)
		subDet, _ := subMatrix.Determinant()
		det += mat.M[0][p] * math.Pow(-1, float64(p)) * subDet
	}
	return det, nil
}

// Minor returns the matrix with the specified row and column removed
func (mat *Matrix) Minor(row, col int) (*Matrix, error) {
	if mat.n <= 1 || mat.m <= 1 {
		return nil, errors.New("matrix too small for minors")
	}
	minor, _ := NewMatrix(mat.n-1, mat.m-1, 0)
	for i := 0; i < mat.n; i++ {
		if i == row {
			continue
		}
		for j := 0; j < mat.m; j++ {
			if j == col {
				continue
			}
			ni, nj := i, j
			if i > row {
				ni--
			}
			if j > col {
				nj--
			}
			minor.M[ni][nj] = mat.M[i][j]
		}
	}
	return minor, nil
}

// Inverse calculates the inverse of a square matrix
func (mat *Matrix) Inverse() (*Matrix, error) {
	det, err := mat.Determinant()
	if err != nil {
		return nil, err
	}
	if det == 0 {
		return nil, errors.New("matrix is singular, no inverse exists")
	}
	cofactors, _ := NewMatrix(mat.n, mat.m, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			minor, _ := mat.Minor(i, j)
			minorDet, _ := minor.Determinant()
			cofactors.M[i][j] = math.Pow(-1, float64(i+j)) * minorDet
		}
	}
	adjugate := cofactors.TransposeMatrix()
	inverse, _ := adjugate.ScalarMultiply(1 / det)
	return inverse, nil
}

// ScalarMultiply multiplies a matrix by a scalar
func (mat *Matrix) ScalarMultiply(scalar float64) (*Matrix, error) {
	result, _ := NewMatrix(mat.n, mat.m, 0)
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			result.M[i][j] = mat.M[i][j] * scalar
		}
	}
	return result, nil
}

// Norm calculates the Euclidean norm of a column vector.
func (mat *Matrix) Norm() (float64, error) {
	if mat.m != 1 {
		return 0, errors.New("norm is defined only for column vectors")
	}
	sum := 0.0
	for i := 0; i < mat.n; i++ {
		sum += mat.M[i][0] * mat.M[i][0]
	}
	return math.Sqrt(sum), nil
}

// Concatenation methods

// HCat concatenates two matrices horizontally
func HCat(A, B *Matrix) (*Matrix, error) {
	if A.n != B.n {
		return nil, errors.New("matrices must have the same number of rows for horizontal concatenation")
	}
	result, _ := NewMatrix(A.n, A.m+B.m, 0)
	for i := 0; i < A.n; i++ {
		copy(result.M[i][:A.m], A.M[i])
		copy(result.M[i][A.m:], B.M[i])
	}
	return result, nil
}

// VCat concatenates two matrices vertically
func VCat(A, B *Matrix) (*Matrix, error) {
	if A.m != B.m {
		return nil, errors.New("matrices must have the same number of columns for vertical concatenation")
	}
	result, _ := NewMatrix(A.n+B.n, A.m, 0)
	for i := 0; i < A.n; i++ {
		copy(result.M[i], A.M[i])
	}
	for i := 0; i < B.n; i++ {
		copy(result.M[A.n+i], B.M[i])
	}
	return result, nil
}

// GetCol returns the specified column as a new matrix
func (mat *Matrix) GetCol(j int) (*Matrix, error) {
	if j >= mat.m || j < 0 {
		return nil, errors.New("column index out of range")
	}
	result, _ := NewMatrix(mat.n, 1, 0)
	for i := 0; i < mat.n; i++ {
		result.M[i][0] = mat.M[i][j]
	}
	return result, nil
}

// GetRow returns the specified row as a new matrix
func (mat *Matrix) GetRow(i int) (*Matrix, error) {
	if i >= mat.n || i < 0 {
		return nil, errors.New("row index out of range")
	}
	result, _ := NewMatrix(1, mat.m, 0)
	for j := 0; j < mat.m; j++ {
		result.M[0][j] = mat.M[i][j]
	}
	return result, nil
}

// GetSize returns the dimensions (rows, columns) of the matrix
func GetSize(mat *Matrix) (int, int) {
	return mat.n, mat.m
}

// GetLen returns the number of rows for a column vector matrix (n x 1)
func GetLen(mat *Matrix) (int, error) {
	if mat.m != 1 {
		return 0, errors.New("GetLen: length can only be returned for column vectors")
	}
	return mat.n, nil
}

// Random matrix generation

// RandMat generates a random matrix with uniform distribution
func RandMat(n, m int) (*Matrix, error) {
	mat, _ := NewMatrix(n, m, 0)
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			mat.M[i][j] = rng.Float64()
		}
	}
	return mat, nil
}

// RandnMat generates a random matrix with normal distribution
func RandnMat(n, m int) (*Matrix, error) {
	mat, _ := NewMatrix(n, m, 0)
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))

	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			mat.M[i][j] = rng.NormFloat64()
		}
	}
	return mat, nil
}

// Utility for printing matrices

// Display prints the matrix in a readable format
func (mat *Matrix) Display() {
	for i := 0; i < mat.n; i++ {
		for j := 0; j < mat.m; j++ {
			fmt.Printf("%f ", mat.M[i][j])
		}
		fmt.Println()
	}
}
