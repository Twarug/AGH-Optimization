package solver

import (
	"errors"
	"math"

	"optimization/matrix"
)

func safeAddMatrix(a, b *matrix.Matrix) *matrix.Matrix {
	result, err := a.AddMatrix(b)
	if err != nil {
		panic(err)
	}
	return result
}

func safeScalarMultiply(a *matrix.Matrix, scalar float64) *matrix.Matrix {
	result, err := a.ScalarMultiply(scalar)
	if err != nil {
		panic(err)
	}
	return result
}

// SolveODE solves the ODE using the 4th-order Runge-Kutta method
func SolveODE(diff func(float64, *matrix.Matrix, *matrix.Matrix, *matrix.Matrix) (*matrix.Matrix, error),
	t0, dt, tend float64, Y0, ud1, ud2 *matrix.Matrix) ([]*matrix.Matrix, error) {

	N := int(math.Floor((tend - t0) / dt))
	if N < 2 {
		return nil, errors.New("time interval is not defined correctly")
	}

	if Y0.M == nil || Y0.M[0] == nil || len(Y0.M[0]) != 1 {
		return nil, errors.New("initial condition must be a column vector")
	}
	n := len(Y0.M)

	T := matrix.NewMatrixScalar(t0)
	Solution, err := matrix.NewMatrix(n, N, 0.0)
	if err != nil {
		return nil, err
	}

	err = Solution.SetCol(Y0, 0)
	if err != nil {
		return nil, err
	}

	for i := 1; i < N; i++ {
		T.M[0][0] = T.M[0][0] + dt
		prevY, err := Solution.GetCol(i - 1)
		if err != nil {
			return nil, err
		}

		k1, err := diff(T.M[0][0], prevY, ud1, ud2)
		if err != nil {
			return nil, err
		}
		k1 = safeScalarMultiply(k1, dt)

		k2, err := diff(T.M[0][0]+0.5*dt, safeAddMatrix(prevY, safeScalarMultiply(k1, 0.5)), ud1, ud2)
		if err != nil {
			return nil, err
		}
		k2 = safeScalarMultiply(k2, dt)

		k3, err := diff(T.M[0][0]+0.5*dt, safeAddMatrix(prevY, safeScalarMultiply(k2, 0.5)), ud1, ud2)
		if err != nil {
			return nil, err
		}
		k3 = safeScalarMultiply(k3, dt)

		k4, err := diff(T.M[0][0]+dt, safeAddMatrix(prevY, k3), ud1, ud2)
		if err != nil {
			return nil, err
		}
		k4 = safeScalarMultiply(k4, dt)

		for j := 0; j < n; j++ {
			newValue := prevY.M[j][0] + (k1.M[j][0]+2*k2.M[j][0]+2*k3.M[j][0]+k4.M[j][0])/6
			Solution.M[j][i] = newValue
		}
	}

	Solution = Solution.TransposeMatrix()

	return []*matrix.Matrix{T, Solution}, nil
}
