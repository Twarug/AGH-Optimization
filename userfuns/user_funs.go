package userfuns

import (
	"math"

	"optimization/matrix"
	"optimization/solver"
)

func Ff0T(x, ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	if ud1 == nil {
		ud1 = matrix.NewMatrixScalar(math.NaN())
	}
	if ud2 == nil {
		ud2 = matrix.NewMatrixScalar(math.NaN())
	}

	result, _ := matrix.NewMatrix(1, 1, 0.0)
	x0Diff := math.Pow(x.M[0][0]-ud1.M[0][0], 2)
	x1Diff := math.Pow(x.M[1][0]-ud1.M[1][0], 2)
	result.M[0][0] = x0Diff + x1Diff

	return result, nil
}

func Ff0R(x, ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	if ud1 == nil {
		ud1 = matrix.NewMatrixScalar(math.NaN())
	}
	if ud2 == nil {
		ud2 = matrix.NewMatrixScalar(math.NaN())
	}

	Y0, _ := matrix.NewMatrix(2, 1, 0.0)
	MT, _ := matrix.NewMatrixFromArray([]float64{matrix.M2D(x), 0.5})

	Y, err := solver.SolveODE(DF0, 0, 0.1, 10, Y0, ud1, MT)
	if err != nil {
		return nil, err
	}

	Solution := Y[1]
	n := len(Y[0].M)
	tetaMax := Solution.M[0][0]
	for i := 1; i < n; i++ {
		if Solution.M[i][0] > tetaMax {
			tetaMax = Solution.M[i][0]
		}
	}

	y := matrix.NewMatrixScalar(math.Abs(tetaMax - matrix.M2D(ud1)))

	return y, nil
}

func DF0(t float64, Y, ud1, ud2 *matrix.Matrix) (*matrix.Matrix, error) {
	m := 1.0
	l := 0.5
	b := 0.5
	g := 9.81
	I := m * math.Pow(l, 2)

	dY, _ := matrix.NewMatrix(2, 1, 0.0)
	dY.SetElement(0, 0, Y.M[1][0])

	force := 0.0
	if t <= ud2.M[1][0] {
		force = ud2.M[0][0]
	}
	angularAcc := (force - m*g*l*math.Sin(Y.M[0][0]) - b*Y.M[1][0]) / I
	dY.SetElement(1, 0, angularAcc)

	return dY, nil
}
