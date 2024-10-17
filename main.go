package main

import (
	"fmt"
	"log"
	"os"

	"optimization/matrix"
	"optimization/optalg"
	"optimization/solution"
	"optimization/solver"
	"optimization/userfuns"
)

func main() {
	defer func() {
		if r := recover(); r != nil {
			fmt.Println("ERROR:")
			log.Println(r)
		}
	}()

	lab0()
}

func lab0() {
	epsilon := 1e-2
	Nmax := 10000
	lb, _ := matrix.NewMatrix(2, 1, -5)
	ub, _ := matrix.NewMatrix(2, 1, 5)
	a, _ := matrix.NewMatrix(2, 1, 0.0)
	a.SetElement(0, 0, -1)
	a.SetElement(1, 0, 2)

	opt, err := optalg.MC(userfuns.Ff0T, 2, lb, ub, epsilon, Nmax, a, nil)
	if err != nil {
		log.Fatalf("Error in MC optimization: %v", err)
	}
	fmt.Println(opt)

	solution.ClearCalls()

	Nmax = 1000
	epsilon = 1e-2
	lb = matrix.NewMatrixScalar(0)
	ub = matrix.NewMatrixScalar(5)
	tetaOpt := matrix.NewMatrixScalar(1)

	opt, err = optalg.MC(userfuns.Ff0R, 1, lb, ub, epsilon, Nmax, tetaOpt, nil)
	if err != nil {
		log.Fatalf("Error in MC optimization for pendulum: %v", err)
	}
	fmt.Println(opt)

	solution.ClearCalls()

	// Save simulation to CSV file
	Y0, _ := matrix.NewMatrix(2, 1, 0.0)
	MT, _ := matrix.NewMatrixFromArray([]float64{matrix.M2D(opt.X), 0.5})

	Y, err := solver.SolveODE(userfuns.DF0, 0, 0.1, 10, Y0, nil, MT)
	if err != nil {
		log.Fatalf("Error solving ODE: %v", err)
	}

	file, err := os.Create("symulacja_lab0.csv")
	if err != nil {
		log.Fatalf("Error creating file: %v", err)
	}
	defer file.Close()

	combinedMatrix, _ := matrix.HCat(Y[0], Y[1])
	file.WriteString(combinedMatrix.String())

	// Clean up
	Y[0] = nil
	Y[1] = nil
}

func lab1() {
	// TODO: Add lab1 code
}

func lab2() {
	// TODO: Add lab2 code
}

func lab3() {
	// TODO: Add lab3 code
}

func lab4() {
	// TODO: Add lab4 code
}

func lab5() {
	// TODO: Add lab5 code
}

func lab6() {
	// TODO: Add lab6 code
}
