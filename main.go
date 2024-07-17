package main

import (
	"log"

	"github.com/hajimehoshi/ebiten/v2"
)

const (
	N = 128
)

func main() {
	fluid := NewFluid(N, 0, 0)

	for i := N/2 - 10; i <= N/2+10; i++ {
		for j := 1; j <= 3; j++ {
			fluid.USource[IX(i, j)] = 1
		}
	}
	for i := N/2 - 10; i <= N/2+10; i++ {
		for j := N - 2; j <= N; j++ {
			fluid.USource[IX(i, j)] = -1
		}
	}
	for i := N/2 - 3; i <= N/2+3; i++ {
		fluid.Sources[IX(i, 1)] = 100
		fluid.Sources[IX(i, N)] = 100
	}

	g := &Game{
		fluid: fluid,
	}

	ebiten.SetWindowSize(900, 900)
	ebiten.SetWindowTitle("Testing")
	if err := ebiten.RunGame(g); err != nil {
		log.Fatal(err)
	}
}

func (g *Game) Update() error {
	g.fluid.add_source()
	g.fluid.diffuse()
	g.fluid.advect()
	g.fluid.diffuse_velocity()
	g.fluid.advect_velocity()
	g.fluid.project()
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	if g.pixels == nil {
		g.pixels = make([]byte, N*N*4)
	}
	var index int
	for i := 1; i <= N; i++ {
		for j := 1; j <= N; j++ {
			value := byte(int(min(g.fluid.Density[IX(i, j)], 1.) * 255))
			index = (i-1)*N*4 + (j-1)*4
			g.pixels[index] = value
			g.pixels[index+1] = value
			g.pixels[index+2] = value
			g.pixels[index+3] = value

		}
	}
	screen.WritePixels(g.pixels)
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (int, int) {
	return N, N
}

type Game struct {
	pixels []byte
	fluid  Fluid
}

type Fluid struct {
	Size       int
	N          int
	Dt         float64
	Diff       float64
	Visc       float64
	Sources    []float64
	Density    []float64
	NewDensity []float64
	U          []float64
	V          []float64
	NewU       []float64
	NewV       []float64
	USource    []float64
	VSource    []float64
	Div        []float64
	P          []float64
}

func NewFluid(size int, diff, visc float64) Fluid {
	var arraySize int = (size + 2) * (size + 2)
	var fluid Fluid = Fluid{
		Size:       size,
		Dt:         float64(1. / 30.),
		Diff:       diff,
		Visc:       visc,
		Sources:    make([]float64, arraySize),
		Density:    make([]float64, arraySize),
		NewDensity: make([]float64, arraySize),
		U:          make([]float64, arraySize),
		V:          make([]float64, arraySize),
		NewU:       make([]float64, arraySize),
		NewV:       make([]float64, arraySize),
		USource:    make([]float64, arraySize),
		VSource:    make([]float64, arraySize),
		Div:        make([]float64, arraySize),
		P:          make([]float64, arraySize),
	}
	return fluid
}

func IX(i, j int) int {
	return i*(N+2) + j
}

func (f *Fluid) add_source() {
	for i, v := range f.Sources {
		f.Density[i] += v * f.Dt
		f.U[i] += f.USource[i] * f.Dt
		f.V[i] += f.VSource[i] * f.Dt
	}
}

func (f *Fluid) diffuse() {
	if f.Diff == 0 {
		return
	}
	a := float64(f.Dt * f.Diff * N * N)
	for k := 0; k < 20; k++ {
		for i := 1; i <= N; i++ {
			for j := 1; j <= N; j++ {
				f.NewDensity[IX(i, j)] = (f.Density[IX(i, j)] + a*(f.NewDensity[IX(i+1, j)]+f.NewDensity[IX(i-1, j)]+f.NewDensity[IX(i, j+1)]+f.NewDensity[IX(i, j-1)])) / (float64(1) + 4*a)
			}
		}
	}
	set_bnd(N, 0, &f.NewDensity)
	var aux []float64
	aux = f.Density
	f.Density = f.NewDensity
	f.NewDensity = aux
}

func (f *Fluid) advect() {
	var i0, j0, i1, j1 int
	var x, y, s0, t0, s1, t1, dt0 float64
	dt0 = f.Dt * float64(f.Size)
	for i := 1; i <= N; i++ {
		for j := 1; j <= N; j++ {
			x = float64(i) - dt0*f.V[IX(i, j)]
			y = float64(j) - dt0*f.U[IX(i, j)]
			if x < 0.5 {
				x = 0.5
			} else if x > N+0.5 {
				x = N + 0.5
			}
			if y < 0.5 {
				y = 0.5
			} else if y > N+0.5 {
				y = N + 0.5
			}
			i0 = int(x)
			i1 = i0 + 1
			j0 = int(y)
			j1 = j0 + 1
			s1 = x - float64(i0)
			s0 = 1 - s1
			t1 = y - float64(j0)
			t0 = 1 - t1
			f.NewDensity[IX(i, j)] = s0*(t0*f.Density[IX(i0, j0)]+t1*f.Density[IX(i0, j1)]) + s1*(t0*f.Density[IX(i1, j0)]+t1*f.Density[IX(i1, j1)])
		}
	}
	var aux []float64
	aux = f.Density
	f.Density = f.NewDensity
	f.NewDensity = aux
}

func (f *Fluid) advect_velocity() {
	var i0, j0, i1, j1 int
	var x, y, s0, t0, s1, t1, dt0 float64
	dt0 = f.Dt * float64(f.Size)
	for i := 1; i <= N; i++ {
		for j := 1; j <= N; j++ {
			x = float64(i) - dt0*f.V[IX(i, j)]
			y = float64(j) - dt0*f.U[IX(i, j)]
			if x < 0.5 {
				x = 0.5
			} else if x > N+0.5 {
				x = N + 0.5
			}
			if y < 0.5 {
				y = 0.5
			} else if y > N+0.5 {
				y = N + 0.5
			}
			i0 = int(x)
			i1 = i0 + 1
			j0 = int(y)
			j1 = j0 + 1
			s1 = x - float64(i0)
			s0 = 1 - s1
			t1 = y - float64(j0)
			t0 = 1 - t1
			f.NewU[IX(i, j)] = s0*(t0*f.U[IX(i0, j0)]+t1*f.U[IX(i0, j1)]) + s1*(t0*f.U[IX(i1, j0)]+t1*f.U[IX(i1, j1)])
			f.NewV[IX(i, j)] = s0*(t0*f.V[IX(i0, j0)]+t1*f.V[IX(i0, j1)]) + s1*(t0*f.V[IX(i1, j0)]+t1*f.V[IX(i1, j1)])
		}
	}
	var aux []float64
	aux = f.U
	f.U = f.NewU
	f.NewU = aux
	aux = f.V
	f.V = f.NewV
	f.NewV = aux
}

func (f *Fluid) project() {
	var i, j, k int
	var h float64 = 1.0 / N
	for i = 1; i <= N; i++ {
		for j = 1; j <= N; j++ {
			f.Div[IX(i, j)] = -0.5 * h * (f.V[IX(i+1, j)] - f.V[IX(i-1, j)] + f.U[IX(i, j+1)] - f.U[IX(i, j-1)])
			f.P[IX(i, j)] = 0
		}
	}
	set_bnd(N, 0, &f.Div)
	set_bnd(N, 0, &f.P)
	for k = 0; k < 20; k++ {
		for i = 1; i <= N; i++ {
			for j = 1; j <= N; j++ {
				f.P[IX(i, j)] = (f.Div[IX(i, j)] + f.P[IX(i-1, j)] + f.P[IX(i+1, j)] + f.P[IX(i, j-1)] + f.P[IX(i, j+1)]) / 4
			}
		}
		set_bnd(N, 0, &f.P)
	}
	for i = 1; i <= N; i++ {
		for j = 1; j <= N; j++ {
			f.V[IX(i, j)] -= 0.5 * (f.P[IX(i+1, j)] - f.P[IX(i-1, j)]) / h
			f.U[IX(i, j)] -= 0.5 * (f.P[IX(i, j+1)] - f.P[IX(i, j-1)]) / h
		}
	}
	set_bnd(N, 1, &f.U)
	set_bnd(N, 2, &f.V)
}

func (f *Fluid) diffuse_velocity() {
	if f.Visc == 0 {
		return
	}
	a := float64(f.Dt * f.Visc * N * N)
	for k := 0; k < 20; k++ {
		for i := 1; i <= N; i++ {
			for j := 1; j <= N; j++ {
				f.NewU[IX(i, j)] = (f.U[IX(i, j)] + a*(f.NewU[IX(i+1, j)]+f.NewU[IX(i-1, j)]+f.NewU[IX(i, j+1)]+f.NewU[IX(i, j-1)])) / (float64(1) + 4*a)
				f.NewV[IX(i, j)] = (f.V[IX(i, j)] + a*(f.NewV[IX(i+1, j)]+f.NewV[IX(i-1, j)]+f.NewV[IX(i, j+1)]+f.NewV[IX(i, j-1)])) / (float64(1) + 4*a)
			}
		}
	}
	set_bnd(N, 0, &f.NewU)
	set_bnd(N, 0, &f.NewV)
	var aux []float64
	aux = f.U
	f.U = f.NewU
	f.NewU = aux
	aux = f.V
	f.V = f.NewV
	f.NewU = aux
}

func set_bnd(N int, b int, x *[]float64) {
	for i := 1; i <= N; i++ {
		if b == 1 {
			(*x)[IX(i, 0)] = -(*x)[IX(i, 1)]
			(*x)[IX(i, N+1)] = -(*x)[IX(i, N)]
			(*x)[IX(0, i)] = (*x)[IX(1, i)]
			(*x)[IX(N+1, i)] = (*x)[IX(N, i)]
		} else if b == 2 {
			(*x)[IX(i, 0)] = (*x)[IX(i, 1)]
			(*x)[IX(i, N+1)] = (*x)[IX(i, N)]
			(*x)[IX(0, i)] = -(*x)[IX(1, i)]
			(*x)[IX(N+1, i)] = -(*x)[IX(N, i)]
		}
	}
	(*x)[IX(0, 0)] = 0.5 * ((*x)[IX(1, 0)] + (*x)[IX(0, 1)])
	(*x)[IX(0, N+1)] = 0.5 * ((*x)[IX(1, N+1)] + (*x)[IX(0, N)])
	(*x)[IX(N+1, 0)] = 0.5 * ((*x)[IX(N, 0)] + (*x)[IX(N+1, 1)])
	(*x)[IX(N+1, N+1)] = 0.5 * ((*x)[IX(N, N+1)] + (*x)[IX(N+1, N)])
}
