package main

import (
	"log"

	"github.com/hajimehoshi/ebiten/v2"
)

type Game struct {
	pixels []byte
	fluid  Fluid
}

type Fluid struct {
	Size int
	N    int
	Dt   float64
	Diff float64
	Visc float64

	Sources []float64
	Density []float64

	U []float64
	V []float64

	USource []float64
	VSource []float64
}

func NewFluid(size int, diff, visc float64) Fluid {
	var fluid Fluid = Fluid{
		Size:    size,
		Dt:      float64(1. / 60.),
		Diff:    diff,
		Visc:    visc,
		Sources: make([]float64, (size+2)*(size+2)),
		Density: make([]float64, (size+2)*(size+2)),
		U:       make([]float64, (size+2)*(size+2)),
		V:       make([]float64, (size+2)*(size+2)),
		USource: make([]float64, (size+2)*(size+2)),
		VSource: make([]float64, (size+2)*(size+2)),
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
	var x []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	copy(x, f.Density)
	a := float64(f.Dt * f.Diff * N * N)
	for k := 0; k < 20; k++ {
		for i := 1; i <= N; i++ {
			for j := 1; j <= N; j++ {
				x[IX(i, j)] = (f.Density[IX(i, j)] + a*(x[IX(i+1, j)]+x[IX(i-1, j)]+x[IX(i, j+1)]+x[IX(i, j-1)])) / (float64(1) + 4*a)
			}
		}
	}
	set_bnd(N, 0, &x)
	copy(f.Density, x)
}

func (f *Fluid) advect() {
	var i0, j0, i1, j1 int
	var x, y, s0, t0, s1, t1, dt0 float64
	var new_density []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	copy(new_density, f.Density)
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
			new_density[IX(i, j)] = s0*(t0*f.Density[IX(i0, j0)]+t1*f.Density[IX(i0, j1)]) + s1*(t0*f.Density[IX(i1, j0)]+t1*f.Density[IX(i1, j1)])
		}
	}
	copy(f.Density, new_density)
}

func (f *Fluid) advect_velocity() {
	var i0, j0, i1, j1 int
	var x, y, s0, t0, s1, t1, dt0 float64
	var new_U []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	var new_V []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	copy(new_U, f.U)
	copy(new_V, f.V)
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
			new_U[IX(i, j)] = s0*(t0*f.U[IX(i0, j0)]+t1*f.U[IX(i0, j1)]) + s1*(t0*f.U[IX(i1, j0)]+t1*f.U[IX(i1, j1)])
			new_V[IX(i, j)] = s0*(t0*f.V[IX(i0, j0)]+t1*f.V[IX(i0, j1)]) + s1*(t0*f.V[IX(i1, j0)]+t1*f.V[IX(i1, j1)])
		}
	}
	copy(f.U, new_U)
	copy(f.V, new_V)
}

func (f *Fluid) project() {
	var i, j, k int
	var h float64
	var div, p []float64
	div = make([]float64, (f.Size*2)*(f.Size*2))
	p = make([]float64, (f.Size*2)*(f.Size*2))
	h = 1.0 / N
	for i = 1; i <= N; i++ {
		for j = 1; j <= N; j++ {
			div[IX(i, j)] = -0.5 * h * (f.V[IX(i+1, j)] - f.V[IX(i-1, j)] + f.U[IX(i, j+1)] - f.U[IX(i, j-1)])
			p[IX(i, j)] = 0
		}
	}
	set_bnd(N, 0, &div)
	set_bnd(N, 0, &p)
	for k = 0; k < 20; k++ {
		for i = 1; i <= N; i++ {
			for j = 1; j <= N; j++ {
				p[IX(i, j)] = (div[IX(i, j)] + p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)]) / 4
			}
		}
		set_bnd(N, 0, &p)
	}
	for i = 1; i <= N; i++ {
		for j = 1; j <= N; j++ {
			f.V[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]) / h
			f.U[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]) / h
		}
	}
	set_bnd(N, 1, &f.U)
	set_bnd(N, 2, &f.V)
}

func (f *Fluid) diffuse_velocity() {
	if f.Visc == 0 {
		return
	}
	var x []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	var y []float64 = make([]float64, (f.Size+2)*(f.Size+2))
	copy(x, f.U)
	copy(y, f.V)
	a := float64(f.Dt * f.Visc * N * N)
	for k := 0; k < 20; k++ {
		for i := 1; i <= N; i++ {
			for j := 1; j <= N; j++ {
				x[IX(i, j)] = (f.U[IX(i, j)] + a*(x[IX(i+1, j)]+x[IX(i-1, j)]+x[IX(i, j+1)]+x[IX(i, j-1)])) / (float64(1) + 4*a)
				y[IX(i, j)] = (f.V[IX(i, j)] + a*(y[IX(i+1, j)]+y[IX(i-1, j)]+y[IX(i, j+1)]+y[IX(i, j-1)])) / (float64(1) + 4*a)
			}
		}
	}
	set_bnd(N, 0, &x)
	copy(f.U, x)
	copy(f.V, y)
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

const (
	N = 256
)

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

func (g *Game) Update() error {
	g.fluid.add_source()
	g.fluid.diffuse()
	g.fluid.advect()
	g.fluid.diffuse_velocity()
	g.fluid.advect_velocity()
	g.fluid.project()
	return nil
}

func main() {
	fluid := NewFluid(N, 0, 0)

	for i := N/2 - 4; i <= N/2+4; i++ {
		for j := 1; j <= 3; j++ {
			fluid.USource[IX(i, j)] = 10
		}
	}
	for i := N/2 - 4; i <= N/2+4; i++ {
		for j := N - 2; j <= N; j++ {
			fluid.USource[IX(i, j)] = -10
		}
	}
	fluid.Sources[IX(N/2, 1)] = 100
	fluid.Sources[IX(N/2, N)] = 100

	g := &Game{
		fluid: fluid,
	}

	ebiten.SetWindowSize(900, 900)
	ebiten.SetWindowTitle("Testing")
	if err := ebiten.RunGame(g); err != nil {
		log.Fatal(err)
	}
}
