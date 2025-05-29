package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"os"
)

type Gene struct {
	ID      string
	Seq     string
	RefPath string
	OutDir  string
	Prefix  string

	Clones map[string]*Clone

	FailClones      int
	MismatchClones  int
	EffectiveClones []string
	Mismatch1Clones []string
	CorrectClones   []string
}

func (g *Gene) CreateRef() error {
	return os.WriteFile(
		g.RefPath,
		fmt.Appendf(nil, ">%s\n%s\n", g.ID, g.Seq),
		0644,
	)
}

type Clone struct {
	ID       string
	GeneID   string
	CloneDir string
	Sangers  []*Sanger

	Variants     []*tracy.Variant
	MatchRegions [][2]int

	Effective bool
	Status    string
}

type Sanger struct {
	Index  int
	Path   string
	Prefix string

	Result      *tracy.Result
	MatchRegion [2]int
}
