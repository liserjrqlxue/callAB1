package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"os"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
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

func (gene *Gene) CloneRun(cloneID string) {
	log.Printf("loop CeneID:[%s:%s]", gene.ID, cloneID)
	clone := gene.Clones[cloneID]
	simpleUtil.CheckErr(os.MkdirAll(clone.CloneDir, 0755))

	clone.Effective = true
	for _, sanger := range clone.Sangers {
		log.Printf("loop sanger:[%s:%s:%d:%s]", gene.ID, clone.ID, sanger.Index, sanger.Path)
		SangerRun(sanger, clone, gene.RefPath)
	}

	clone.MatchRegions = MergeIntervals(clone.MatchRegions)
	clone.Variants = MergeVariants(clone.Variants)
	if clone.Effective {
		if len(clone.MatchRegions) > 1 {
			clone.Effective = false
			clone.Status += "CoverGap"
		} else if clone.MatchRegions[0][0] > TermTrim+1 || clone.MatchRegions[0][1] < len(gene.Seq)-TermTrim {
			clone.Effective = false
			clone.Status += "CoverFail"
			log.Printf("CoverFail [%d,%d] vs %d[%d,%d]", clone.MatchRegions[0][0], clone.MatchRegions[0][1], len(gene.Seq), TermTrim+1, len(gene.Seq)-TermTrim)
		}
	}

	if clone.Effective {
		gene.EffectiveClones = append(gene.EffectiveClones, clone.ID)
		if len(clone.Variants) > 0 {
			gene.MismatchClones++
			if len(clone.Variants) == 1 {
				gene.Mismatch1Clones = append(gene.Mismatch1Clones, clone.ID)
			}
		} else {
			gene.CorrectClones = append(gene.CorrectClones, clone.ID)
		}
	} else {
		gene.FailClones++
	}
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
