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

	// stats
	FailClones      int
	MismatchClones  int
	EffectiveClones []string
	Mismatch1Clones []string
	CorrectClones   []string
}

func (gene *Gene) CreateRef() error {
	return os.WriteFile(
		gene.RefPath,
		fmt.Appendf(nil, ">%s\n%s\n", gene.ID, gene.Seq),
		0644,
	)
}

func (gene *Gene) CloneRun(cloneID string) {
	log.Printf("loop CeneID:[%s:%s]", gene.ID, cloneID)

	var clone, ok = gene.Clones[cloneID]
	if !ok {
		log.Fatalf("cloneID[%s] not exists for Gene[%s]", cloneID, gene.ID)
	}

	clone.Run()

	gene.UpdateStats(clone)
}

func (gene *Gene) UpdateStats(clone *Clone) {
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
	RefPath  string
	CloneDir string
	Sangers  []*Sanger

	Variants     []*tracy.Variant
	MatchRegions [][2]int

	Effective bool
	Status    string
}

func (clone *Clone) Run() {
	simpleUtil.CheckErr(os.MkdirAll(clone.CloneDir, 0755))

	clone.Effective = true
	clone.SangersRun()

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
}

func (clone *Clone) SangersRun() {
	for _, sanger := range clone.Sangers {
		log.Printf("loop sanger:[%s:%d:%s]", clone.ID, sanger.Index, sanger.Path)
		clone.SangerRun(sanger)
	}
}

func (clone *Clone) SangerRun(sanger *Sanger) {
	sanger.Run(clone.RefPath)

	if !sanger.Result.Pass {
		clone.Effective = false
		clone.Status += sanger.Result.Status
	}

	log.Printf("sanger.MatchRegions %+v,clone.MatchRegions %+v", sanger.MatchRegion, clone.MatchRegions)
	clone.MatchRegions = append(clone.MatchRegions, sanger.MatchRegion)
	if sanger.Result.Variants != nil {
		for _, v := range sanger.Result.Variants.Variants {
			if v.Pos >= sanger.Result.AlignResult.MatchRegion[0] && v.Pos <= sanger.Result.AlignResult.MatchRegion[1] {
				fmt.Println(v.String())
				clone.Variants = append(clone.Variants, v)
			}
		}
	}
}

type Sanger struct {
	Index  int
	Path   string
	Prefix string

	Result      *tracy.Result
	MatchRegion [2]int
}

func (sanger *Sanger) Run(ref string) {
	result := simpleUtil.HandleError(tracy.RunSingle(*bin, ref, sanger.Path, sanger.Prefix, false))
	sanger.Result = &result
	sanger.MatchRegion = sanger.Result.AlignResult.MatchRegion
}
