package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"sync"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

type Genes struct {
	GeneInfo map[string]*Gene
	GeneList []string
	OutDir   string
}

func (genes *Genes) LoadInput(excel string) {
	var xlsx = simpleUtil.HandleError(excelize.OpenFile(excel))

	var sheet1Data = GetRows2MapArray(xlsx, "Sheet1")
	genes.GetGenes(sheet1Data)

	var sheet2Data = GetRows2MapArray(xlsx, "Sheet2")
	genes.LoadSangerFromGlob(sheet2Data)

}

func (genes *Genes) GetGenes(data []map[string]string) {
	for i := range data {
		item := data[i]

		gene := &Gene{
			ID:  item["基因名称"],
			Seq: item["目标序列"],
			// Prefix: item["测序结果名"],
			Clones: make(map[string]*Clone),
			OutDir: genes.OutDir,
		}
		gene.RefPath = filepath.Join(gene.OutDir, "ref", gene.ID+".ref.fa")

		genes.GeneInfo[gene.ID] = gene
		genes.GeneList = append(genes.GeneList, gene.ID)
	}
}

func (genes *Genes) LoadSangerFromDataArray(data []map[string]string) {
	for i := range data {
		item := data[i]
		geneID := item["基因名称"]
		cloneID := item["基因名称"] + "_C" + item["克隆号"]
		sangerPath := filepath.Join(*seqDir, item["文件名"])
		sanger := &Sanger{
			Path: sangerPath,
		}

		gene, ok := genes.GeneInfo[geneID]
		if !ok {
			log.Fatalf("GeneID not exists:[%s]", geneID)
		}
		clone, ok := gene.Clones[cloneID]
		if !ok {
			clone = &Clone{
				ID:         cloneID,
				GeneID:     geneID,
				RefPath:    gene.RefPath,
				GeneLength: len(gene.Seq),
				CloneDir:   filepath.Join(gene.OutDir, geneID, cloneID),
			}
			gene.Clones[cloneID] = clone
		}
		clone.Sangers = append(clone.Sangers, sanger)

		sanger.Index = len(clone.Sangers)
		sanger.Prefix = filepath.Join(clone.CloneDir, strconv.Itoa(sanger.Index))
	}
}

func (genes *Genes) LoadSangerFromGlob(data []map[string]string) {
	for i := range data {
		item := data[i]

		geneID := item["基因名称"]
		gene, ok := genes.GeneInfo[geneID]
		if !ok {
			log.Fatalf("GeneID not exists:[%s]", geneID)
		}

		prefix := item["测序结果名"]
		reg := regexp.MustCompile(prefix + `-(\d+)`)
		pattern := filepath.Join(*seqDir, prefix+"*.ab1")
		sangerFiles, err := filepath.Glob(pattern)
		if err != nil {
			log.Fatalf("can not glob ab1 file:[%s][%v]", pattern, err)
		}

		for _, file := range sangerFiles {
			sanger := &Sanger{
				Path: file,
			}
			baseName := filepath.Base(file)
			match := reg.FindStringSubmatch(baseName)
			if len(match) < 2 {
				log.Fatalf("can not parse cloneID:[reg:%s][name:%s][match:%+v][%s]", reg, baseName, match, file)
			}
			cloneID := prefix + "-" + match[1]
			clone, ok := gene.Clones[cloneID]
			if !ok {
				clone = &Clone{
					ID:         cloneID,
					GeneID:     geneID,
					RefPath:    gene.RefPath,
					GeneLength: len(gene.Seq),
					CloneDir:   filepath.Join(gene.OutDir, geneID, cloneID),
				}
				gene.Clones[cloneID] = clone
			}
			clone.Sangers = append(clone.Sangers, sanger)

			sanger.Index = len(clone.Sangers)
			sanger.Prefix = filepath.Join(clone.CloneDir, strconv.Itoa(sanger.Index))
		}
	}
}

func (genes *Genes) GeneRun() {
	var wg sync.WaitGroup
	simpleUtil.CheckErr(os.MkdirAll(filepath.Join(genes.OutDir, "ref"), 0755))

	for _, geneID := range genes.GeneList {
		log.Printf("loop GeneID:[%s]", geneID)
		gene, ok := genes.GeneInfo[geneID]
		if !ok {
			log.Fatalf("geneID[%s] not exists", geneID)
		}
		wg.Add(1)
		go func() {
			defer wg.Done()
			gene.Run()
		}()
	}

	wg.Wait()
}

func (gene *Genes) CreateResult() {
	var xlsx = excelize.NewFile()

	title := []string{
		"基因名称",
		"目标序列",
		"长度",
		"测序有效克隆数",
		"测序失败克隆数",
		"测序无法匹配克隆数", "有效克隆号",
		"正确克隆号",
		"1处错误克隆号",
	}
	simpleUtil.CheckErr(
		xlsx.SetSheetRow(
			"Sheet1",
			"A1",
			&title,
		),
	)

	for i := range gene.GeneList {
		geneID := gene.GeneList[i]
		gene := gene.GeneInfo[geneID]
		simpleUtil.CheckErr(
			xlsx.SetSheetRow(
				"Sheet1",
				"A"+strconv.Itoa(i+2),
				&[]any{
					gene.ID,
					gene.Seq,
					len(gene.Seq),
					len(gene.EffectiveClones),
					gene.FailClones,
					gene.MismatchClones,
					strings.Join(gene.EffectiveClones, " "),
					strings.Join(gene.CorrectClones, " "),
					strings.Join(gene.Mismatch1Clones, " "),
				},
			),
		)

		for cloneID := range gene.Clones {
			clone := gene.Clones[cloneID]
			fmt.Printf("%s\t%s\t%t\t%s\n", geneID, cloneID, clone.Effective, clone.Status)
		}
	}
	simpleUtil.CheckErr(xlsx.SaveAs(*result), *result)
}

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

func (gene *Gene) Run() {
	simpleUtil.CheckErr(gene.CreateRef())

	for cloneID := range gene.Clones {
		gene.CloneRun(cloneID)
	}
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

	GeneLength int

	Sangers []*Sanger

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
		} else if clone.MatchRegions[0][0] > TermTrim+1 || clone.MatchRegions[0][1] < clone.GeneLength-TermTrim {
			clone.Effective = false
			clone.Status += "CoverFail"
			log.Printf("CoverFail [%d,%d] vs %d[%d,%d]", clone.MatchRegions[0][0], clone.MatchRegions[0][1], clone.GeneLength, TermTrim+1, clone.GeneLength-TermTrim)
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
