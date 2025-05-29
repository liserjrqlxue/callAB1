package main

import (
	"callAB1/pkg/tracy"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// flag

var (
	bin = flag.String(
		"tracy",
		"tracy",
		"https://github.com/gear-genomics/tracy binary path",
	)
	inputDir = flag.String(
		"i",
		"",
		"input dir",
	)
	refXlsx = flag.String(
		"r",
		"目标序列.xlsx",
		"输入信息",
	)
	seqDir = flag.String(
		"s",
		"测序目录",
		"ab1 dir",
	)
	outDir = flag.String(
		"d",
		"",
		"output dir, default same as -i",
	)
	result = flag.String(
		"o",
		"result.xlsx",
		"result file name",
	)
)

func init() {
	flag.Parse()
	if *inputDir == "" || *result == "" {
		flag.Usage()
		log.Fatal("-i/-o required")
	}
	if *outDir == "" {
		*outDir = *inputDir
	}

	*refXlsx = filepath.Join(*inputDir, *refXlsx)
	*seqDir = filepath.Join(*inputDir, *seqDir)
	*result = filepath.Join(*outDir, *result)
}

type Gene struct {
	ID      string
	Seq     string
	RefPath string
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
	ID      string
	GeneID  string
	Sangers []*Sanger

	Variants     []*tracy.Variant
	MatchRegions [][2]int

	Effective bool
	Status    string
}

type Sanger struct {
	Index int
	Path  string

	Result      *tracy.Result
	MatchRegion [2]int
}

var (
	TermTrim = 20
)

func main() {
	var (
		xlsx       = simpleUtil.HandleError(excelize.OpenFile(*refXlsx))
		resultXlsx = excelize.NewFile()

		geneInfo map[string]*Gene
		geneList []string
	)
	sheet1Data := GetRows2MapArray(xlsx, "Sheet1")
	geneInfo, geneList = CreateGeneInfoFromDataArray(sheet1Data)

	sheet2Data := GetRows2MapArray(xlsx, "Sheet2")
	LoadSangerFromDataArray(geneInfo, sheet2Data)

	simpleUtil.CheckErr(os.MkdirAll(filepath.Join(*outDir, "ref"), 0755))

	for _, geneID := range geneList {
		log.Printf("loop GeneID:[%s]", geneID)
		gene := geneInfo[geneID]
		gene.RefPath = filepath.Join(*outDir, "ref", geneID+".ref.fa")
		simpleUtil.CheckErr(gene.CreateRef())

		for cloneID := range gene.Clones {
			log.Printf("loop CeneID:[%s:%s]", geneID, cloneID)
			clone := gene.Clones[cloneID]
			cloneDir := filepath.Join(*outDir, geneID, cloneID)
			simpleUtil.CheckErr(os.MkdirAll(cloneDir, 0755))

			clone.Effective = true
			for _, sanger := range clone.Sangers {
				log.Printf("loop sanger:[%s:%s:%d:%s]", geneID, cloneID, sanger.Index, sanger.Path)
				prefix := filepath.Join(cloneDir, strconv.Itoa(sanger.Index))
				result := simpleUtil.HandleError(tracy.RunSingle(*bin, gene.RefPath, sanger.Path, prefix, false))
				sanger.Result = &result
				if !result.Pass {
					clone.Effective = false
					clone.Status += result.Status
				}
				sanger.MatchRegion = result.AlignResult.MatchRegion

				log.Printf("sanger.MatchRegions %+v,clone.MatchRegions %+v", sanger.MatchRegion, clone.MatchRegions)
				clone.MatchRegions = append(clone.MatchRegions, sanger.MatchRegion)
				if result.Variants != nil {
					for _, v := range result.Variants.Variants {
						if v.Pos >= result.AlignResult.MatchRegion[0] && v.Pos <= result.AlignResult.MatchRegion[1] {
							fmt.Println(v.String())
							clone.Variants = append(clone.Variants, v)
						}
					}
				}
			}
			log.Printf("clone.MatchRegions %+v", clone.MatchRegions)
			clone.MatchRegions = MergeIntervals(clone.MatchRegions)
			log.Printf("clone.MatchRegions %+v", clone.MatchRegions)
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
	}
	resultTitle := []string{
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
		resultXlsx.SetSheetRow(
			"Sheet1",
			"A1",
			&resultTitle,
		),
	)
	for i := range geneList {
		geneID := geneList[i]
		gene := geneInfo[geneID]
		simpleUtil.CheckErr(
			resultXlsx.SetSheetRow(
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
	simpleUtil.CheckErr(resultXlsx.SaveAs(*result))

}
