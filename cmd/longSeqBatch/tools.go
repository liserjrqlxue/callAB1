package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// mergeIntervals 合并有交集的区间
func MergeIntervals(intervals [][2]int) [][2]int {
	if len(intervals) < 2 {
		return intervals
	}

	// 按起点排序
	sort.Slice(intervals, func(i, j int) bool {
		return intervals[i][0] < intervals[j][0]
	})

	merged := make([][2]int, 0)
	current := intervals[0]

	for _, interval := range intervals[1:] {
		// fmt.Printf("\t%d-%d\n", interval[0], interval[1])
		if interval[0] <= current[1] { // 有交集
			// 合并区间
			if interval[1] > current[1] {
				current[1] = interval[1]
			}
		} else {
			// 没有交集，保存当前区间并更新
			merged = append(merged, current)
			current = interval
		}
	}

	// 添加最后一个区间
	merged = append(merged, current)

	return merged
}

func MergeVariants(variants []*tracy.Variant) []*tracy.Variant {
	if len(variants) < 2 {
		return variants
	}

	// 按起点排序
	sort.Slice(variants, func(i, j int) bool {
		if variants[i].Pos == variants[j].Pos {
			if variants[i].Ref == variants[j].Ref {
				return variants[i].Alt < variants[j].Alt
			}
			return variants[i].Ref < variants[j].Ref

		}
		return variants[i].Pos < variants[j].Pos
	})

	merged := make([]*tracy.Variant, 0)
	current := variants[0]

	for _, variant := range variants[1:] {
		if variant.Pos == current.Pos && variant.Ref == current.Ref && variant.Alt == current.Alt { // 相同
			if variant.Type == "hom. ALT" && variant.Qual > current.Qual {
				current = variant
			}
		} else { // 不同
			merged = append(merged, current)
			current = variant
		}
	}

	// 添加最后一个
	merged = append(merged, current)

	return merged
}

func GetRows2MapArray(xlsx *excelize.File, sheet string) (data []map[string]string) {
	var (
		rows  = simpleUtil.HandleError(xlsx.GetRows(sheet))
		title []string
	)
	for i := range rows {
		if i == 0 {
			title = rows[i]
			continue
		}
		var item = make(map[string]string)
		for j := range rows[i] {
			item[title[j]] = rows[i][j]
		}
		data = append(data, item)
	}
	return
}

func CreateGeneInfoFromDataArray(data []map[string]string, outDir string) (geneInfo map[string]*Gene, geneList []string) {
	geneInfo = make(map[string]*Gene)

	for i := range data {
		item := data[i]

		gene := &Gene{
			ID:  item["基因名称"],
			Seq: item["目标序列"],
			// Prefix: item["测序结果名"],
			Clones: make(map[string]*Clone),
			OutDir: outDir,
		}
		gene.RefPath = filepath.Join(gene.OutDir, "ref", gene.ID+".ref.fa")

		geneInfo[gene.ID] = gene
		geneList = append(geneList, gene.ID)
	}

	return
}

func LoadSangerFromGlob(geneInfo map[string]*Gene, data []map[string]string) {
	for i := range data {
		item := data[i]

		geneID := item["基因名称"]
		gene, ok := geneInfo[geneID]
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
			cloneID := item["基因名称"] + "_C" + match[1]
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

func LoadSangerFromDataArray(geneInfo map[string]*Gene, data []map[string]string) {
	for i := range data {
		item := data[i]
		geneID := item["基因名称"]
		cloneID := item["基因名称"] + "_C" + item["克隆号"]
		sangerPath := filepath.Join(*seqDir, item["文件名"])
		sanger := &Sanger{
			Path: sangerPath,
		}

		gene, ok := geneInfo[geneID]
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

func CreateResult(geneInfo map[string]*Gene, geneList []string) {
	var resultXlsx = excelize.NewFile()

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
	simpleUtil.CheckErr(resultXlsx.SaveAs(*result), *result)
}
