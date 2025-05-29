package main

import (
	"callAB1/pkg/tracy"
	"log"
	"path/filepath"
	"sort"

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

func CreateGeneInfoFromDataArray(data []map[string]string) (geneInfo map[string]*Gene, geneList []string) {
	geneInfo = make(map[string]*Gene)

	for i := range data {
		item := data[i]
		gene := &Gene{
			ID:     item["基因名称"],
			Seq:    item["目标序列"],
			Prefix: item["测序结果名"],
			Clones: make(map[string]*Clone),
		}
		geneInfo[gene.ID] = gene
		geneList = append(geneList, gene.ID)
	}

	return
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
				ID:     cloneID,
				GeneID: geneID,
			}
			gene.Clones[cloneID] = clone
		}
		clone.Sangers = append(clone.Sangers, sanger)
		sanger.Index = len(clone.Sangers)
	}

}
