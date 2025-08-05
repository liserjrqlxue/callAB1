package main

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/samber/lo"
	"github.com/xuri/excelize/v2"
)

// 添加 拼接引物板
func addSplicedPrimerPlate(xlsx *excelize.File, sheet, path string, primerACC map[string][2]float64) {
	var (
		orderXlsx = simpleUtil.HandleError(excelize.OpenFile(path))
		orderRows = simpleUtil.HandleError(orderXlsx.GetRows(sheet))
	)

	simpleUtil.HandleError(xlsx.NewSheet(sheet))

	for i, row := range orderRows {
		for j, cell := range row {
			xlsx.SetCellValue(
				sheet,
				CoordinatesToCellName(j+1, i+1),
				cell,
			)
			pID := strings.Split(cell, "\n")[0]
			accs, ok := primerACC[pID]
			if ok {
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(orderRows)+1),
					accs[0],
				)
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(orderRows)*2+2),
					accs[1],
				)
			} else {
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(orderRows)+1),
					cell,
				)
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(orderRows)*2+2),
					cell,
				)
			}
		}
		xlsx.SetCellValue(
			sheet,
			CoordinatesToCellName(1, 1),
			"板位",
		)
		xlsx.SetCellValue(
			sheet,
			CoordinatesToCellName(1, 1+len(orderRows)+1),
			"平均准确率(%)",
		)
		xlsx.SetCellValue(
			sheet,
			CoordinatesToCellName(1, 1+len(orderRows)*2+2),
			"参考单步准确率(%)",
		)
	}
}

type Verification struct {
	SnvRatio float64
	InsRatio float64
	DelRatio float64
	Status   string
}

func (vf *Verification) BatchValidation(snvRatios, insRatios, delRatios []float64) {
	vf.SnvRatio = lo.Mean(snvRatios)
	vf.InsRatio = lo.Mean(insRatios)
	vf.DelRatio = lo.Mean(delRatios)

	if vf.SnvRatio >= SnvRatio || vf.InsRatio >= InsRatio || vf.DelRatio >= DelRatio {
		vf.Status = "不合格"
	}
}

// 添加 批次统计
func addBatchStats(xlsx *excelize.File, sheet string, batchVerif *Verification) {
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A1", &BatchTitle)

	xlsx.SetSheetRow(
		sheet, "A2",
		&[]any{
			batchVerif.SnvRatio,
			batchVerif.InsRatio,
			batchVerif.DelRatio,
			batchVerif.Status,
		},
	)

}

// 添加 测序结果板位图
func AddSequencingResultPlate(xlsx *excelize.File, sheet string, geneList []string, geneMap map[string]*Seq) {
	var (
		bgColor1 = "#8CDDFA"
		bgColor2 = "#FFC000"
		bgColor3 = "#FFD966"

		bgStyle1 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor1},
				Pattern: 1,
			},
		}
		bgStyle2 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor2},
				Pattern: 1,
			},
		}
		bgStyle3 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor3},
				Pattern: 1,
			},
		}

		style1 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle1))
		style2 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle2))
		style3 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle3))
	)

	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A3", &[]string{"基因名称", "长度", "节数"})
	xlsx.SetSheetRow(sheet, "E3", &PlateCols)
	xlsx.SetSheetRow(sheet, "S3", &PlateCols)
	xlsx.SetCellStr(sheet, "D2", "批次号")
	xlsx.SetSheetCol(sheet, "D4", &PlateRows)
	xlsx.SetCellStr(sheet, "R2", "批次号")
	xlsx.SetSheetCol(sheet, "R4", &PlateRows)
	xlsx.SetCellStyle(sheet, "B4", "D11", style1)
	xlsx.SetCellStyle(sheet, "R4", "R11", style1)
	xlsx.SetCellStyle(sheet, "E3", "P3", style2)
	xlsx.SetCellStyle(sheet, "S3", "AD3", style2)

	for i, geneID := range geneList {
		gene := geneMap[geneID]
		xlsx.SetCellStr(sheet, fmt.Sprintf("A%d", 4+i), geneID)
		xlsx.SetCellStr(sheet, fmt.Sprintf("B%d", 4+i), fmt.Sprintf("%dbp", len(gene.Seq)))
		xlsx.SetCellInt(sheet, fmt.Sprintf("C%d", 4+i), len(gene.SubSeq))

		for j, segment := range gene.SubSeq {
			cellName := CoordinatesToCellName(19+j, 4+i)
			xlsx.SetCellInt(sheet, cellName, segment.End-segment.Start)
			cellName = CoordinatesToCellName(5+j, 4+i)
			xlsx.SetCellStr(sheet, cellName, "N")
			xlsx.SetCellStyle(sheet, cellName, cellName, style3)
			if len(segment.rIDs) > 0 {
				xlsx.SetCellStr(sheet, cellName, segment.ID+"-"+segment.rIDs[0])
				xlsx.SetCellStr(sheet, cellName, segment.ID+"-"+strings.Join(segment.rIDs, "、"))
			}
		}
	}
}

// 添加 片段结果
func addSegmentResult(xlsx *excelize.File, sheet string, segmentList []string, segmentLines map[string][]any) {
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A1", &SeqTitle)

	var row = 2
	for _, id := range segmentList {
		line := segmentLines[id]
		xlsx.SetSheetRow(sheet, fmt.Sprintf("A%d", row), &line)
		row++
	}
}

// 添加 Sanger结果
func addSangerResult(xlsx *excelize.File, sheet string, segmentList []string, resultLines map[string][][]any) (primerACC map[string][2]float64) {
	primerACC = make(map[string][2]float64)

	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A1", &ResultTitle)

	var row = 2
	for _, id := range segmentList {
		for _, line := range resultLines[id] {
			pID := line[1].(string)
			valid := line[7].(int)
			acc := line[19].(float64)
			geomMeanAcc := line[21].(float64)

			if valid > 0 {
				primerACC[pID] = [2]float64{acc, geomMeanAcc}
			}
			xlsx.SetSheetRow(sheet, fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	return
}

func WriteSliceSheet(xlsx *excelize.File, sheet string, list, title []string, data map[string][][]any) {
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A1", &title)

	row := 2
	for _, v := range list {
		for _, s := range data[v] {
			xlsx.SetSheetRow(sheet, "A"+strconv.Itoa(row), &s)
			row++
		}
	}
}
