package main

import (
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
