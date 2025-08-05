package main

import (
	"strings"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// 添加 拼接引物板
func addSplicedPrimerPlate(xlsx *excelize.File, path string, primerACC map[string][2]float64) {
	var (
		sheet = "拼接引物板"
	)
	inputOderXlsx := simpleUtil.HandleError(excelize.OpenFile(path))
	rows := simpleUtil.HandleError(inputOderXlsx.GetRows(sheet))
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	for i, row := range rows {
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
					CoordinatesToCellName(j+1, i+1+len(rows)+1),
					accs[0],
				)
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(rows)*2+2),
					accs[1],
				)
			} else {
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(rows)+1),
					cell,
				)
				xlsx.SetCellValue(
					sheet,
					CoordinatesToCellName(j+1, i+1+len(rows)*2+2),
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
			CoordinatesToCellName(1, 1+len(rows)+1),
			"平均准确率(%)",
		)
		xlsx.SetCellValue(
			sheet,
			CoordinatesToCellName(1, 1+len(rows)*2+2),
			"参考单步准确率(%)",
		)
	}
}
