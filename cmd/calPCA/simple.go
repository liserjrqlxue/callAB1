package main

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"

	"github.com/liserjrqlxue/DNA/pkg/util"
	"github.com/liserjrqlxue/PrimerDesigner/v2/pkg/batch"
	"github.com/liserjrqlxue/PrimerDesigner/v2/pkg/cy0130"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// 化学补充
func runFix(prefix, template string, segmentList []string, segmentMap map[string]*Seq) {
	prefix = simpleUtil.HandleError(filepath.Abs(prefix))
	var (
		listPath  = prefix + "-化补2清单.xlsx"
		orderPath = prefix + "-化补2引物订购单.xlsx"

		listXlsx    = excelize.NewFile()
		listSheet   = "化补2清单"
		plateSheet  = "化补2引物板位图"
		panelHeight = 13
		bgStyleMap  = batch.CreateStyles(listXlsx)

		orderXlsx       = simpleUtil.HandleError(excelize.OpenFile(template))
		orderSheet      = "引物订购单"
		orderColOffset  = 4
		orderRowOffset  = 17
		orderPanelHight = 96
		orderRowSkip    = 8

		chemicalSuppList []*util.Seq
		index            = -1
	)

	// 化补2清单
	simpleUtil.HandleError(listXlsx.NewSheet(listSheet))
	listXlsx.SetSheetRow(listSheet, "A1", &[]string{"片段名称", "序列", "长度", "设计成功"})
	listXlsx.SetColWidth(listSheet, "B", "B", 100)
	// "化补2引物板位图"
	simpleUtil.HandleError(listXlsx.NewSheet(plateSheet))
	listXlsx.SetColWidth(plateSheet, "A", "N", 10)

	for _, id := range segmentList {
		segment, ok := segmentMap[id]
		if ok && segment.CloneHit == 0 {
			id := segment.ID + "C"
			seq := util.NewSeq(id, segment.Seq, filepath.Join(*outputDir, segment.ID, id), true, false)
			chemicalSuppList = append(chemicalSuppList, seq)
		}
	}

	for i, seq := range chemicalSuppList {
		pass := runSimple(seq)

		cell1 := fmt.Sprintf("A%d", i+2)
		cell2 := fmt.Sprintf("D%d", i+2)
		listXlsx.SetSheetRow(listSheet, cell1, &[]any{seq.Name, seq.Seq, len(seq.Seq), pass})

		if !pass {
			continue
		}
		index++
		panelIndex := index / 6
		rowOffset := panelIndex * panelHeight
		panelCol := index % 6 * 2
		if panelCol == 0 {
			panelInit(listXlsx, plateSheet, rowOffset, bgStyleMap[-1])
		}

		simpleUtil.CheckErr(listXlsx.SetCellStyle(listSheet, cell1, cell2, bgStyleMap[index%3]))
		halfBreak := (seq.PrimerPairCount + 1) / 2
		for j, pair := range seq.SegmentPairs {
			row := j/halfBreak*4 + j%halfBreak

			cellName1 := CoordinatesToCellName(3+panelCol, row+6+rowOffset)
			cellName2 := CoordinatesToCellName(3+panelCol+1, row+6+rowOffset)
			listXlsx.SetSheetRow(plateSheet, cellName1, &[]string{pair.Left.Name, pair.Right.Name})
			simpleUtil.CheckErr(listXlsx.SetCellStyle(plateSheet, cellName1, cellName2, bgStyleMap[index%3]))

			batch.WritePrimerOrder(
				orderXlsx, orderSheet,
				orderColOffset, orderRowOffset+row+panelCol*orderRowSkip+panelIndex*orderPanelHight,
				1, pair.Left, false,
			)
			batch.WritePrimerOrder(
				orderXlsx, orderSheet,
				orderColOffset, orderRowOffset+row+(panelCol+1)*orderRowSkip+panelIndex*orderPanelHight,
				1, pair.Right, true,
			)

		}
	}

	simpleUtil.CheckErr(listXlsx.DeleteSheet("Sheet1"))
	log.Printf("SaveAs(%s)", listPath)
	simpleUtil.CheckErr(listXlsx.SaveAs(listPath))

	log.Printf("SaveAs(%s)", orderPath)
	simpleUtil.CheckErr(orderXlsx.SaveAs(orderPath))

	// zip
	var baseName = filepath.Base(prefix)
	var patterns = []string{
		baseName + ".Sanger结果.xlsx",
		baseName + "-化补2清单.xlsx",
		baseName + "-化补2引物订购单.xlsx",
	}
	simpleUtil.CheckErr(cy0130.CompressArchive(filepath.Dir(prefix), prefix+".calPCA.zip", patterns))

}

func runSimple(seq *util.Seq) (status bool) {
	defer func() {
		if e := recover(); e != nil {
			status = false
		}
	}()
	// slog.SetLogLoggerLevel(slog.LevelDebug)
	// log.Printf("%+v\n", seq)
	// logFile := osUtil.Create("PrimerDesigner.log")
	// slog.SetDefault(slog.New(slog.NewTextHandler(logFile, nil)))
	simpleUtil.CheckErr(os.MkdirAll(filepath.Dir(seq.OutputPrefix), 0755))

	// log
	// util.InitLog(seq)

	// slog.Info("Step1 Done")
	// seq.Step1()
	seq.CalAll()

	if !seq.FindPrimerPairs() || seq.PrimerPairCount > 16 {
		slog.Error("SegmentSplit fail", "seq.Name", seq.Name, "seq.PrimerPairCount", seq.PrimerPairCount)
		return false
	}
	// output after design
	// slog.Info("Step3 Start")
	// seq.Step3()
	return true
}

func panelInit(xlsx *excelize.File, sheet string, rowOffset, styleID int) {
	var ()
	// "化补2引物板位图"
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetCol(
		sheet,
		fmt.Sprintf("A%d", 1+rowOffset),
		&[]string{"", "金唯智引物订购订单号：", "", "金唯智引物订购单板位图："},
	)
	xlsx.SetSheetCol(
		sheet,
		fmt.Sprintf("B%d", 6+rowOffset),
		&[]string{"A", "B", "C", "D", "E", "F", "G", "H"},
	)
	xlsx.SetSheetRow(
		sheet,
		fmt.Sprintf("C%d", 5+rowOffset),
		&[]int{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
	)

	simpleUtil.CheckErr(xlsx.SetCellStyle(sheet, fmt.Sprintf("B%d", 5+rowOffset), fmt.Sprintf("N%d", 13+rowOffset), styleID))
}
