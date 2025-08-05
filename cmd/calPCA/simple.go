package main

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"

	"github.com/liserjrqlxue/DNA/pkg/util"
	"github.com/liserjrqlxue/PrimerDesigner/v2/pkg/batch"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// 化学补充
func runFix(prefix, template string, segmentList []string, segmentMap map[string]*Seq) {
	var (
		listPath  = prefix + "化补2清单.xlsx"
		orderPath = prefix + "化补2引物订购单.xlsx"

		xlsx   = excelize.NewFile()
		sheet1 = "化补2清单"
		sheet2 = "化补2引物板位图"

		orderXlsx       = simpleUtil.HandleError(excelize.OpenFile(template))
		orderSheet      = "引物订购单"
		orderColOffset  = 4
		orderRowOffset  = 17
		orderPanelHight = 96
		orderRowSkip    = 8

		chemicalSuppList []*util.Seq

		panelHeight = 13
		index       = -1
	)

	// 化补2清单
	simpleUtil.HandleError(xlsx.NewSheet(sheet1))
	xlsx.SetSheetRow(sheet1, "A1", &[]string{"片段名称", "序列", "长度", "设计成功"})
	// "化补2引物板位图"
	simpleUtil.HandleError(xlsx.NewSheet(sheet2))

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

		xlsx.SetSheetRow(sheet1, fmt.Sprintf("A%d", i+2), &[]any{seq.Name, seq.Seq, len(seq.Seq), pass})

		if !pass {
			continue
		}
		index++
		panelIndex := index / 6
		rowOffset := panelIndex * panelHeight
		panelInit(xlsx, sheet2, rowOffset)

		// fmt.Printf("%s\t%s\t%d\n", seq.Name, string(rune('A'+i*2)), seq.PrimerPairCount)
		halfBreak := (seq.PrimerPairCount + 1) / 2
		for j, pair := range seq.SegmentPairs {
			// row := string(rune('A' + j/halfBreak*4 + j%halfBreak))
			// fmt.Printf(
			// 	"\t%d\t%s%d:%s\t%s%d:%s\n",
			// 	j+1,
			// 	row, index*2+1,
			// 	pair.Left.Name,
			// 	row, index*2+2,
			// 	pair.Right.Name,
			// )
			row := j/halfBreak*4 + j%halfBreak
			cellName := simpleUtil.HandleError(excelize.CoordinatesToCellName(3+index*2, row+6+rowOffset))
			xlsx.SetSheetRow(sheet2, cellName, &[]string{pair.Left.Name, pair.Right.Name})

			batch.WritePrimerOrder(
				orderXlsx, orderSheet,
				orderColOffset, orderRowOffset+row+(index*2+1)*orderRowSkip+panelIndex*orderPanelHight,
				1, pair.Left, false,
			)
			batch.WritePrimerOrder(
				orderXlsx, orderSheet,
				orderColOffset, orderRowOffset+row+index*2*orderRowSkip+panelIndex*orderPanelHight,
				1, pair.Right, true,
			)

		}
	}

	simpleUtil.CheckErr(xlsx.DeleteSheet("Sheet1"))
	log.Printf("SaveAs(%s)", listPath)
	simpleUtil.CheckErr(xlsx.SaveAs(listPath))

	log.Printf("SaveAs(%s)", orderPath)
	simpleUtil.CheckErr(orderXlsx.SaveAs(orderPath))
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

func panelInit(xlsx *excelize.File, sheet string, rowOffset int) {
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
}
